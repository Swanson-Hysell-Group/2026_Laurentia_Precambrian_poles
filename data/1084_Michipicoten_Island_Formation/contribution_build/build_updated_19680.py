"""
Produce the merged, updated Fairchild et al. (2017) MagIC contribution (19680)
carrying BOTH of its poles.

The published contribution earthref.org/MagIC/19680 archives the sites, samples,
specimens, measurements, and ages of Fairchild et al. (2017) for three late-stage
Midcontinent Rift successions but has NO location-level pole result for either of
its paleomagnetic poles. This script downloads 19680 and adds, without altering
anything else, the two location-level poles of the study together with the small
number of literature site means each pole requires that are not already in 19680:

  Michipicoten Island Formation pole  (on the 'Michipicoten Island' location)
    - adds the two Palmer & Davis (1987) cooling-unit-mean SITES
      (P_QuebecHarbour, P_Davieaux), result_type='a'
    - sets the pole: 17.0 degN, 174.7 degE, A95 4.4, k 47.9, N 23
    - the 21 South Shore sites in the pole are already in 19680
    (built by MIF_build_contribution.py -> ../sites.txt, ../locations.txt)

  Schroeder-Lutsen basalts pole       (on the 'Two Island River' location)
    - adds the 10 Tauxe & Kodama (2009) flow SITES (ns006-ns015), result_type='i'
    - sets the pole: 27.1 degN, 187.8 degE, A95 3.0, k 45.4, N 50
    - the 40 Two Island River 'mag' tilt-corrected sites in the pole are already
      in 19680
    (built by SLB_build_contribution.py ->
     ../../1090_Schroeder_Lutsen_Basalts/{sites,locations}.txt)

All other tables of 19680 -- including the full measurements table -- pass
through unchanged. Output: a single merged MagIC upload file written into the
1084 pole folder, ready to upload to MagIC as an update of contribution 19680
that contains both Fairchild et al. (2017) poles.

Note: this downloads ~2.4 MB (incl. ~12,000 measurement rows) from MagIC; run
with a network connection. Do not commit the large merged measurements file.
"""

import shutil
import tempfile
from pathlib import Path

import pandas as pd
import pmagpy.ipmag as ipmag

MAGIC_ID = '19680'
HERE = Path(__file__).parent
MIF_DIR = HERE.parent                                  # data/1084_Michipicoten_Island_Formation
SLB_DIR = MIF_DIR.parent / '1090_Schroeder_Lutsen_Basalts'

# 19680 location names that carry each pole
MIF_LOCATION = 'Michipicoten Island'
SLB_LOCATION = 'Two Island River'

CIT_TAUXE = '10.1016/j.pepi.2009.07.006'               # Tauxe & Kodama (2009) sites

POLE_COLS = ['result_type', 'result_name', 'method_codes', 'citations',
             'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95',
             'pole_k', 'pole_n_sites', 'sites', 'description']


def read_magic(path):
    return pd.read_csv(path, sep='\t', skiprows=1, dtype=str).fillna('')


def write_magic(path, table_type, df):
    with open(path, 'w') as f:
        f.write(f'tab delimited\t{table_type}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')


def set_pole(loc, location_name, pole_row):
    """Write the pole_* result onto the matching locations row (in place)."""
    for col in POLE_COLS:
        if col not in loc.columns:
            loc[col] = ''
    mask = loc['location'] == location_name
    if not mask.any():
        raise SystemExit(f"location '{location_name}' not found in {MAGIC_ID}")
    for col in POLE_COLS:
        loc.loc[mask, col] = pole_row[col]
    return loc


def main():
    work = Path(tempfile.mkdtemp(prefix='magic19680_merge_'))
    ok, msg = ipmag.download_magic_from_id(MAGIC_ID, directory=str(work))
    if not ok:
        raise SystemExit(f'could not download contribution {MAGIC_ID}: {msg}')
    ipmag.download_magic(infile=f'magic_contribution_{MAGIC_ID}.txt',
                         dir_path=str(work), input_dir_path=str(work),
                         print_progress=False)

    # ---- sites: append the literature site means each pole needs ----------
    sites = read_magic(work / 'sites.txt')
    n_19680_sites = len(sites)

    # MIF: the two Palmer & Davis (1987) cooling-unit means (result_type='a')
    mif_clean = read_magic(MIF_DIR / 'sites.txt')
    palmer = mif_clean[mif_clean['site'].str.startswith('P_')].copy()
    palmer['location'] = MIF_LOCATION

    # SLB: the 10 Tauxe & Kodama (2009) flow sites (result_type='i'); the 40
    # Fairchild Two Island River 'mag' tc sites are already in 19680.
    slb_clean = read_magic(SLB_DIR / 'sites.txt')
    tauxe = slb_clean[slb_clean['citations'] == CIT_TAUXE].copy()
    tauxe['location'] = SLB_LOCATION

    sites = pd.concat([sites, palmer, tauxe], ignore_index=True).fillna('')
    write_magic(work / 'sites.txt', 'sites', sites)

    # ---- locations: set both poles on their respective rows --------------
    loc = read_magic(work / 'locations.txt')
    mif_pole = read_magic(MIF_DIR / 'locations.txt').iloc[0]
    slb_pole = read_magic(SLB_DIR / 'locations.txt').iloc[0]
    loc = set_pole(loc, MIF_LOCATION, mif_pole)
    loc = set_pole(loc, SLB_LOCATION, slb_pole)
    write_magic(work / 'locations.txt', 'locations', loc)

    # samples / specimens / measurements / ages pass through unchanged.
    print(f'Merged contribution {MAGIC_ID}:')
    print(f'  sites: {n_19680_sites} (19680) + {len(palmer)} Palmer means '
          f'+ {len(tauxe)} Tauxe & Kodama = {len(sites)}')
    print(f'  MIF pole on "{MIF_LOCATION}": {mif_pole["pole_lat"]}N '
          f'{mif_pole["pole_lon"]}E A95={mif_pole["pole_alpha95"]} '
          f'N={mif_pole["pole_n_sites"]}')
    print(f'  SLB pole on "{SLB_LOCATION}": {slb_pole["pole_lat"]}N '
          f'{slb_pole["pole_lon"]}E A95={slb_pole["pole_alpha95"]} '
          f'N={slb_pole["pole_n_sites"]}')

    result = ipmag.upload_magic(dir_path=str(work), input_dir_path=str(work))
    if not result[0]:
        raise SystemExit(f'upload_magic failed: {result[1]}')
    # upload_magic names the file by the contribution's locations
    # (e.g. 'Michipicoten-Island_Two-Island-River_<date>.txt'); prefix it.
    dest = MIF_DIR / ('Fairchild2017_updated_' + Path(result[0]).name)
    shutil.copy(result[0], dest)
    print(f'Wrote merged contribution: {dest} ({dest.stat().st_size/1e6:.1f} MB)')
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
