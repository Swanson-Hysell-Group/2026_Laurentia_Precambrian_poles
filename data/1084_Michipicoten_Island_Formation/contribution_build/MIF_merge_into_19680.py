"""
Produce the merged, updated Fairchild et al. (2017) MagIC contribution (19680).

Starts from the published contribution earthref.org/MagIC/19680 (downloaded with
ipmag.download_magic_from_id) and adds, without altering anything else:
  - the two Palmer & Davis (1987) cooling-unit-mean SITES (P_QuebecHarbour,
    P_Davieaux), flagged result_type='a' (averages of multiple site means), and
  - the location-level Michipicoten Island Formation mean POLE on the
    'Michipicoten Island' location (17.0 degN, 174.7 degE, A95 4.4, k 47.9,
    N = 23), matching Fairchild et al. (2017).

All other tables of 19680 -- including the full measurements table -- pass
through unchanged. The 21 South Shore sites that enter the pole are already in
19680. The added site rows and the pole are taken from the sibling files
../sites.txt and ../locations.txt built by MIF_build_contribution.py.

Output: a single merged MagIC upload file written into the pole folder, ready to
upload to MagIC as an update of contribution 19680.

Note: this downloads ~2.4 MB (incl. 12,095 measurement rows) from MagIC; run
with a network connection.
"""

import shutil
import tempfile
from pathlib import Path

import pandas as pd
import pmagpy.ipmag as ipmag

MAGIC_ID = '19680'
MIF_LOCATION_IN_19680 = 'Michipicoten Island'   # 19680's location name for the MIF sites
HERE = Path(__file__).parent
POLE_DIR = HERE.parent


def read_magic(path):
    return pd.read_csv(path, sep='\t', skiprows=1, dtype=str).fillna('')


def write_magic(path, table_type, df):
    with open(path, 'w') as f:
        f.write(f'tab delimited\t{table_type}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')


def main():
    work = Path(tempfile.mkdtemp(prefix='mif_merge_'))
    ok, msg = ipmag.download_magic_from_id(MAGIC_ID, directory=str(work))
    if not ok:
        raise SystemExit(f'could not download contribution {MAGIC_ID}: {msg}')
    ipmag.download_magic(infile=f'magic_contribution_{MAGIC_ID}.txt',
                         dir_path=str(work), input_dir_path=str(work),
                         print_progress=False)

    # --- sites: append the two Palmer cooling-unit means (result_type='a') ---
    sites = read_magic(work / 'sites.txt')
    clean_sites = read_magic(POLE_DIR / 'sites.txt')
    palmer = clean_sites[clean_sites['site'].str.startswith('P_')].copy()
    palmer['location'] = MIF_LOCATION_IN_19680            # align to 19680 location
    sites = pd.concat([sites, palmer], ignore_index=True).fillna('')
    write_magic(work / 'sites.txt', 'sites', sites)
    n_palmer = len(palmer)

    # --- locations: put the MIF mean pole on the 'Michipicoten Island' row ---
    loc = read_magic(work / 'locations.txt')
    pole = read_magic(POLE_DIR / 'locations.txt').iloc[0]
    pole_cols = ['result_type', 'result_name', 'method_codes',
                 'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95',
                 'pole_k', 'pole_n_sites', 'sites', 'description', 'citations']
    for col in pole_cols:
        if col not in loc.columns:
            loc[col] = ''
    mask = loc['location'] == MIF_LOCATION_IN_19680
    if not mask.any():
        raise SystemExit(f"location '{MIF_LOCATION_IN_19680}' not found in 19680")
    loc.loc[mask, 'result_type'] = 'a'
    for col in ['result_name', 'method_codes', 'dir_tilt_correction', 'pole_lat',
                'pole_lon', 'pole_alpha95', 'pole_k', 'pole_n_sites', 'sites',
                'description', 'citations']:
        loc.loc[mask, col] = pole[col]
    write_magic(work / 'locations.txt', 'locations', loc)

    # samples / specimens / measurements / ages pass through unchanged.
    print(f'Merged: 19680 sites ({len(sites) - n_palmer}) + {n_palmer} Palmer '
          f'means; pole on "{MIF_LOCATION_IN_19680}" = '
          f'{pole["pole_lat"]}N {pole["pole_lon"]}E A95={pole["pole_alpha95"]} '
          f'N={pole["pole_n_sites"]}')

    result = ipmag.upload_magic(dir_path=str(work), input_dir_path=str(work))
    if not result[0]:
        raise SystemExit(f'upload_magic failed: {result[1]}')
    dest = POLE_DIR / ('Fairchild2017_updated_' + Path(result[0]).name)
    shutil.copy(result[0], dest)
    print(f'Wrote merged contribution: {dest} ({dest.stat().st_size/1e6:.1f} MB)')
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
