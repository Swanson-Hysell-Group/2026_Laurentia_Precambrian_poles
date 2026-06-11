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
    - sets the pole: 27.1 degN, 187.8 degE, A95 3.0, k 46.4, N 50
    - the 40 Two Island River 'mag' tilt-corrected sites in the pole are already
      in 19680
    (built by SLB_build_contribution.py ->
     ../../1090_Schroeder_Lutsen_Basalts/{sites,locations}.txt)

  Schroeder-Lutsen SITE-MEAN CORRECTION (8 sites)
    The site means archived in 19680 differ from the published Fairchild et al.
    (2017) study values for 8 Two Island River sites (SLB08, SLB10, SLB15, SLB23,
    SLB27, SLB28, SLB31, SLB32): in 19680 each retained one additional specimen
    (n+1) that the study rejected as a directional outlier (a near-antipodal /
    excursional sample direction, e.g. SLB08.5a = 153/-13 vs the flow mean
    ~282/55). Including it collapses the site precision (e.g. SLB08 k 617 -> 6)
    and shifts the site mean by ~0.6-8 deg, moving the published 27.1/187.8 pole
    to 26.9/188.0. This build corrects those 8 sites to the published study means
    (all four mag/hem x tilt-0/100 rows, from
    ../../1090_Schroeder_Lutsen_Basalts/contribution_build/Fairchild2017_published_site_means.csv)
    and flags the 8 rejected specimens result_quality='b' so the site means are
    consistent with the specimen data. The excluded sample at each site is
    identified reproducibly: it is the one whose removal reproduces the published
    site mean direction and precision exactly (off by <0.05 deg; verified).

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
HERE = Path(__file__).resolve().parent
MIF_DIR = HERE.parent                                  # data/1084_Michipicoten_Island_Formation
SLB_DIR = MIF_DIR.parent / '1090_Schroeder_Lutsen_Basalts'

# 19680 location names that carry each pole
MIF_LOCATION = 'Michipicoten Island'
SLB_LOCATION = 'Two Island River'

CIT_TAUXE = '10.1016/j.pepi.2009.07.006'               # Tauxe & Kodama (2009) sites

# Published Fairchild et al. (2017) Schroeder-Lutsen site means (all mag/hem x
# tilt-0/100), archived from the study analysis file (pmag_results.csv).
PUB_MEANS = (SLB_DIR / 'contribution_build' / 'Fairchild2017_published_site_means.csv')

# The 8 Two Island River sites whose 19680 site means retain one rejected
# specimen, with the outlier SAMPLE the published study excluded. Each outlier is
# the one sample whose removal reproduces the published site mean (dir + k + n)
# exactly; its specimen is <sample>+'a'. See module docstring.
SLB_OUTLIER_SAMPLES = {
    'SLB08': 'SLB08.5', 'SLB10': 'SLB10.9', 'SLB15': 'SLB15.7', 'SLB23': 'SLB23.2',
    'SLB27': 'SLB27.1', 'SLB28': 'SLB28.7', 'SLB31': 'SLB31.1', 'SLB32': 'SLB32.5',
}
SLB_LOCATION_TIR = 'Two Island River'

POLE_COLS = ['result_type', 'result_name', 'method_codes', 'citations',
             'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95',
             'pole_k', 'pole_n_sites', 'sites', 'description']


def _drop_from_colon_list(value, sample):
    """Remove a sample and its specimen (<sample>a) from a ':'-delimited list."""
    keep = [p for p in str(value).split(':')
            if p and p != sample and p != sample + 'a']
    return ':'.join(keep)


def correct_slb_site_means(sites):
    """Correct the 8 Schroeder-Lutsen site means to the published study values.

    Updates the mag/hem x tilt-0/100 rows of each affected Two Island River site
    (direction, k, alpha95, n, dir_r, VGP, and the samples/specimens lists) to the
    published Fairchild et al. (2017) means, dropping the rejected outlier sample.
    Edits ``sites`` in place and returns it.
    """
    pub = pd.read_csv(PUB_MEANS).set_index(['site', 'comp', 'tilt'])
    n_rows = 0
    for site, outlier in SLB_OUTLIER_SAMPLES.items():
        for comp in ('mag', 'hem'):
            for tilt in (0, 100):
                row = pub.loc[(site, comp, tilt)]
                mask = ((sites['site'] == site)
                        & (sites['dir_comp_name'] == comp)
                        & (sites['dir_tilt_correction'] == str(tilt)))
                if mask.sum() != 1:
                    raise SystemExit(
                        f'expected 1 site row for {site}/{comp}/tilt{tilt}, '
                        f'found {mask.sum()}')
                n = int(row['n'])
                k = float(row['k'])
                R = n - (n - 1) / k if k > 0 else float(n)
                sites.loc[mask, 'dir_dec'] = '%.1f' % row['dec']
                sites.loc[mask, 'dir_inc'] = '%.1f' % row['inc']
                sites.loc[mask, 'dir_k'] = '%.0f' % k
                sites.loc[mask, 'dir_alpha95'] = '%.1f' % row['alpha95']
                sites.loc[mask, 'dir_n_samples'] = str(n)
                sites.loc[mask, 'dir_n_specimens'] = str(n)
                sites.loc[mask, 'dir_r'] = '%.4f' % R
                sites.loc[mask, 'vgp_lat'] = '%.1f' % row['vgp_lat']
                sites.loc[mask, 'vgp_lon'] = '%.1f' % row['vgp_lon']
                sites.loc[mask, 'vgp_dp'] = '%.1f' % row['vgp_dp']
                sites.loc[mask, 'vgp_dm'] = '%.1f' % row['vgp_dm']
                sites.loc[mask, 'samples'] = sites.loc[mask, 'samples'].apply(
                    lambda v: _drop_from_colon_list(v, outlier))
                sites.loc[mask, 'specimens'] = sites.loc[mask, 'specimens'].apply(
                    lambda v: _drop_from_colon_list(v, outlier))
                n_rows += 1
    return sites, n_rows


def flag_slb_outlier_specimens(specimens):
    """Flag the 8 rejected Schroeder-Lutsen outlier specimens result_quality='b'.

    Edits ``specimens`` in place and returns (specimens, n_specimens_flagged).
    """
    out_specs = [s + 'a' for s in SLB_OUTLIER_SAMPLES.values()]
    mask = specimens['specimen'].isin(out_specs)
    specimens.loc[mask, 'result_quality'] = 'b'
    return specimens, specimens.loc[mask, 'specimen'].nunique()


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

    # ---- sites: correct the 8 SLB site means, append literature means -----
    sites = read_magic(work / 'sites.txt')
    n_19680_sites = len(sites)

    # Correct the 8 Two Island River site means to the published study values
    # (and drop the rejected outlier sample from each).
    sites, n_corrected = correct_slb_site_means(sites)

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

    # ---- specimens: flag the 8 rejected SLB outlier specimens 'b' ---------
    specimens = read_magic(work / 'specimens.txt')
    specimens, n_flagged = flag_slb_outlier_specimens(specimens)
    write_magic(work / 'specimens.txt', 'specimens', specimens)

    # ---- locations: set both poles on their respective rows --------------
    loc = read_magic(work / 'locations.txt')
    mif_pole = read_magic(MIF_DIR / 'locations.txt').iloc[0]
    slb_pole = read_magic(SLB_DIR / 'locations.txt').iloc[0]
    loc = set_pole(loc, MIF_LOCATION, mif_pole)
    loc = set_pole(loc, SLB_LOCATION, slb_pole)
    write_magic(work / 'locations.txt', 'locations', loc)

    # samples / measurements / ages pass through unchanged.
    print(f'Merged contribution {MAGIC_ID}:')
    print(f'  SLB site means corrected: {n_corrected} rows '
          f'({len(SLB_OUTLIER_SAMPLES)} sites x mag/hem x tilt-0/100); '
          f'{n_flagged} outlier specimens flagged result_quality=b')
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
