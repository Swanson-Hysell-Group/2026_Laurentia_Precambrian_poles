"""
Build a MagIC 3.0 sites.txt + locations.txt for the Schroeder-Lutsen basalts
(SLB) ca. 1090 Ma pole, reproducing the careful pole of Fairchild et al. (2017):
187.8 degE, 27.1 degN, A95 3.0, N = 50 (validated; this build = 27.1 / 187.8,
A95 3.0, k 46.5, N 50, matching the published value).

The pole is the Fisher mean of 50 lava-flow site VGPs:
    - 40 Schroeder-Lutsen basalt sites measured along the Two Island River near
      Schroeder, Minnesota by Fairchild et al. (2017), magnetite ('mag')
      component, tilt-corrected (SLB01-SLB40). The site MEAN DIRECTIONS are the
      published values from the study analysis file (pmag_results.csv in the
      Swanson-Hysell-Group/2017_Late_Rift study repo), archived here as
      Fairchild2017_published_site_means.csv. The site coordinates and method
      codes are carried from the published MagIC contribution 19680 sites table
      (Fairchild2017_19680_sites.txt; location 'Two Island River',
      dir_comp_name=='mag', dir_tilt_correction==100). Citation 10.1130/L580.1.

      NOTE: the site means archived in contribution 19680 differ from the
      published study values for 8 of the 40 sites (SLB08, SLB10, SLB15, SLB23,
      SLB27, SLB28, SLB31, SLB32): in 19680 each of those site means retains one
      additional specimen (n+1) that the study rejected as an outlier, which
      collapses the site precision (e.g. SLB08 k 617 -> 6, SLB32 k 271 -> 8) and
      shifts the site mean by ~0.6-8 deg. Using the published study means (here)
      reproduces the published pole exactly; using the 19680-archived means gives
      26.9/188.0 (a ~0.2 deg offset). The updated MagIC contribution corrects
      those 8 site rows and flags the rejected specimens (see
      ../../1084_Michipicoten_Island_Formation/contribution_build/build_updated_19680.py).
    - 10 Schroeder-Lutsen basalt flow sites of Tauxe & Kodama (2009)
      (ns006-ns015, tilt-corrected), the 'nsl' (above-NSVG) sequence in the
      compiled North Shore Volcanic Group data set. AF + thermal demagnetization;
      VGPs recomputed here from the tilt-corrected site mean directions.
      Citation 10.1016/j.pepi.2009.07.006.

Following the published selection (and the Michipicoten precedent of matching the
paper rather than a looser compilation), the 15 Books (1972) Schroeder sites used
in the APWP pole_means.csv (N=65 -> 28.3/187.6) are NOT included: Fairchild et
al. (2017) treated the older single-AF-step data separately and built the pole
from their 40 flows + Tauxe & Kodama's 10. Each site enters as one VGP per
cooling unit (each site is an individual lava flow), result_type='i'.

Age: the SLB are not directly dated. They unconformably overlie the North Shore
Volcanic Group and postdate the 1091.61 +/- 0.14 Ma Silver Bay aplite dike of the
Beaver Bay Complex (Fairchild et al., 2017; their best maximum age constraint),
and are likely older than the youngest dated rift volcanism (Michipicoten Island
Formation, 1083.52 +/- 0.23 Ma). age_low/age_high bracket this interval; the pole
is treated as ca. 1090 Ma.

GPMDB: 9915.  MagIC data model: v3.0.  Source study contribution: MagIC 19680.
"""

from pathlib import Path
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

LOCATION = 'Schroeder-Lutsen basalts'
CIT_FAIRCHILD = '10.1130/L580.1'
CIT_TAUXE = '10.1016/j.pepi.2009.07.006'
AGE_LOW, AGE_HIGH = '1083.52', '1091.61'

TAUXE_NSL = ['ns006', 'ns007', 'ns008', 'ns009', 'ns010',
             'ns011', 'ns012', 'ns013', 'ns014', 'ns015']

SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95',
    'dir_n_samples', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description',
]
LOC_COLS = [
    'location', 'location_type', 'result_name', 'result_type', 'result_quality',
    'method_codes', 'citations', 'geologic_classes', 'lithologies',
    'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k',
    'pole_n_sites', 'sites', 'description',
]


def main():
    d = Path(__file__).parent
    out = d.parent
    rows = []

    # --- Fairchild Two Island River SLB sites (mag, tilt-corrected) ----------
    # Direction/k/n: published study means (pmag_results.csv -> archived here).
    # Coordinates + method codes: carried from the 19680 sites table.
    pub = pd.read_csv(d / 'Fairchild2017_published_site_means.csv')
    pub = pub[(pub['comp'] == 'mag') & (pub['tilt'] == 100)].sort_values('site')
    meta = pd.read_csv(d / 'Fairchild2017_19680_sites.txt', sep='\t', header=1)
    meta = meta[(meta['location'] == 'Two Island River')
                & (meta['dir_comp_name'] == 'mag') & (meta['dir_tilt_correction'] == 100)]
    coord = meta.set_index('site')[['lat', 'lon', 'method_codes']].to_dict('index')
    for _, s in pub.iterrows():
        site = s['site']
        dec, inc, a95 = float(s['dec']), float(s['inc']), float(s['alpha95'])
        lat, lon = float(coord[site]['lat']), float(coord[site]['lon'])
        vlat, vlon = float(s['vgp_lat']), float(s['vgp_lon'])
        # dp/dm recomputed from the published mean direction at the site
        _, _, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
        rows.append({
            'site': site, 'location': LOCATION, 'result_type': 'i',
            'result_quality': 'g',
            'method_codes': str(coord[site]['method_codes']),
            'citations': CIT_FAIRCHILD, 'geologic_classes': 'Igneous',
            'geologic_types': 'Lava Flow', 'lithologies': 'Basalt',
            'lat': f'{lat:.4f}', 'lon': f'{lon:.4f}',
            'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '100', 'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}',
            'dir_k': f'{float(s["k"]):.1f}', 'dir_alpha95': f'{a95:.1f}',
            'dir_n_samples': str(int(s['n'])),
            'vgp_lat': f'{vlat:.1f}', 'vgp_lon': f'{vlon:.1f}',
            'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}',
            'description': ('Schroeder-Lutsen basalt lava flow, Two Island River '
                           'section; magnetite (low-temperature) component, '
                           'tilt-corrected; published site mean of Fairchild et '
                           'al. (2017).'),
        })
    n_fc = len(rows)

    # --- Tauxe & Kodama (2009) Schroeder-Lutsen flow sites (tilt-corrected) ---
    tk = pd.read_csv(d / 'TauxeKodama2009_sites.txt', sep='\t', header=1)
    tk = tk[(tk['site'].isin(TAUXE_NSL)) & (tk['dir_tilt_correction'] == 100)]
    tk = tk.sort_values('site')
    for _, s in tk.iterrows():
        dec, inc, a95 = float(s['dir_dec']), float(s['dir_inc']), float(s['dir_alpha95'])
        lat, lon = float(s['lat']), float(s['lon'])
        vlon, vlat, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
        rows.append({
            'site': s['site'], 'location': LOCATION, 'result_type': 'i',
            'result_quality': 'g',
            'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM',
            'citations': CIT_TAUXE, 'geologic_classes': 'Igneous',
            'geologic_types': 'Lava Flow', 'lithologies': 'Basalt',
            'lat': f'{lat:.4f}', 'lon': f'{lon:.4f}',
            'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '100', 'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}',
            'dir_k': f'{float(s["dir_k"]):.1f}', 'dir_alpha95': f'{a95:.1f}',
            'dir_n_samples': str(int(float(s['dir_n_samples']))),
            'vgp_lat': f'{vlat:.1f}', 'vgp_lon': f'{vlon:.1f}',
            'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}',
            'description': ('Schroeder-Lutsen basalt lava flow (above the North '
                           'Shore Volcanic Group); tilt-corrected site mean of '
                           'Tauxe & Kodama (2009); VGP recomputed here from the '
                           'site mean direction.'),
        })
    n_tk = len(rows) - n_fc

    sites_path = out / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in SITE_COLS) + '\n')
    print(f'Wrote {sites_path}: {len(rows)} site rows '
          f'({n_fc} Fairchild Two Island River, {n_tk} Tauxe & Kodama)')

    pole = ipmag.fisher_mean([float(r['vgp_lon']) for r in rows],
                             [float(r['vgp_lat']) for r in rows])
    lats = [float(r['lat']) for r in rows]; lons = [float(r['lon']) for r in rows]
    loc = {
        'location': LOCATION, 'location_type': 'Outcrop',
        'result_name': 'Schroeder-Lutsen basalts ca. 1090 Ma pole',
        'result_type': 'a', 'result_quality': 'g',
        'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP',
        'citations': f'{CIT_FAIRCHILD}:{CIT_TAUXE}',
        'geologic_classes': 'Igneous', 'lithologies': 'Basalt',
        'lat_s': f'{min(lats):.4f}', 'lat_n': f'{max(lats):.4f}',
        'lon_w': f'{min(lons):.4f}', 'lon_e': f'{max(lons):.4f}',
        'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': '100',
        'pole_lat': f'{pole["inc"]:.1f}', 'pole_lon': f'{pole["dec"]:.1f}',
        'pole_alpha95': f'{pole["alpha95"]:.1f}', 'pole_k': f'{pole["k"]:.1f}',
        'pole_n_sites': str(pole['n']),
        'sites': ':'.join(r['site'] for r in rows),
        'description': (
            'Schroeder-Lutsen basalts mean pole (Fairchild et al., 2017). Fisher '
            'mean of 50 lava-flow site VGPs: 40 Two Island River flows of '
            'Fairchild et al. (2017) (magnetite component, tilt-corrected) and 10 '
            'Schroeder-Lutsen flows of Tauxe & Kodama (2009). The older Books '
            '(1972) Schroeder sites are excluded, following the published '
            'selection. The VGP population is non-Fisherian, comprising two '
            'irregular clusters (one more northerly, one more southerly) at '
            'repeated stratigraphic levels; the mean pole nonetheless provides a '
            'robust late-stage rift constraint. N = 50.'),
    }
    loc_path = out / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(LOC_COLS) + '\n')
        f.write('\t'.join(str(loc[c]) for c in LOC_COLS) + '\n')
    print(f'Wrote {loc_path}: pole {loc["pole_lat"]}N {loc["pole_lon"]}E '
          f'A95={loc["pole_alpha95"]} k={loc["pole_k"]} N={loc["pole_n_sites"]}')

    result = ipmag.upload_magic(dir_path=str(out), input_dir_path=str(out))
    print(f'Created {result[0]}' if result[0] else f'Validation issue: {result[1]}')


if __name__ == '__main__':
    main()
