"""
Build a MagIC 3.0 sites.txt + locations.txt for the Michipicoten Island
Formation (MIF) ca. 1084 Ma pole, reproducing the careful, date-bracketed pole
of Fairchild et al. (2017): 174.7 degE, 17.0 degN, A95 4.4, N = 23 (validated).

This is the enhancement to be merged into the Fairchild et al. (2017) MagIC
contribution (earthref.org/MagIC/19680): it adds the location-level Michipicoten
Island Formation pole and the underlying Palmer & Davis (1987) cooling-unit
means that the pole requires.

Site selection (from Fairchild et al., 2017; reproduced from the data
repository github.com/Swanson-Hysell-Group/2017_Late_Rift):
    The pole uses only flows stratigraphically BRACKETED by the two dated
    horizons -- the West Sand Bay Member tuff (1084.35 Ma) below and the
    Davieaux Island Member rhyolite (1083.52 Ma) above -- i.e. the Quebec
    Harbour, South Shore, and Davieaux Island Members. The lower Cuesta and
    Channel Lake Members (below the tuff) are excluded.

    - 21 South Shore Member basalt sites: Fairchild et al. (2017),
      high-temperature component, tilt-corrected. These are taken directly from
      the published MagIC contribution 19680 (its sites table downloaded with
      ipmag.download_magic_from_id('19680') and saved here as
      Fairchild2017_19680_sites.txt; sites SS* of location 'Michipicoten
      Island'). Citation 10.1130/L580.1.
    - 2 Palmer & Davis (1987) cooling-unit means (each combining their
      site means that fall within a single cooling unit), from the recalculated
      table Palmer1987_combined_sites.csv:
        * Quebec Harbour Member = Fisher mean of Palmer sites 2 and 4
        * Davieaux Island Member = Fisher mean of Palmer sites 15 and 17
      Citation 10.1016/0301-9268(87)90077-5.
    Palmer's own South Shore sites (5, 14, 16) are not used (superseded by the
    new Fairchild South Shore data), and Palmer's Cuesta/Channel Lake sites are
    excluded as below the dated tuff.

Age: bracketed by U-Pb (CA-ID-TIMS, 206Pb/238U) zircon dates of 1084.35 +/- 0.20
    Ma (West Sand Bay tuff) and 1083.52 +/- 0.23 Ma (Davieaux Island rhyolite)
    (Fairchild et al., 2017); age_low/age_high used.

GPMDB: 9916.  MagIC data model: v3.0.
"""

from pathlib import Path
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

LOCATION = 'Michipicoten Island Formation'
CIT_FAIRCHILD = '10.1130/L580.1'
CIT_PALMER = '10.1016/0301-9268(87)90077-5'
AGE_LOW, AGE_HIGH = '1083.52', '1084.35'

# Palmer cooling-unit means used in the pole (member -> selection), from
# Palmer1987_combined_sites.csv. Each is a Fisher mean of the Palmer & Davis
# (1987) site means that fall within a single cooling unit (per the mapping of
# Fairchild et al., 2017, after Annells, 1974). The description spells out the
# derivation; the VGP is recomputed here from the combined mean direction.
PALMER_MEMBERS = {
    'Quebec_Harbour_Member': ('P_QuebecHarbour', 'Andesite',
        'Quebec Harbour Member cooling-unit mean (a combined result, not an '
        'original field site). Fisher mean of the tilt-corrected site-mean '
        'directions of Palmer & Davis (1987) sites 2 (Dec=303.0, Inc=35.0) and '
        '4 (Dec=312.0, Inc=43.0), interpreted by Fairchild et al. (2017; after '
        'Annells, 1974) to belong to a single cooling unit; combined mean '
        'Dec=307.2, Inc=39.1 (k=116.7, a95=23.0). The VGP (vgp_lat/vgp_lon) is '
        'recomputed from this mean direction at the locality (47.71N, 274.06E). '
        'Combined so the formation pole does not include multiple VGPs from one '
        'cooling unit.'),
    'Davieaux_Island_Member': ('P_Davieaux', 'Rhyolite',
        'Davieaux Island Member cooling-unit mean (a combined result, not an '
        'original field site). Fisher mean of the tilt-corrected site-mean '
        'directions of Palmer & Davis (1987) sites 15 (Dec=291.0, Inc=5.0) and '
        '17 (Dec=283.0, Inc=34.0), which belong to the single Davieaux Island '
        'rhyolite cooling unit; combined mean Dec=287.1, Inc=19.1 (k=133.0, '
        'a95=21.8). The VGP (vgp_lat/vgp_lon) is recomputed from this mean '
        'direction at the locality (47.70N, 274.21E). Combined so the formation '
        'pole does not include multiple VGPs from one cooling unit.'),
}

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

    # --- Fairchild South Shore Member sites (HT, tilt-corrected) from MagIC 19680 ---
    fc = pd.read_csv(d / 'Fairchild2017_19680_sites.txt', sep='\t', header=1)
    fc = fc[(fc['location'] == 'Michipicoten Island')
            & (fc['dir_comp_name'] == 'HT') & (fc['dir_tilt_correction'] == 100)
            & (fc['site'].astype(str).str.startswith('SS'))]
    for _, s in fc.iterrows():
        dec, inc, a95 = float(s['dir_dec']), float(s['dir_inc']), float(s['dir_alpha95'])
        lat, lon = float(s['lat']), float(s['lon'])
        _, _, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
        rows.append({
            'site': s['site'], 'location': LOCATION, 'result_type': 'i',
            'result_quality': 'g', 'method_codes': str(s.get('method_codes', 'LP-DIR-T:DE-BFL:DE-FM')),
            'citations': CIT_FAIRCHILD, 'geologic_classes': 'Igneous',
            'geologic_types': 'Lava Flow', 'lithologies': str(s.get('lithologies', 'Basalt')),
            'lat': f'{lat:.3f}', 'lon': f'{lon:.3f}',
            'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '100', 'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}',
            'dir_k': f'{float(s["dir_k"]):.1f}', 'dir_alpha95': f'{a95:.1f}',
            'dir_n_samples': str(int(float(s['dir_n_samples']))),
            'vgp_lat': f'{float(s["vgp_lat"]):.1f}', 'vgp_lon': f'{float(s["vgp_lon"]):.1f}',
            'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}',
            'description': 'South Shore Member; Fairchild et al. (2017)',
        })

    # --- Palmer & Davis (1987) cooling-unit means (Quebec Harbour, Davieaux) ---
    pa = pd.read_csv(d / 'Palmer1987_combined_sites.csv')
    pa = pa[pa['formation'] == 'Michipicoten_Island']
    for member, (name, lith, desc) in PALMER_MEMBERS.items():
        s = pa[pa['member'] == member].iloc[0]
        dec, inc, a95 = float(s['dec_tc']), float(s['inc_tc']), float(s['a95'])
        lat = float(s['site_lat']); lon = float(s['site_lon']) % 360
        vlon, vlat, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
        rows.append({
            'site': name, 'location': LOCATION, 'result_type': 'a',
            'result_quality': 'g', 'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP',
            'citations': CIT_PALMER, 'geologic_classes': 'Igneous',
            'geologic_types': 'Lava Flow', 'lithologies': lith,
            'lat': f'{lat:.3f}', 'lon': f'{lon:.3f}',
            'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '100', 'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}',
            'dir_k': f'{float(s["site_kappa"]):.1f}', 'dir_alpha95': f'{a95:.1f}',
            'dir_n_samples': str(int(float(s['n']))),
            'vgp_lat': f'{vlat:.1f}', 'vgp_lon': f'{vlon:.1f}',
            'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}',
            'description': desc,
        })

    sites_path = out / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in SITE_COLS) + '\n')
    print(f'Wrote {sites_path}: {len(rows)} site rows '
          f'({sum(r["citations"]==CIT_FAIRCHILD for r in rows)} Fairchild SS, '
          f'{sum(r["citations"]==CIT_PALMER for r in rows)} Palmer cooling-unit means)')

    pole = ipmag.fisher_mean([float(r['vgp_lon']) for r in rows],
                             [float(r['vgp_lat']) for r in rows])
    lats = [float(r['lat']) for r in rows]; lons = [float(r['lon']) for r in rows]
    loc = {
        'location': LOCATION, 'location_type': 'Outcrop',
        'result_name': 'Michipicoten Island Formation ca. 1084 Ma pole',
        'result_type': 'a', 'result_quality': 'g',
        'method_codes': 'LP-DIR-T:DE-BFL:DE-FM:DE-VGP',
        'citations': f'{CIT_FAIRCHILD}:{CIT_PALMER}',
        'geologic_classes': 'Igneous', 'lithologies': 'Basalt:Andesite:Rhyolite',
        'lat_s': f'{min(lats):.3f}', 'lat_n': f'{max(lats):.3f}',
        'lon_w': f'{min(lons):.3f}', 'lon_e': f'{max(lons):.3f}',
        'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': '100',
        'pole_lat': f'{pole["inc"]:.1f}', 'pole_lon': f'{pole["dec"]:.1f}',
        'pole_alpha95': f'{pole["alpha95"]:.1f}', 'pole_k': f'{pole["k"]:.1f}',
        'pole_n_sites': str(pole['n']),
        'sites': ':'.join(r['site'] for r in rows),
        'description': (
            'Michipicoten Island Formation mean pole (Fairchild et al., 2017). '
            'Uses only flows bracketed by the U-Pb dated horizons: the West Sand '
            'Bay Member tuff (1084.35 Ma) below and the Davieaux Island Member '
            'rhyolite (1083.52 Ma) above. Combines 21 South Shore Member basalt '
            'sites (Fairchild et al., 2017, high-temperature, tilt-corrected) '
            'with two Palmer & Davis (1987) cooling-unit means: the Quebec '
            'Harbour Member (their site means 2 and 4) and the Davieaux Island '
            'Member (their site means 15 and 17). The lower Cuesta and Channel '
            'Lake Members (below the dated tuff) and Palmer\'s own South Shore '
            'sites are excluded. N = 23.'),
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
