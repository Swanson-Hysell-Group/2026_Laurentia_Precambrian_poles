"""
Build a MagIC 3.0 sites.txt + locations.txt for the ca. 583 Ma Baie des Moutons
(Mutton Bay) syenite, Quebec, from the site-level data of McCausland et al.
(2011, Precambrian Research 187, 58-78; doi:10.1016/j.precamres.2011.02.004).

The study reports two characteristic remanence components carried by different
parts of the nested-cone alkalic complex, giving two paleomagnetic poles:

    ChRM A -- 8 'steep' syenite sites, easterly and steep, PSD magnetite.
        Mean direction D = 98.6, I = 78.0 (a95 = 6.6, k = 71.7) at the locality;
        published pole 42.6 degN, 332.7 degE (dp = 11.7, dm = 12.4). GPMDB 9364.
    ChRM B -- 6 'shallow' sites (feldspar-porphyry / aplite / carbonatite dykes
        plus the late red-syenite site MB20, whose direction is reversed and is
        inverted for inclusion). Mean direction D = 163.1, I = 6.0 (a95 = 21.7,
        k = 10.5); published pole -34.2 degN, 321.5 degE (dp = 10.9, dm = 21.8).
        GPMDB 9365.

There is no measurement-/specimen-level MagIC contribution for this study, so the
contribution is assembled from the published Table 1 site means. Both components
are in geographic coordinates -- the complex is a set of nested ring intrusions
whose varied dyke orientations indicate it has not been significantly tilted
since emplacement (McCausland et al., 2011), so no tilt correction is applied
(dir_tilt_correction = 0).

Site coordinates: Table 1 gives location as decimal minutes added to 50 deg N
and 58 deg W (NAD-27); here lat = 50 + min/60, and lon (0-360 E) = 360 - (58 +
min/60). dir_n_samples is set to the number of ChRM-bearing (endpoint + great-
circle) specimens used in each site mean (Table 1 'e+g'); the VGP and its dp/dm
are computed from the site mean direction with pmag.dia_vgp. Site MB15 has only
great-circle fits (no listed a95); its a95 is estimated as 140/sqrt(k*n).

Age: 583.4 +/- 2.0 Ma, weighted mean of replicate hornblende 40Ar-39Ar plateau
ages from unit-1 syenite (site MB44) and the late unit-3 syenite (site MB20)
(McCausland et al., 2011). Baked-contact tests at six sites were inconclusive
(no stable, non-VRM remanence in the host syenite).

MagIC data model v3.0. Two locations are written -- one per component/pole --
since the A and B site sets are disjoint and each carries its own pole result.
"""

from pathlib import Path
import math
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

CIT = '10.1016/j.precamres.2011.02.004'
AGE, AGE_SIGMA = '583.4', '2.0'
LOC_A = 'Baie des Moutons syenite (ChRM A)'
LOC_B = 'Baie des Moutons syenite (ChRM B)'
DEMAG = 'LP-DIR-AF:LP-DIR-T:LT-LT-Z:DE-BFL:DE-BFP:DE-FM'

# Table 1 of McCausland et al. (2011). Fields:
#   site, dec, inc, a95 (None if great-circle only), k, n (endpoint+great-circle
#   specimens used), lat_min (+50 N), lonW_min (+58 W), lithology, geologic_type,
#   rocktype description, intrusive unit.
SITES_A = [
    ('MB03', 109.2, 80.8, 7.2, 41.5, 11, 48.44, 58.62, 'Syenite', 'Pluton', 'foliated grey syenite', 'unit 2'),
    ('MB04',  82.9, 79.4, 7.1, 54.7,  9, 48.65, 58.61, 'Syenite', 'Pluton', 'foliated grey syenite', 'unit 2'),
    ('MB05', 141.5, 82.8, 9.3, 68.8,  5, 48.68, 58.63, 'Syenite', 'Pluton', 'pink feldspar syenite', 'unit 1'),
    ('MB15',  94.5, 64.4, None, 10.8, 5, 47.75, 59.62, 'Syenite', 'Pluton', 'grey syenite (great-circle site mean)', 'unit 2'),
    ('MB19',  46.4, 73.6, 11.4, 46.1, 5, 51.64, 59.67, 'Aplite', 'Volcanic Dike', 'aplite dyke', 'unit 1'),
    ('MB25',  78.7, 79.7, 9.7, 62.9,  5, 54.61, 57.53, 'Aplite', 'Volcanic Dike', 'aplite dyke and host', 'unit 2'),
    ('MB26', 136.0, 80.7, 8.5, 81.5,  5, 54.61, 57.53, 'Mafic Dike', 'Volcanic Dike', 'mafic dyke', 'unit 3'),
    ('MB44', 124.0, 69.2, 7.3, 46.0, 10, 49.69, 57.75, 'Syenite', 'Pluton', 'green syenite (Ar-Ar dated)', 'unit 1'),
]
SITES_B = [
    ('MB06', 183.8, -5.9, 5.9,  80.1, 9, 48.42, 59.49, 'Aplite', 'Volcanic Dike', 'aplite dyke', 'unit 1'),
    ('MB11', 162.8, -3.7, 4.4, 325.0, 5, 47.83, 59.35, 'Syenite', 'Volcanic Dike', 'red feldspar-porphyry dyke', 'dyke'),
    ('MB12', 163.2, 23.7, 9.1,  92.1, 5, 47.83, 59.35, 'Syenite', 'Volcanic Dike', 'red feldspar-porphyry dyke', 'dyke'),
    ('MB20', 321.1, 22.0, 6.7,  59.6, 9, 50.20, 58.56, 'Syenite', 'Pluton', 'red syenite (reversed; Ar-Ar dated)', 'unit 3'),
    ('MB36', 173.1,  9.7, 10.3, 25.8, 9, 48.09, 59.17, 'Syenite', 'Volcanic Dike', 'red feldspar-porphyry dyke', 'dyke'),
    ('MB43', 151.4, 33.4, 5.8, 108.0, 7, 50.48, 57.96, 'Carbonatite', 'Volcanic Dike', 'ribbon carbonatite dyke', 'dyke'),
]

SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_sigma', 'age_unit',
    'dir_tilt_correction', 'dir_comp_name', 'dir_dec', 'dir_inc', 'dir_k',
    'dir_alpha95', 'dir_n_samples', 'dir_n_specimens',
    'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description',
]
LOC_COLS = [
    'location', 'location_type', 'result_name', 'result_type', 'result_quality',
    'method_codes', 'citations', 'geologic_classes', 'lithologies',
    'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age', 'age_sigma', 'age_unit',
    'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k',
    'pole_n_sites', 'sites', 'description',
]


def site_latlon(lat_min, lonW_min):
    lat = 50.0 + lat_min / 60.0
    lon = 360.0 - (58.0 + lonW_min / 60.0)
    return lat, lon


def build_rows(sites, comp, location):
    rows = []
    for (site, dec, inc, a95, k, n, lat_min, lonW_min,
         lith, gtype, rock, unit) in sites:
        lat, lon = site_latlon(lat_min, lonW_min)
        if a95 is None:                       # great-circle-only site mean
            a95 = 140.0 / math.sqrt(k * n)
        vlon, vlat, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
        rows.append({
            'site': site, 'location': location, 'result_type': 'i',
            'result_quality': 'g', 'method_codes': DEMAG, 'citations': CIT,
            'geologic_classes': 'Igneous', 'geologic_types': gtype,
            'lithologies': lith,
            'lat': f'{lat:.4f}', 'lon': f'{lon:.4f}',
            'age': AGE, 'age_sigma': AGE_SIGMA, 'age_unit': 'Ma',
            'dir_tilt_correction': '0', 'dir_comp_name': comp,
            'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}', 'dir_k': f'{k:.1f}',
            'dir_alpha95': f'{a95:.1f}', 'dir_n_samples': str(n),
            'dir_n_specimens': str(n),
            'vgp_lat': f'{vlat:.1f}', 'vgp_lon': f'{vlon:.1f}',
            'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}',
            'description': (f'Baie des Moutons {rock} ({unit}); ChRM {comp}, '
                           f'geographic coordinates (McCausland et al., 2011).'),
        })
    return rows


def pole_row(rows, location, comp, result_name, target, description,
             pole_lat_sign=1):
    # Fisher mean of the site VGPs with polarity unified (combine antipodes),
    # matching pt.compute_mean_pole; needed for ChRM B, where the late
    # red-syenite site MB20 carries the reversed direction. The polarity-
    # unified mean is reported in the hemisphere McCausland et al. (2011) chose
    # for each component (A: northern, B: southern), flipping 180 deg if needed.
    vgp_block = ipmag.make_di_block([float(r['vgp_lon']) for r in rows],
                                    [float(r['vgp_lat']) for r in rows])
    vgp_block = pmag.flip(vgp_block, combine=True)
    pole = ipmag.fisher_mean(di_block=vgp_block)
    if (pole['inc'] >= 0) != (pole_lat_sign >= 0):
        vgp_block = ipmag.do_flip(di_block=vgp_block)
        pole = ipmag.fisher_mean(di_block=vgp_block)
    lats = [float(r['lat']) for r in rows]; lons = [float(r['lon']) for r in rows]
    loc = {
        'location': location, 'location_type': 'Outcrop',
        'result_name': result_name, 'result_type': 'a', 'result_quality': 'g',
        'method_codes': DEMAG + ':GM-ARAR', 'citations': CIT,
        'geologic_classes': 'Igneous', 'lithologies': 'Syenite',
        'lat_s': f'{min(lats):.4f}', 'lat_n': f'{max(lats):.4f}',
        'lon_w': f'{min(lons):.4f}', 'lon_e': f'{max(lons):.4f}',
        'age': AGE, 'age_sigma': AGE_SIGMA, 'age_unit': 'Ma',
        'dir_tilt_correction': '0',
        'pole_lat': f'{pole["inc"]:.1f}', 'pole_lon': f'{pole["dec"]:.1f}',
        'pole_alpha95': f'{pole["alpha95"]:.1f}', 'pole_k': f'{pole["k"]:.1f}',
        'pole_n_sites': str(pole['n']), 'sites': ':'.join(r['site'] for r in rows),
        'description': description,
    }
    print(f'  {comp}: pole {loc["pole_lat"]}N {loc["pole_lon"]}E '
          f'A95={loc["pole_alpha95"]} k={loc["pole_k"]} N={loc["pole_n_sites"]} '
          f'(target {target})')
    return loc


def main():
    out = Path(__file__).parent.parent
    rows_a = build_rows(SITES_A, 'A', LOC_A)
    rows_b = build_rows(SITES_B, 'B', LOC_B)
    rows = rows_a + rows_b

    sites_path = out / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in SITE_COLS) + '\n')
    print(f'Wrote {sites_path}: {len(rows)} sites ({len(rows_a)} ChRM A, '
          f'{len(rows_b)} ChRM B)')

    loc_a = pole_row(rows_a, LOC_A, 'A',
                     'Baie des Moutons syenite ChRM A pole (ca. 583 Ma)',
                     '42.6N/332.7E dp11.7 dm12.4',
                     ('Baie des Moutons syenite ChRM A mean pole (McCausland '
                      'et al., 2011). Fisher mean of the 8 steep, easterly '
                      'syenite-site VGPs (PSD magnetite). Geographic coordinates '
                      '(no tilt correction; the nested-cone complex is '
                      'interpreted as untilted). Baked-contact tests were '
                      'inconclusive. N = 8.'))
    loc_b = pole_row(rows_b, LOC_B, 'B',
                     'Baie des Moutons syenite ChRM B pole (ca. 583 Ma)',
                     '-34.2N/321.5E dp10.9 dm21.8',
                     ('Baie des Moutons syenite ChRM B mean pole (McCausland '
                      'et al., 2011). Fisher mean of 6 shallow site VGPs '
                      '(feldspar-porphyry, aplite and carbonatite dykes plus the '
                      'late red-syenite site MB20, whose reversed direction is '
                      'brought to common polarity), implying at least one '
                      'reversal during ChRM B acquisition. Geographic '
                      'coordinates; baked-contact tests inconclusive. N = 6.'),
                     pole_lat_sign=-1)

    loc_path = out / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(LOC_COLS) + '\n')
        for loc in (loc_a, loc_b):
            f.write('\t'.join(str(loc[c]) for c in LOC_COLS) + '\n')
    print(f'Wrote {loc_path}: 2 poles')

    result = ipmag.upload_magic(dir_path=str(out), input_dir_path=str(out))
    print(f'Created {result[0]}' if result[0] else f'Validation issue: {result[1]}')


if __name__ == '__main__':
    main()
