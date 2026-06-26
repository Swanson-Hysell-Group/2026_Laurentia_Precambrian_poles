"""
Build a MagIC 3.0 sites.txt + locations.txt for the ca. 1107 Ma Coldwell Complex
paleomagnetic pole, from the site-mean directions of Kulakov, Smirnov & Diehl
(2014, J. Geophys. Res. 119, 8633-8654; Table 2). There is no measurement-level
MagIC contribution for this study, so the contribution is assembled from the
published Table 2 site means.

The multiphase alkaline Coldwell Complex (north shore of Lake Superior, Ontario)
was emplaced during Midcontinent Rift magmatism and is grouped by Kulakov et al.
into three intrusive centers (after Lewchuk & Symons, 1990):
  - Center A (eastern gabbro + ferroaugite syenite): 14 accepted sites, REVERSED ChRM.
  - Center B (western gabbro + syenite): 10 accepted sites, NORMAL polarity ChRM.
  - Center C (central biotite gabbro, nepheline syenite + syenite, Geordie Lakes
    area): 16 accepted sites, REVERSED ChRM.
40 sites (238 samples) met the acceptance criteria. The reversed directions of
Centers A and C are statistically indistinguishable and are combined into the
preferred pole (Pole CCr): D = 114.8, I = -63.7 (a95 3.6, k 54, N = 30) ->
VGP 47.2 N, 206.5 E (A95 4.8, K 31). Center B (normal) gives a separate group
mean direction (D = 298.0, I = 56.9; pole 44.9 N, 193.2 E); the normal and
reversed groups can be compared with a reversal test.

Polarity convention here: Table 2 lists inclinations as positive magnitudes. The
reversed Center A/C ChRM is stored with NEGATIVE inclination (up, the published
compilation lists I = -63.7); Center B normal with POSITIVE inclination. VGPs are
computed with pmag.dia_vgp from the signed directions (reversed -> southern VGP);
the pole is the polarity-unified Fisher mean, reproducing 47.2 N / 206.5 E.

Age: 1108 +/- 1 Ma and 1107 +5/-1 Ma U-Pb zircon/baddeleyite on the eastern
gabbro (Heaman & Machado, 1992); the complex spans the ~1102-1105 Ma reversal.
GPMDB 9838.  MagIC data model v3.0.
"""

from pathlib import Path
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

LOCATION = 'Coldwell Complex'
CIT = 'Kulakov2014a'                 # 10.1002/2014JB011463
AGE, AGE_LOW, AGE_HIGH = '1107', '1105', '1109'

# Table 2 accepted site means (Kulakov et al., 2014).
# (center, site, lat_N, lon_W, n_samples, n_measured, method, dec, inc_magnitude, a95, k)
# Centers A and C are REVERSED (polarity 'r'); Center B is NORMAL ('n').
SITES = [
    # --- Center A (eastern gabbro + ferrosyenite); reversed ---
    ('A', 'CLD1', 48.713, 86.322, 3, 7, 'AF+T', 116.0, 61.5, 7.4, 187),
    ('A', 'CLD3', 48.719, 86.673, 8, 9, 'AF+T', 125.7, 73.8, 8.6, 37),
    ('A', 'CLD4', 48.727, 86.328, 14, 17, 'AF+T', 103.7, 66.4, 4.4, 76),
    ('A', 'CLD5', 48.729, 86.329, 7, 11, 'AF+T', 116.9, 73.3, 4.4, 76),
    ('A', 'CLD6', 48.742, 86.334, 9, 10, 'AF+T', 106.3, 67.2, 8.5, 34),
    ('A', 'NMG1', 48.847, 86.484, 9, 10, 'AF+T', 99.1, 51.9, 12.3, 32),
    ('A', 'SL10', 48.768, 86.374, 5, 6, 'AF+T', 118.5, 70.7, 6.4, 114),
    ('A', 'SL14', 48.780, 86.420, 6, 11, 'AF+T', 125.5, 49.2, 11.0, 32),
    ('A', 'SL15', 48.773, 86.393, 5, 12, 'T', 118.4, 68.5, 11.0, 43),
    ('A', 'SL33', 48.785, 86.428, 8, 9, 'AF+T', 118.6, 66.7, 3.9, 174),
    ('A', 'SL34', 48.726, 86.327, 4, 7, 'AF+T', 138.5, 75.6, 9.1, 73),
    ('A', 'SL5-6', 48.732, 86.330, 8, 8, 'AF+T', 98.0, 63.5, 6.7, 61),
    ('A', 'SL7', 48.756, 86.310, 5, 6, 'AF+T', 113.7, 68.8, 6.6, 110),
    ('A', 'SL7-8', 48.755, 86.311, 4, 6, 'AF', 127.4, 62.4, 11.0, 52),
    # --- Center B (central syenites); normal ---
    ('B', 'CCW1', 48.814, 86.673, 5, 11, 'T', 267.0, 53.2, 14.3, 23),
    ('B', 'CCW2', 48.817, 86.684, 8, 11, 'AF+T', 293.3, 51.0, 5.9, 77),
    ('B', 'CLD7', 48.770, 86.549, 4, 12, 'AF', 318.0, 57.3, 9.3, 73),
    ('B', 'SL13E', 48.796, 86.453, 8, 8, 'T', 303.5, 62.6, 4.8, 120),
    ('B', 'SL18', 48.799, 86.651, 6, 8, 'AF+T', 299.2, 49.8, 6.0, 107),
    ('B', 'SL19', 48.796, 86.646, 5, 7, 'AF+T', 295.5, 58.5, 10.0, 46),
    ('B', 'SL20', 48.800, 86.636, 6, 8, 'AF+T', 308.2, 59.3, 14.0, 21),
    ('B', 'SL37', 48.802, 86.626, 6, 8, 'T', 321.3, 54.4, 4.9, 154),
    ('B', 'SL47', 48.817, 86.682, 7, 8, 'AF', 287.7, 58.7, 10.8, 28),
    ('B', 'SL49', 48.797, 86.651, 6, 8, 'AF', 289.6, 55.1, 5.2, 142),
    # --- Center C (western quartz syenites/syenites/granites, Geordie Lakes); reversed ---
    ('C', 'CLD11', 48.784, 86.596, 4, 7, 'AF', 129.4, 49.2, 7.7, 107),
    ('C', 'GLS', 48.824, 86.485, 5, 5, 'AF+T', 92.1, 57.2, 14.9, 23),
    ('C', 'GLG', 48.821, 86.485, 4, 6, 'AF+T', 131.5, 43.8, 8.4, 90),
    ('C', 'GLS2', 48.824, 86.485, 3, 4, 'AF', 110.2, 49.6, 7.9, 162),
    ('C', 'CLD8', 48.765, 86.525, 3, 6, 'AF', 131.2, 64.5, 6.3, 259),
    ('C', 'SL11', 48.792, 86.487, 5, 10, 'AF+T', 77.1, 68.9, 4.8, 205),
    ('C', 'SL13', 48.793, 86.463, 6, 6, 'AF+T', 120.4, 60.5, 4.4, 191),
    ('C', 'SL14-13', 48.779, 86.420, 8, 9, 'AF+T', 89.8, 69.1, 3.1, 284),
    ('C', 'SL21', 48.789, 86.611, 6, 8, 'AF+T', 120.7, 57.5, 7.0, 77),
    ('C', 'SL25', 48.793, 86.501, 4, 6, 'AF', 110.6, 57.0, 1.7, 2315),
    ('C', 'SL28', 48.802, 86.626, 6, 12, 'AF+T', 118.6, 45.4, 7.7, 64),
    ('C', 'SL32', 48.795, 86.495, 8, 8, 'AF+T', 99.4, 75.2, 6.4, 67),
    ('C', 'SL41', 48.797, 86.445, 6, 8, 'AF', 132.2, 67.9, 8.3, 56),
    ('C', 'SL44', 48.764, 86.526, 3, 9, 'AF+T', 114.0, 59.7, 7.5, 56),
    ('C', 'SL45', 48.772, 86.513, 8, 10, 'AF', 108.5, 69.7, 9.1, 38),
    ('C', 'SL46', 48.795, 86.497, 3, 9, 'AF', 114.1, 76.3, 5.9, 292),
]

CENTER_NAME = {
    'A': 'Center A (eastern gabbro and ferroaugite syenite); reversed ChRM',
    'B': 'Center B (western gabbro and syenite); normal-polarity ChRM',
    'C': 'Center C (central biotite gabbro, nepheline syenite and syenite, Geordie Lakes area); reversed ChRM',
}
METHOD_CODES = {
    'AF': 'LP-DIR-AF:DE-BFL:DA-DIR-GEO',
    'T': 'LP-DIR-T:DE-BFL:DA-DIR-GEO',
    'AF+T': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DA-DIR-GEO',
}

SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'dir_comp_name', 'dir_dec', 'dir_inc', 'dir_polarity',
    'dir_k', 'dir_alpha95', 'dir_n_samples',
    'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description',
]
LOC_COLS = [
    'location', 'location_type', 'result_name', 'result_type', 'result_quality',
    'method_codes', 'citations', 'geologic_classes', 'lithologies',
    'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k',
    'pole_n_sites', 'sites', 'description',
]


def main():
    out = Path(__file__).parent.parent
    rows = []
    for center, site, lat, lonw, n, nmeas, method, dec, inc_mag, a95, k in SITES:
        lon = 360.0 - lonw
        polarity = 'r' if center in ('A', 'C') else 'n'
        inc = -inc_mag if polarity == 'r' else inc_mag      # signed inclination
        vlon, vlat, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
        rows.append({
            'site': site, 'location': LOCATION, 'result_type': 'i',
            'result_quality': 'g', 'method_codes': METHOD_CODES[method],
            'citations': CIT, 'geologic_classes': 'Igneous',
            'geologic_types': 'Intrusion', 'lithologies': 'Alkaline Igneous Rock',
            'lat': f'{lat:.3f}', 'lon': f'{lon:.3f}',
            'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '0', 'dir_comp_name': 'ChRM',
            'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}', 'dir_polarity': polarity,
            'dir_k': f'{k:.1f}', 'dir_alpha95': f'{a95:.1f}', 'dir_n_samples': str(n),
            'vgp_lat': f'{vlat:.1f}', 'vgp_lon': f'{vlon:.1f}',
            'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}',
            'description': CENTER_NAME[center],
        })

    sites_path = out / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in SITE_COLS) + '\n')
    n_a = sum(r['description'].startswith('Center A') for r in rows)
    n_b = sum(r['description'].startswith('Center B') for r in rows)
    n_c = sum(r['description'].startswith('Center C') for r in rows)
    n_samp = sum(int(r['dir_n_samples']) for r in rows)
    print(f'Wrote {sites_path}: {len(rows)} sites '
          f'(A={n_a} rev, B={n_b} norm, C={n_c} rev; {n_samp} samples)')

    # preferred pole = Centers A+C (reversed), polarity-unified VGP Fisher mean
    ac = [r for r in rows if r['dir_polarity'] == 'r']
    vgp_block = ipmag.make_di_block([float(r['vgp_lon']) for r in ac],
                                    [float(r['vgp_lat']) for r in ac])
    vgp_block = pmag.flip(vgp_block, combine=True)
    pole = ipmag.fisher_mean(di_block=vgp_block)
    lats = [float(r['lat']) for r in rows]
    lons = [float(r['lon']) for r in rows]
    loc = {
        'location': LOCATION, 'location_type': 'Outcrop',
        'result_name': 'Coldwell Complex ca. 1107 Ma pole (Centers A+C reversed, Pole CCr)',
        'result_type': 'a', 'result_quality': 'g',
        'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DA-DIR-GEO:DE-VGP',
        'citations': CIT, 'geologic_classes': 'Igneous',
        'lithologies': 'Alkaline Igneous Rock',
        'lat_s': f'{min(lats):.3f}', 'lat_n': f'{max(lats):.3f}',
        'lon_w': f'{min(lons):.3f}', 'lon_e': f'{max(lons):.3f}',
        'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': '0',
        'pole_lat': f'{pole["inc"]:.1f}', 'pole_lon': f'{pole["dec"]:.1f}',
        'pole_alpha95': f'{pole["alpha95"]:.1f}', 'pole_k': f'{pole["k"]:.1f}',
        'pole_n_sites': str(pole['n']),
        'sites': ':'.join(r['site'] for r in ac),
        'description': (
            'Coldwell Complex mean pole (Kulakov et al., 2014, Pole CCr). Fisher '
            'mean of the 30 reversed-polarity site VGPs of Centers A and C '
            '(eastern gabbro/ferrosyenite + western syenites/granites), which give '
            'statistically indistinguishable group mean directions (D = 114.8, '
            'I = -63.7, a95 3.6, k 54). The 10 normal-polarity Center B sites '
            'define a separate group mean (D = 298.0, I = 56.9; pole 44.9 N, '
            '193.2 E). Geographic coordinates (intrusive, no tilt correction). '
            'N = 30 sites.'),
    }
    loc_path = out / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(LOC_COLS) + '\n')
        f.write('\t'.join(str(loc[c]) for c in LOC_COLS) + '\n')
    print(f'Wrote {loc_path}: pole {loc["pole_lat"]}N {loc["pole_lon"]}E '
          f'A95={loc["pole_alpha95"]} k={loc["pole_k"]} N={loc["pole_n_sites"]} '
          f'(published CCr 47.2N/206.5E A95 4.8 K 31)')


if __name__ == '__main__':
    main()
