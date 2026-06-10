"""
Build a MagIC 3.0 sites.txt + locations.txt for the ca. 615 Ma Long Range dyke
swarm, southeast Labrador, from the dyke-mean directions of Murthy et al. (1992,
Can. J. Earth Sci. 29, 1224-1234).

Six northeast-trending mafic (diabase) dykes of the Long Range suite were
studied. Murthy et al. interpret **three dykes (2, 4, 6)** as carrying a primary
remanence (steep, southeasterly and down), giving a combined mean of D = 124.8,
I = 55.5 (k = 48, a95 = 18.0; their "Mean characteristic remanence, N = 3 dykes")
and a paleopole at 10.8 degS, 164.3 degE. The other **three dykes (1, 3, 5)**
carry "anomalous" (disturbed) remanences and are excluded from the pole, but are
retained here (result_quality='b') so the notebook can show them for exposition.

This contribution is assembled from the published Table 1 dyke means (there is no
measurement-level MagIC contribution). All directions are in geographic
coordinates (dir_tilt_correction = 0; steeply-dipping dykes, no bedding
correction). The swarm locality is ~53.7 degN, 56.7 degW (303.3 degE); per-dyke
VGPs and dp/dm are computed with pmag.dia_vgp. dir_n_samples is the total number
of samples averaged in each dyke (Table 1 'N' summed over its sampling sites);
the three pole dykes (2, 4, 6) total 69 samples.

NOTE on a prior error corrected here: the earlier extraction listed dyke 6 as
D = 124.8, I = 55.5 -- that is actually the *combined* 3-dyke mean, not dyke 6.
Dyke 6's own mean (its two sites PR86-008 = 117.0/66.0 and PR86-015 = 132.6/68.6)
is D = 124.6, I = 67.5. With the correct dyke 6, the Fisher mean of dykes 2, 4, 6
reproduces the paper exactly (124.8 / 55.5, k 48, a95 18.0).

Age: 615 +/- 2 Ma (U-Pb zircon + baddeleyite on the longest dyke, dyke 1; Kamo
et al., 1989). K-Ar minimum ages of 514 +/- 8 Ma and 553 +/- 22 Ma come from two
other dykes. A positive baked-contact test was obtained for dyke 2.
GPMDB 6934-6936.  MagIC data model v3.0.
"""

from pathlib import Path
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

LOCATION = 'Long Range Dykes'
CIT = 'Murthy1992'                # no DOI in references.bib for this 1992 CJES paper
LAT, LON = 53.7, 303.3            # swarm locality (SE Labrador)
DEMAG = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-BFP:DE-FM'
AGE, AGE_LOW, AGE_HIGH = '615', '613', '617'
POLE_DYKES = {'2', '4', '6'}      # primary remanence -> the pole

SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95',
    'dir_n_samples', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description',
]
LOC_COLS = [
    'location', 'location_type', 'result_name', 'result_type', 'result_quality',
    'method_codes', 'citations', 'geologic_classes', 'lithologies',
    'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k',
    'pole_n_sites', 'sites', 'description',
]


def main():
    d = Path(__file__).parent
    out = d.parent
    src = pd.read_csv(d / 'Murthy1992_dyke_means.csv')

    rows = []
    for _, s in src.iterrows():
        dyke = str(s['dyke']).strip()
        dec, inc, k, a95 = float(s['dec']), float(s['inc']), float(s['k']), float(s['a95'])
        n = int(s['n_samples'])
        primary = dyke in POLE_DYKES
        vlon, vlat, dp, dm = pmag.dia_vgp(dec, inc, a95, LAT, LON)
        note = (f'Dyke {dyke} mean ({int(s["n_sites"])} sampling site(s), '
                f'{n} samples); Murthy et al. (1992) Table 1. ')
        note += ('Primary remanence (interpreted), used in the pole.' if primary
                 else 'Anomalous (disturbed) remanence; excluded from the pole.')
        rows.append({
            'site': f'Dyke {dyke}', 'location': LOCATION, 'result_type': 'i',
            'result_quality': 'g' if primary else 'b',
            'method_codes': DEMAG, 'citations': CIT,
            'geologic_classes': 'Igneous', 'geologic_types': 'Volcanic Dike',
            'lithologies': 'Diabase',
            'lat': f'{LAT:.2f}', 'lon': f'{LON:.2f}',
            'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '0', 'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}',
            'dir_k': f'{k:.1f}', 'dir_alpha95': f'{a95:.1f}', 'dir_n_samples': str(n),
            'vgp_lat': f'{vlat:.1f}', 'vgp_lon': f'{vlon:.1f}',
            'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}', 'description': note,
        })
    # write in dyke order 1..6 for readability
    rows.sort(key=lambda r: int(r['site'].split()[-1]))

    sites_path = out / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in SITE_COLS) + '\n')
    n_pole = sum(r['result_quality'] == 'g' for r in rows)
    print(f'Wrote {sites_path}: {len(rows)} dykes ({n_pole} primary -> pole, '
          f'{len(rows) - n_pole} anomalous -> excluded)')

    pole_rows = [r for r in rows if r['result_quality'] == 'g']
    vgp_block = ipmag.make_di_block([float(r['vgp_lon']) for r in pole_rows],
                                    [float(r['vgp_lat']) for r in pole_rows])
    vgp_block = pmag.flip(vgp_block, combine=True)
    # report in the southern hemisphere, matching Murthy et al.'s published
    # paleopole polarity (10.8 S, 164.3 E) rather than its northern antipode
    vgp_block = ipmag.do_flip(di_block=vgp_block)
    pole = ipmag.fisher_mean(di_block=vgp_block)
    lats = [float(r['lat']) for r in pole_rows]; lons = [float(r['lon']) for r in pole_rows]
    loc = {
        'location': LOCATION, 'location_type': 'Outcrop',
        'result_name': 'Long Range dykes ca. 615 Ma pole (primary dykes 2, 4, 6)',
        'result_type': 'a', 'result_quality': 'g',
        'method_codes': DEMAG + ':DE-VGP', 'citations': CIT,
        'geologic_classes': 'Igneous', 'lithologies': 'Diabase',
        'lat_s': f'{min(lats):.2f}', 'lat_n': f'{max(lats):.2f}',
        'lon_w': f'{min(lons):.2f}', 'lon_e': f'{max(lons):.2f}',
        'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': '0',
        'pole_lat': f'{pole["inc"]:.1f}', 'pole_lon': f'{pole["dec"]:.1f}',
        'pole_alpha95': f'{pole["alpha95"]:.1f}', 'pole_k': f'{pole["k"]:.1f}',
        'pole_n_sites': str(pole['n']),
        'sites': ':'.join(r['site'] for r in pole_rows),
        'description': (
            'Long Range dykes mean pole (Murthy et al., 1992). Fisher mean of the '
            'three primary-remanence dyke VGPs (dykes 2, 4, 6); the combined dyke '
            'mean direction is D = 124.8, I = 55.5 (k = 48, a95 = 18.0). Geographic '
            'coordinates. The anomalous dykes 1, 3, 5 are excluded. Reported in the '
            'southern hemisphere, matching the published paleopole at '
            '10.8S, 164.3E. N = 3 dykes.'),
    }
    loc_path = out / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(LOC_COLS) + '\n')
        f.write('\t'.join(str(loc[c]) for c in LOC_COLS) + '\n')
    print(f'Wrote {loc_path}: pole {loc["pole_lat"]}N {loc["pole_lon"]}E '
          f'A95={loc["pole_alpha95"]} k={loc["pole_k"]} N={loc["pole_n_sites"]} '
          f'(paper antipode 10.8S/164.3E)')

    result = ipmag.upload_magic(dir_path=str(out), input_dir_path=str(out))
    print(f'Created {result[0]}' if result[0] else f'Validation issue: {result[1]}')


if __name__ == '__main__':
    main()
