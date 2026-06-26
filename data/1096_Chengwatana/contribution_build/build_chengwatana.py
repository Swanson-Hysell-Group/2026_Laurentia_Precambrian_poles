"""Build the Chengwatana Volcanics MagIC contribution (data/1096_Chengwatana/)
from Kean, Williams & Feeney (1997), Geophys. Res. Lett. 24, 1523-1526
(doi:10.1029/97gl00993), Table 1.

Six sites of Keweenawan-age Chengwatana plateau basalts near St. Croix Falls,
Polk County, Wisconsin. Sites 1 and 2 each carry a normal (C1) and a reversed
(C2) ChRM direction (a short reversal in the otherwise normal sequence); sites
3-6 are normal only. The combined-polarity pole (8 site directions, 60
specimens) is 30.9 N / 186.1 E, A95 8.0 (Kean et al., 1997, "average N and R").
Directions are in situ (not corrected for the regional ~10-15 deg SW dip; the
paper notes the correction makes little difference). Per-site VGPs are recomputed
here with pmag.dia_vgp and reproduce Table 1.
"""
import os
import pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
DOI = '10.1029/97gl00993'

# Kean et al. (1997) Table 1: site, dec, inc, k, a95, N(specimens), site_lat(N), site_lon(W), polarity
TABLE1 = [
    ('SC1n', 288.9,  55.0, 49.5,  6.9, 10, 45.33, 92.65, 'n'),
    ('SC1r', 114.8, -51.8,  5.2, 29.2,  7, 45.33, 92.65, 'r'),
    ('SC2n', 274.0,  51.1, 13.0, 17.0,  7, 45.37, 92.67, 'n'),
    ('SC2r', 106.7, -41.5,  6.7, 23.1,  8, 45.37, 92.67, 'r'),
    ('SC3',  280.0,  33.2, 14.0, 20.8,  5, 45.39, 92.66, 'n'),
    ('SC4',  285.1,  31.0, 16.3, 14.1,  8, 45.39, 92.66, 'n'),
    ('SC5',  306.0,  45.5, 11.2, 18.8,  7, 45.36, 92.64, 'n'),
    ('SC6',  287.9,  42.6, 25.4, 11.2,  8, 45.37, 92.63, 'n'),
]

rows = []
for site, dec, inc, k, a95, n, lat, lonW, pol in TABLE1:
    lon = round(360 - lonW, 3)                    # 0-360 deg east
    plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
    rows.append({
        'site': site, 'location': 'Chengwatana Volcanics', 'result_type': 'i',
        'result_quality': 'g', 'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP',
        'citations': DOI, 'geologic_classes': 'Extrusive:Igneous',
        'geologic_types': 'Lava Flow', 'lithologies': 'Basalt',
        'lat': lat, 'lon': lon, 'age': 1095, 'age_low': 1093, 'age_high': 1097, 'age_unit': 'Ma',
        'dir_tilt_correction': 0, 'dir_dec': dec, 'dir_inc': inc, 'dir_alpha95': a95,
        'dir_k': k, 'dir_n_samples': n, 'dir_comp_name': 'ChRM', 'dir_polarity': pol,
        'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
        'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1),
        'description': f'Chengwatana flow site ({"normal" if pol == "n" else "reversed"} ChRM); Kean et al. (1997) Table 1'})

sites = pd.DataFrame(rows)

# pole from the 8 site VGPs (polarity unified) for the locations table
import pmagpy.ipmag as ipmag
blk = pmag.flip(ipmag.make_di_block(sites['vgp_lon'].tolist(), sites['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)

locs = pd.DataFrame([{
    'location': 'Chengwatana Volcanics', 'location_type': 'Outcrop',
    'result_name': 'Chengwatana Volcanics ca. 1095 Ma pole', 'result_type': 'a',
    'result_quality': 'g', 'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP:RT',
    'citations': DOI, 'geologic_classes': 'Extrusive:Igneous', 'lithologies': 'Basalt',
    'lat_s': 45.33, 'lat_n': 45.39, 'lon_w': 267.33, 'lon_e': 267.37,
    'age': 1095, 'age_low': 1093, 'age_high': 1097, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1), 'pole_alpha95': round(p['alpha95'], 1),
    'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'description': 'Chengwatana Volcanics combined-polarity pole, Fisher mean of 8 site VGPs (sites 1-2 dual polarity, 3-6 normal); 60 specimens; in situ (uncorrected for the regional SW dip). Reproduces Kean et al. (1997) 30.9N/186.1E, A95 8.0.'}])


def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f"-I- wrote sites.txt ({len(sites)} sites, {int(sites['dir_n_samples'].sum())} specimens), "
      f"locations.txt; pole {p['inc']:.1f}/{p['dec']:.1f} A95 {p['alpha95']:.1f} N {int(p['n'])}")
