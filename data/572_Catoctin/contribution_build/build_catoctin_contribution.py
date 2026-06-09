"""Build a MagIC contribution for the Catoctin volcanic province from Meert, Van
der Voo & Payne (1994), JGR 99, 4625-4641 (doi:10.1029/93JB01723).

All site-level data are transcribed from the paper's Tables 1 (C and B components)
and 2 (A component). Three magnetization components:
  A = primary ca. 572 Ma (magnetite) -- passes fold + reversal + baked-contact tests
  B = Late Cambrian ~505 Ma remagnetization (hematite) -- passes fold test
  C = Taconic ~450 Ma remagnetization (magnetite) -- fails fold test
Each component is written with in-situ (dir_tilt_correction=0) and, where the
paper reports them, tilt-corrected (=100) site-mean rows, with the per-site VGPs
from the tables (S latitudes negative). The C component table reports only in-situ
per-site directions/VGPs. n=1 sites (A: 6, 16) are flagged result_quality='b'.

Site coordinates are not tabulated per site; the study locality 38.5 N / 281.8 E
(Blue Ridge, central Virginia) is used (the pole uses the tabulated VGPs directly).
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
LAT, LON = 38.5, 281.8

# --- A component (Table 2): site, n, Ds, Is, Dc, Ic, k, a95, Plts, Plgs, Pltc, Plgc
A = [
    ('2',  3, 191, -81, 256, -74, 504, 6,  -56, 118, -39, 151),
    ('3',  6, 289,  52, 281,  80, 154, 6,  -33,  40, -40,  87),
    ('6',  1,  41, -79, 133, -80, None, None, -21,  97, -50,  90),
    ('8',  6, 164,  89, 164,  89,  43, 8,  -36, 113, -36, 113),
    ('9',  5,  55,  80,  86,  76, 334, 5,  -48, 137, -36, 145),
    ('15', 4, 194, -76, 203, -86,  97, 9,  -63, 125, -45, 116),
    ('16', 1, 297,  82, 159,  82, None, None, -44,  91, -24, 117),
    ('28', 4, 189, -75, 248, -75, 730, 4,  -66, 123, -44, 147),
]
# --- B component (Table 1): site, n, Ds, Is, Dc, Ic, k, a95, Plts, Plgs, Pltc, Plgc
B = [
    ('1',  5, 100,  43, 102,  28,  87,  8,  -8, 176,  0, 183),
    ('4',  6, 103,  29,  94,  31,  13, 19,  -1, 182, -7, 186),
    ('12', 5,  70,  29,  79,  18,  65, 12, -26, 202, -15, 202),
    ('17', 4, 107,  -4, 106, -16,  20, 25, -14,  13, -17,  18),
    ('18', 6, 251,  -8, 252,  -6, 100,  8, -18, 211, -16, 211),
    ('20', 4, 102,  -3, 101,  16,  20, 22, -11,  15,  -3,   9),
    ('25', 7,  51,  49,  77,  20,  27, 10, -48, 198, -17, 202),
    ('26', 5, 117, -64, 104, -35,  26, 15, -12, 152,  -2, 178),
    ('30', 8, 281, -30, 279,  -5,  21, 12,  -1, 183,  -4,  14),
]
# --- C component (Table 1, in-situ only): site, n, Ds, Is, k, a95, Plt, Plg
C = [
    ('1',  3, 150,  72,  55, 13,   9, 308),
    ('4',  3, 152,  34,  24, 25,  27, 140),
    ('6',  4, 111,  55,  33, 17,   8, 342),
    ('10', 3, 150,  31,  70, 10,  28, 144),
    ('11', 8, 164,  28,  51,  7,  35, 131),
    ('13', 4, 131,  59,  53, 24,   0, 147),
    ('16', 3, 137,  55,  36, 21,   6, 146),
    ('19', 4, 157,  29, 154,  9,  32, 138),
    ('21', 3, 146,  30,  21, 26,  27, 148),
    ('22', 17, 143, 58,  27,  7,  21, 168),
    ('23', 8, 155,  44,  50,  7,  21, 136),
    ('24', 6, 159,  45,  30, 14,  22, 132),
    ('29', 6, 129,  21,  58,  7,  21, 168),
    ('30', 6, 144,  53,  17, 16,  10, 142),
    ('31', 6, 157,  37,  24, 14,  27, 136),
]

rows = []
def add(site, comp, tc, dec, inc, k, a95, n, vlat, vlon, quality='g'):
    rows.append({'site': f'{site}{comp}', 'location': 'Catoctin volcanic province',
                 'dir_comp_name': comp, 'result_type': 'i', 'result_quality': quality,
                 'method_codes': 'LP-DIR-T:LP-DIR-AF:DE-BFL:DE-FM:DE-VGP',
                 'citations': '10.1029/93JB01723', 'geologic_classes': 'Extrusive:Igneous',
                 'geologic_types': 'Lava Flow', 'lithologies': 'Basalt',
                 'lat': LAT, 'lon': LON, 'dir_tilt_correction': tc,
                 'dir_dec': dec, 'dir_inc': inc, 'dir_k': k, 'dir_alpha95': a95,
                 'dir_n_samples': n, 'vgp_lat': vlat, 'vgp_lon': vlon,
                 'age': 568 if comp == 'A' else (505 if comp == 'B' else 450),
                 'age_unit': 'Ma'})

for s, n, ds, is_, dc, ic, k, a95, plts, plgs, pltc, plgc in A:
    q = 'b' if n == 1 else 'g'                 # exclude n=1 sites 6, 16 from the A pole
    add(s, 'A', 0, ds, is_, k, a95, n, plts, plgs, q)
    add(s, 'A', 100, dc, ic, k, a95, n, pltc, plgc, q)
for s, n, ds, is_, dc, ic, k, a95, plts, plgs, pltc, plgc in B:
    add(s, 'B', 0, ds, is_, k, a95, n, plts, plgs)
    add(s, 'B', 100, dc, ic, k, a95, n, pltc, plgc)
for s, n, ds, is_, k, a95, plt, plg in C:
    add(s, 'C', 0, ds, is_, k, a95, n, plt, plg)

sites = pd.DataFrame(rows)

locs = pd.DataFrame([
    {'location': 'Catoctin volcanic province', 'location_type': 'Region',
     'result_name': 'Catoctin A component (primary) ca. 572 Ma pole', 'result_type': 'a',
     'sites': '2A:3A:8A:9A:15A:28A', 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM:DE-VGP:ST-C:FT-LT:RT',
     'citations': '10.1029/93JB01723', 'geologic_classes': 'Extrusive:Igneous', 'lithologies': 'Basalt',
     'lat_s': LAT, 'lat_n': LAT, 'lon_w': LON, 'lon_e': LON,
     'age': 568, 'age_low': 555, 'age_high': 577, 'age_unit': 'Ma', 'dir_tilt_correction': 100,
     'pole_lat': -42.4, 'pole_lon': 126.9, 'pole_alpha95': 16.5, 'pole_k': 17.4, 'pole_n_sites': 6,
     'description': 'Primary A-component pole (magnetite, ca. 572 Ma); Fisher mean of 6 tilt-corrected site VGPs (n>1). Reported by Meert et al. (1994) as 43 S/128 E (= 42.4 N/306.9 E). Positive fold, reversal, and baked-contact tests.'},
    {'location': 'Catoctin volcanic province', 'location_type': 'Region',
     'result_name': 'Catoctin B component Late Cambrian (~505 Ma) remagnetization pole', 'result_type': 'a',
     'sites': '1B:4B:12B:17B:18B:20B:25B:26B:30B', 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM:DE-VGP:FT-LT',
     'citations': '10.1029/93JB01723', 'geologic_classes': 'Extrusive:Igneous', 'lithologies': 'Basalt',
     'lat_s': LAT, 'lat_n': LAT, 'lon_w': LON, 'lon_e': LON,
     'age': 505, 'age_low': 500, 'age_high': 515, 'age_unit': 'Ma', 'dir_tilt_correction': 100,
     'pole_lat': -4.0, 'pole_lon': 193.0, 'pole_alpha95': 13.0, 'pole_k': 16.0, 'pole_n_sites': 9,
     'description': 'B-component pole (hematite), tilt-corrected mean D=92/I=+17; pole 4 S/193 E. A Late Cambrian (Carolinian-Penobscottian) remagnetization that passes the fold test; secondary, not the Catoctin pole.'},
    {'location': 'Catoctin volcanic province', 'location_type': 'Region',
     'result_name': 'Catoctin C component Taconic (~450 Ma) remagnetization pole', 'result_type': 'a',
     'sites': '1C:4C:6C:10C:11C:13C:16C:19C:21C:22C:23C:24C:29C:30C:31C',
     'method_codes': 'LP-DIR-T:DE-BFL:DE-FM:DE-VGP', 'citations': '10.1029/93JB01723',
     'geologic_classes': 'Extrusive:Igneous', 'lithologies': 'Basalt',
     'lat_s': LAT, 'lat_n': LAT, 'lon_w': LON, 'lon_e': LON,
     'age': 450, 'age_low': 450, 'age_high': 470, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
     'pole_lat': 19.0, 'pole_lon': 143.0, 'pole_alpha95': 9.0, 'pole_k': 21.0, 'pole_n_sites': 15,
     'description': 'C-component pole (magnetite), in-situ mean D=147/I=+44; pole 19 N/143 E. A Taconic (~450 Ma) greenschist remagnetization that FAILS the fold test; secondary, not the Catoctin pole.'},
])

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)} rows: '
      f"A={sum(sites.dir_comp_name=='A')}, B={sum(sites.dir_comp_name=='B')}, "
      f"C={sum(sites.dir_comp_name=='C')}), locations.txt ({len(locs)} poles)")
