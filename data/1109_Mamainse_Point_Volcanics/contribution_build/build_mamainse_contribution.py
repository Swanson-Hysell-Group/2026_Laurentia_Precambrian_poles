"""Assemble the Mamainse Point Volcanics contribution (four poles).

The Mamainse Point succession (eastern Lake Superior Basin) is the most continuous
Keweenawan record of rift volcanism and geomagnetic reversals (Swanson-Hysell et
al., 2014a, Geology, doi 10.1130/G35271.1; MagIC contribution 16333). Paleomagnetic
poles are calculated from stratigraphic subsets of the site VGPs, following
Swanson-Hysell et al. (2009, 2014a) and the APWP_StratModels compilation, giving
four poles (oldest to youngest), spanning three geomagnetic reversals:

- Mamainse lower R1 — lowermost ~600 m; older Alona Bay reversed-polarity zone.
  height < 600 m. N = 24. Target 227.0 degE / 49.5 degN, A95 5.3.
- Mamainse lower R2 — younger Alona Bay reversed-polarity zone.
  1070 < height < 1350 m. N = 14. Target 205.2 degE / 37.5 degN, A95 4.5.
- Mamainse Flour Bay (lower N + upper R) — the Flour Bay normal- and
  reversed-polarity zones combined (latent stage). (1350<h<1810) + (1860<h<2100).
  N = 24. Target 189.7 degE / 36.1 degN, A95 4.9. Dual polarity.
- Mamainse upper N — Portage Lake normal-polarity zone, above the ~300 m "Great
  Conglomerate". height > 2400 m. N = 34. Target 183.2 degE / 31.2 degN, A95 2.5.

The three flows between ~1810 and ~1860 m (heights 1820, 1833, 1858) are the
Flour Bay normal->reversed transition and are excluded from the stable-polarity
poles. The Flour Bay pole combines a normal and a reversed zone whose VGPs are
antipodal, so polarities are unified before averaging (pmag.flip(combine=True)).

Outputs ../sites.txt (96 pole sites, location = pole grouping) and ../locations.txt
(four poles), validating each pole against its target.
"""
import os
import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'SH2014a_sites_source.txt')   # mirrors MagIC 16333

POLES = [   # (location, height-mask fn, dir_polarity, age, age_high, age_low, target)
    ('Mamainse lower R1',  lambda h: h < 600,                              'r', 1109, 1112, 1106,  (227.0, 49.5)),
    ('Mamainse lower R2',  lambda h: (h > 1070) & (h < 1350),              'r', 1105, 1109, 1100.4, (205.2, 37.5)),
    ('Mamainse Flour Bay', lambda h: ((h > 1350) & (h < 1810)) | ((h > 1860) & (h < 2100)), 'mixed', 1100.36, 1100.61, 1100.11, (189.7, 36.1)),
    ('Mamainse upper N',   lambda h: h > 2400,                             'n', 1094, 1100, 1090,  (183.2, 31.2)),
]

mp = pd.read_csv(SRC, sep='\t', header=1)
mp['citations'] = '10.1130/G35271.1'

groups, loc_rows = [], []
for name, mask, pol, age, ahi, alo, target in POLES:
    g = mp[mask(mp['height'])].copy()
    g['location'] = name
    if pol != 'mixed':
        g['dir_polarity'] = pol
    g['age'] = age; g['age_high'] = ahi; g['age_low'] = alo; g['age_unit'] = 'Ma'
    g['geologic_classes'] = 'Extrusive:Igneous'; g['geologic_types'] = 'Lava Flow'
    g['lithologies'] = 'Basalt'; g['result_quality'] = 'g'; g['result_type'] = 'i'
    groups.append(g)

    blk = pmag.flip(ipmag.make_di_block(g['vgp_lon'].tolist(), g['vgp_lat'].tolist()),
                    combine=True)
    p = pmag.fisher_mean(blk)
    print('   %-20s plon=%.1f plat=%.1f A95=%.1f K=%.1f N=%d  (target %.1f/%.1f)'
          % (name, p['dec'], p['inc'], p['alpha95'], p['k'], p['n'], target[0], target[1]))
    loc_rows.append({
        'location': name, 'location_type': 'Outcrop',
        'age': age, 'age_high': ahi, 'age_low': alo, 'age_unit': 'Ma',
        'citations': '10.1130/G35271.1', 'geologic_classes': 'Extrusive:Igneous',
        'lithologies': 'Basalt', 'lat_n': g['lat'].max(), 'lat_s': g['lat'].min(),
        'lon_e': g['lon'].max(), 'lon_w': g['lon'].min(), 'dir_tilt_correction': 100,
        'result_type': 'a', 'result_quality': 'g',
        'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
        'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
        'pole_n_sites': int(p['n']), 'result_name': f'{name} pole',
    })

sites = pd.concat(groups, ignore_index=True)
sites['dir_tilt_correction'] = 100
out_sites = os.path.join(OUT, 'sites.txt')
with open(out_sites, 'w') as f:
    f.write('tab\tsites\n')
sites.to_csv(out_sites, sep='\t', index=False, mode='a')
print(f'-I- wrote {out_sites} ({len(sites)} sites)')

out_locs = os.path.join(OUT, 'locations.txt')
with open(out_locs, 'w') as f:
    f.write('tab\tlocations\n')
pd.DataFrame(loc_rows).to_csv(out_locs, sep='\t', index=False, mode='a')
print(f'-I- wrote {out_locs} ({len(loc_rows)} poles)')
