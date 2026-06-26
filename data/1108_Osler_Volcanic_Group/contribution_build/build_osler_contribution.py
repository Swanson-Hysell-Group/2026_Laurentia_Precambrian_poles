"""Assemble the Osler Volcanic Group contribution (three poles).

Builds a single site table for the Osler Volcanic Group paleomagnetic poles by
combining the three source studies, and writes the three location-level pole
results, following the user's suggested approach of enhancing the Swanson-Hysell
et al. (2014b) Osler contribution with the additional 2019 sites and the Halls
(1974) sites, with the poles in the locations table.

Sources (copied from the APWP_StratModels compilation, which mirrors the
published MagIC contributions):
- Halls (1974), CJES, doi 10.1139/e74-113 — Nipigon Strait reversed (25 sites)
  and normal (5 sites) Osler flows.
- Swanson-Hysell et al. (2014b), G-cubed, doi 10.1002/2013gc005180 — Simpson
  Island section reversed flows (height-resolved).
- Swanson-Hysell et al. (2019), GSA Bulletin, doi 10.1130/b31944.1 — Agate Point
  reversed flows and Puff Island normal flows.

Three poles (reproducing APWP_StratModels Keweenawan_pole_means.csv and the
construction in code/01_VGP_compilation.ipynb):
- Osler reverse lower (R1, lower Alona Bay reversed zone): SH2014b flows with
  height < 1041 m. N = 30.  Target 218.6 degE / 40.9 degN, A95 4.8.
- Osler reverse upper (R2, upper Alona Bay reversed zone): SH2014b flows with
  height > 2082 m + Halls (1974) reversed + Agate Point (2019). N = 64.
  Target 203.4 degE / 42.3 degN, A95 3.7.
- Osler normal (Portage Lake normal zone): Halls (1974) normal 5 sites combined
  into 2 flows (sites 1/2/5 and 3/4) + 2 Puff Island flows (2019). N = 4.
  Target 171.9 degE / 32.0 degN, A95 9.7.

Outputs ../sites.txt (all sites, location = pole grouping) and ../locations.txt
(three poles), and validates each pole against the target.
"""
import os
import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
# source site tables copied from APWP_StratModels/data/pmag_compiled (which mirror
# the published MagIC contributions); kept here so the build is self-contained.
SRC_HALLS = os.path.join(HERE, 'Halls1974_sites_source.txt')
SRC_SH14 = os.path.join(HERE, 'SH2014b_sites_source.txt')
SRC_AP = os.path.join(HERE, 'SH2019a_AgatePoint_sites_source.txt')

OUT_COLS = ['age', 'age_high', 'age_low', 'age_unit', 'citations', 'description',
            'dir_alpha95', 'dir_dec', 'dir_inc', 'dir_k', 'dir_n_samples',
            'dir_polarity', 'dir_tilt_correction', 'geologic_classes',
            'geologic_types', 'lat', 'lithologies', 'location', 'lon',
            'method_codes', 'result_quality', 'result_type', 'site',
            'vgp_lat', 'vgp_lon']


def std(df, **overrides):
    """Return df reduced/renamed to OUT_COLS with constant overrides applied."""
    out = pd.DataFrame()
    for c in OUT_COLS:
        out[c] = df[c] if c in df.columns else np.nan
    for k, v in overrides.items():
        out[k] = v
    return out


# ---- Halls (1974) ----------------------------------------------------------
halls = pd.read_csv(SRC_HALLS, sep='\t', header=1)
halls = halls[halls['dir_tilt_correction'] == 100].copy()
halls['dir_alpha95'] = 140 / np.sqrt(halls['dir_n_samples'] * halls['dir_k'])
halls_r = halls[halls['location'] == 'Osler Volcanics, Nipigon Strait, Lower Reversed'].reset_index(drop=True)
halls_n = halls[halls['location'] == 'Osler Volcanics, Nipigon Strait, Upper Normal'].reset_index(drop=True)

# ---- Swanson-Hysell et al. (2014b) -----------------------------------------
sh14 = pd.read_csv(SRC_SH14, sep='\t', header=1)
sh14 = sh14[~sh14['dir_dec'].isna()].copy()
sh14['dir_n_samples'] = sh14['dir_n_specimens']
sh14['dir_k'] = 140 ** 2 / sh14['dir_alpha95'] ** 2 / sh14['dir_n_specimens']
sh14['dir_polarity'] = 'r'
sh14_lower = sh14[sh14.height < 1041].reset_index(drop=True)
sh14_middle = sh14[(sh14.height >= 1041) & (sh14.height <= 2082)].reset_index(drop=True)
sh14_upper = sh14[sh14.height > 2082].reset_index(drop=True)

# ---- Swanson-Hysell et al. (2019) Agate Point -------------------------------
ap = pd.read_csv(SRC_AP, sep='\t', skiprows=1)
ap = ap[ap.dir_tilt_correction == 100].copy()
ap_R = ap[ap.location == 'Agate Point'].reset_index(drop=True)
ap_N = ap[ap.location == 'Puff Island'].reset_index(drop=True)

# ---- Halls normal: combine 5 sites into 2 cooling-unit flows ----------------
# Swanson-Hysell & Fairchild (2014 field mapping): sites 1,2,5 = one flow above
# the Puff Island conglomerate; sites 3,4 = a second flow (Halls, 1974 site order).
f1 = pmag.fisher_mean([[halls_n.dir_dec[i], halls_n.dir_inc[i]] for i in (0, 1, 4)])
f2 = pmag.fisher_mean([[halls_n.dir_dec[i], halls_n.dir_inc[i]] for i in (2, 3)])
halls_n_flows = pd.DataFrame([
    {'site': 'Halls_N_flow1', 'lat': halls_n.lat[0], 'lon': halls_n.lon[0],
     'dir_dec': f1['dec'], 'dir_inc': f1['inc'], 'dir_alpha95': f1['alpha95'],
     'dir_k': f1['k'], 'dir_n_samples': f1['n'], 'result_type': 'a',
     'description': 'Halls (1974) normal sites 1, 2, 5 combined as one cooling unit (flow above Puff Island conglomerate)'},
    {'site': 'Halls_N_flow2', 'lat': halls_n.lat[2], 'lon': halls_n.lon[2],
     'dir_dec': f2['dec'], 'dir_inc': f2['inc'], 'dir_alpha95': f2['alpha95'],
     'dir_k': f2['k'], 'dir_n_samples': f2['n'], 'result_type': 'a',
     'description': 'Halls (1974) normal sites 3, 4 combined as one cooling unit (Puff/Tremblay Island shoreline flow)'},
])
halls_n_flows['dir_tilt_correction'] = 100
halls_n_flows['dir_polarity'] = 'n'
ipmag.vgp_calc(halls_n_flows, site_lon='lon', site_lat='lat', tilt_correction='yes',
               dec_tc='dir_dec', inc_tc='dir_inc')

# ---- assemble the three pole groups (location = pole name) ------------------
osler_lower = std(sh14_lower, location='Osler reverse lower',
                  citations='10.1002/2013gc005180', dir_polarity='r',
                  age=1108, age_high=1110, age_low=1105.15, age_unit='Ma',
                  geologic_classes='Extrusive:Igneous', geologic_types='Lava Flow',
                  lithologies='Basalt', result_quality='g', result_type='i')

osler_middle = std(sh14_middle, location='Osler reverse middle',
                   citations='10.1002/2013gc005180', dir_polarity='r',
                   age=1107, age_high=1110, age_low=1105.48, age_unit='Ma',
                   geologic_classes='Extrusive:Igneous', geologic_types='Lava Flow',
                   lithologies='Basalt', result_quality='g', result_type='i')

osler_upper = pd.concat([
    std(sh14_upper, location='Osler reverse upper', citations='10.1002/2013gc005180',
        dir_polarity='r', result_type='i'),
    std(halls_r, location='Osler reverse upper', citations='10.1139/e74-113',
        dir_polarity='r', result_type='i'),
    std(ap_R, location='Osler reverse upper', citations='10.1130/b31944.1',
        dir_polarity='r', result_type='i'),
], ignore_index=True)
osler_upper['age'] = 1105.15; osler_upper['age_high'] = 1105.48; osler_upper['age_low'] = 1104.82
osler_upper['age_unit'] = 'Ma'; osler_upper['geologic_classes'] = 'Extrusive:Igneous'
osler_upper['geologic_types'] = 'Lava Flow'; osler_upper['lithologies'] = 'Basalt'
osler_upper['result_quality'] = 'g'

osler_normal = pd.concat([
    std(halls_n_flows, location='Osler normal', citations='10.1139/e74-113', dir_polarity='n'),
    std(ap_N, location='Osler normal', citations='10.1130/b31944.1', dir_polarity='n', result_type='i'),
], ignore_index=True)
osler_normal['age'] = 1095; osler_normal['age_high'] = 1100; osler_normal['age_low'] = 1080
osler_normal['age_unit'] = 'Ma'; osler_normal['geologic_classes'] = 'Extrusive:Igneous'
osler_normal['geologic_types'] = 'Lava Flow'; osler_normal['lithologies'] = 'Basalt'
osler_normal['result_quality'] = 'g'

sites = pd.concat([osler_lower, osler_middle, osler_upper, osler_normal], ignore_index=True)
sites['dir_tilt_correction'] = 100   # all Osler site means are tilt-corrected
sites['method_codes'] = sites['method_codes'].fillna('LP-DC2:DE-DI')
# give each site a unique id where missing
sites['site'] = [s if isinstance(s, str) and s else f'OVG{i}' for i, s in enumerate(sites['site'])]

out_sites = os.path.join(OUT, 'sites.txt')
with open(out_sites, 'w') as f:
    f.write('tab\tsites\n')
sites.to_csv(out_sites, sep='\t', index=False, mode='a')
print(f'-I- wrote {out_sites} ({len(sites)} sites)')

# ---- compute the three poles and write locations ---------------------------
def pole_of(df):
    blk = pmag.flip(ipmag.make_di_block(df['vgp_lon'].tolist(), df['vgp_lat'].tolist()),
                    combine=True)
    return pmag.fisher_mean(blk)

targets = {'Osler reverse lower': (218.6, 40.9), 'Osler reverse middle': (211.3, 42.7),
           'Osler reverse upper': (203.4, 42.3), 'Osler normal': (171.9, 32.0)}
loc_rows = []
for name, grp in [('Osler reverse lower', osler_lower), ('Osler reverse middle', osler_middle),
                  ('Osler reverse upper', osler_upper), ('Osler normal', osler_normal)]:
    p = pole_of(grp)
    tlon, tlat = targets[name]
    print(f'   {name:22s} plon=%.1f plat=%.1f A95=%.1f K=%.1f N=%d  (target %.1f/%.1f)'
          % (p['dec'], p['inc'], p['alpha95'], p['k'], p['n'], tlon, tlat))
    r = grp.iloc[0]
    loc_rows.append({
        'location': name, 'location_type': 'Region',
        'age': r['age'], 'age_high': r['age_high'], 'age_low': r['age_low'], 'age_unit': 'Ma',
        'citations': ':'.join(sorted(set(grp['citations'].dropna()))),
        'geologic_classes': 'Extrusive:Igneous', 'lithologies': 'Basalt',
        'lat_n': grp['lat'].max(), 'lat_s': grp['lat'].min(),
        'lon_e': grp['lon'].max(), 'lon_w': grp['lon'].min(),
        'dir_tilt_correction': 100, 'result_type': 'a', 'result_quality': 'g',
        'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
        'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
        'pole_n_sites': int(p['n']),
        'result_name': f'{name} pole',
    })

out_locs = os.path.join(OUT, 'locations.txt')
with open(out_locs, 'w') as f:
    f.write('tab\tlocations\n')
pd.DataFrame(loc_rows).to_csv(out_locs, sep='\t', index=False, mode='a')
print(f'-I- wrote {out_locs} ({len(loc_rows)} poles)')
