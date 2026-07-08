"""Build the St. Francois Mountains MagIC contribution (sites.txt + locations.txt).

The pole notebook (``pole_notebooks/1466_Francois.ipynb``) is hand-maintained and
edited directly; there is no notebook builder.

Source: a student MagIC contribution (id 20670) for Meert & Stuckey (2002),
Tectonics 21(2), 1007 (doi:10.1029/2000TC001265), audited against the paper.

Audited fixes applied here (from the instructor review):
- Per-site VGPs were entered as unsigned magnitudes; recomputed signed VGPs from
  the tilt-corrected site directions with pmag.dia_vgp (the pole is 13.2 S).
- Location pole_n_sites 23 -> 18 (the mean excludes the 6-CT conglomerate test,
  the 1330 Ma dike / baked-contact rows, and reversed site 15).
- Trailing whitespace stripped from site / dir_comp_name values.
- Conglomerate-test note placed on the 6-CT row.

Second audit pass (2026-06, ready-for-upload):
- 6-CT dir_n_samples filled (13; the blank integer crashed MagIC validation).
- Location bounding box recomputed from the site coordinates so it encloses all
  sites (lat_s 37.48 -> 37.452 [17-BM]; lon_e 269.65 -> 269.7 [6-CT]).
- Age model: SFM host rocks adopt du Bray et al. (2021) 1.48-1.45 Ga episode,
  recorded as age_low/age_high 1448/1484 Ma (was 1476 +/- 16 Ma, Van Schmus 1993);
  the dike / baked-contact component rows (3-SG dike, 4-GM large dike, 4-GM baked
  contact) carry age_low/age_high 1280/1330 Ma (du Bray's younger 1.33-1.28 Ga
  magmatic pulse); the Cambrian conglomerate 6-CT keeps 500 +/- 50 Ma. du Bray
  (10.3133/pp1866) added to citations on every row whose age derives from it (SFM
  host + dike/baked-contact rows) and the location; result_name -> "ca. 1466 Ma".
- Location description documents the VGP recomputation, the 15-TS exclusion, and
  the age basis.

Published poles (Table 1): Mean-TC (from the mean direction 233.4/+36.9, k=27,
a95=6.8) = 13.2 S, 219.0 E, dp 4.7, dm 8.0 -- the abstract's headline pole, used
in the location row. Mean-VGP (Fisher mean of the 18 site VGPs) = 12.1 S, 219.1 E,
k=34, A95=6.0. Both reproduce from the contribution. Positive conglomerate,
inverse baked-contact, and fold tests.
"""
import os, io
import numpy as np
import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Meert2002_StFrancois_magic_20670_source.txt')

data = open(SRC, encoding='latin-1').read()
tables = {}
for b in data.split('>>>>>>>>>>'):
    b = b.strip()
    if not b:
        continue
    lines = b.splitlines()
    kind = lines[0].split('\t')[1].strip()
    tables[kind] = pd.read_csv(io.StringIO('\n'.join(lines[1:])), sep='\t')

sites = tables['sites'].copy()
locs = tables['locations'].copy()

# strip trailing whitespace
for c in ('site', 'dir_comp_name'):
    sites[c] = sites[c].astype(str).str.strip().replace({'nan': np.nan})

# recompute signed VGPs from each tilt-corrected site direction
def signed_vgp(row):
    if row['dir_tilt_correction'] != 100 or pd.isna(row['dir_dec']):
        return pd.Series({'vgp_lat': np.nan, 'vgp_lon': np.nan})
    a95 = row['dir_alpha95'] if not pd.isna(row['dir_alpha95']) else 5.0
    plon, plat, dp, dm = pmag.dia_vgp(row['dir_dec'], row['dir_inc'], a95, row['lat'], row['lon'])
    return pd.Series({'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1)})

vgp = sites.apply(signed_vgp, axis=1)
sites['vgp_lat'] = vgp['vgp_lat']
sites['vgp_lon'] = vgp['vgp_lon']

MEERT = '10.1029/2000TC001265'
DUBRAY = '10.3133/pp1866'  # du Bray et al. (2021) USGS PP 1866 (verified DOI)


def _is_dike(site, comp):
    """ca. 1.33-1.28 Ga dike / baked-contact component rows (excluded from the pole)."""
    site, comp = str(site).strip(), str(comp)
    return (site == '3-SG' and comp == 'b') or (site == '4-GM' and comp in ('d', 'e'))


# --- paper-audited geology (Meert & Stuckey, 2002, Table 1 footnote + text) ---
# The SFM volcanic units are ash-flow tuffs / ignimbrites, NOT lava flows: the
# Grassy Mountain unit is an ignimbrite emplaced by rheoignimbritic flow (Sides,
# 1976), and the Lake Killarney (site 8) and Taum Sauk / Bell's Mountain units at
# sites 16 and 18 are explicitly called "ash flow tuffs" (Meert & Stuckey, 2002,
# section 4). The 'b' (site 3) and 'd' (site 4) component rows are the ca. 1330 Ma
# mafic dikes, not their granite/ignimbrite host. 6-CT is the Cambrian boulder
# conglomerate (clasts of "both granitic and volcanic material").
GRANITE_UNITS = {'BH', 'KG', 'SG', 'SM'}         # plutonic granites
IGNIMBRITE_UNITS = {'GM'}                          # Grassy Mountain ignimbrite
RHYOLITE_UNITS = {'FR', 'RG', 'TS', 'BM', 'LM'}    # ash-flow tuff rhyolites


def geology_for(site, comp):
    """Return (geologic_classes, geologic_types, lithologies) audited to the paper."""
    site, comp = str(site).strip(), str(comp).strip()
    unit = site.split('-')[-1]
    if (site == '3-SG' and comp == 'b') or (site == '4-GM' and comp == 'd'):
        return ('Intrusive', 'Volcanic Dike', 'Diabase')   # ca. 1330 Ma mafic dike
    if site == '6-CT':
        return ('Sedimentary', 'Conglomerate', 'Conglomerate')
    if unit in GRANITE_UNITS:
        return ('Intrusive', 'Pluton', 'Granite')
    if unit in IGNIMBRITE_UNITS:                            # incl. 4-GM 'c'/'e' (baked)
        return ('Extrusive', 'Pyroclastic Flow', 'Ignimbrite')
    if unit in RHYOLITE_UNITS:
        return ('Extrusive', 'Pyroclastic Flow', 'Rhyolite')
    raise ValueError(f'unmapped unit for site {site!r}')


_geo = sites.apply(lambda r: pd.Series(
    geology_for(r['site'], r['dir_comp_name']),
    index=['geologic_classes', 'geologic_types', 'lithologies']), axis=1)
sites[['geologic_classes', 'geologic_types', 'lithologies']] = _geo


# --- 6-CT sample count (was blank; paper Table 1 reports 13/13) ---
m6 = sites['site'].astype(str).str.strip() == '6-CT'
sites.loc[m6 & sites['dir_n_samples'].isna(), 'dir_n_samples'] = 13

# --- 6-CT: geology set by geology_for() to the Cambrian boulder conglomerate
#     (Sedimentary / Conglomerate); its clasts are "both granitic and volcanic
#     material" (Meert & Stuckey, 2002). Age stays 500 +/- 50 Ma (Cambrian
#     depositional / conglomerate-test limiting age). ---
sites.loc[m6, 'description'] = (
    'Conglomerate-test target: Cambrian boulder conglomerate whose clasts are '
    'both granitic and volcanic material (Meert & Stuckey, 2002; depositional '
    'age ~500 Ma). The high-temperature clast directions are randomly distributed '
    '(k~1, a95~60 deg), a positive conglomerate test constraining the SFM '
    'magnetization to be older than Late Cambrian.')

# --- age model: adopt du Bray et al. (2021) 1.48-1.45 Ga SFM episode ---
pos = sites.columns.get_loc('age_sigma') + 1
sites.insert(pos, 'age_low', np.nan)
sites.insert(pos + 1, 'age_high', np.nan)

def _set_age(row):
    site = str(row['site']).strip()
    if site == '6-CT':                      # Cambrian boulder conglomerate
        return pd.Series({'age': 500.0, 'age_sigma': 50.0, 'age_low': np.nan, 'age_high': np.nan})
    if _is_dike(site, row['dir_comp_name']):  # younger 1.33-1.28 Ga magmatic pulse
        return pd.Series({'age': np.nan, 'age_sigma': np.nan, 'age_low': 1280.0, 'age_high': 1330.0})
    # SFM host rocks: nominal 1466 Ma, bounds 1448-1484 (du Bray et al., 2021) --
    # matches the notebook (nominal_age=1466, lomag/himag 1448/1484) and the Nordic
    # summary export
    return pd.Series({'age': 1466.0, 'age_sigma': np.nan, 'age_low': 1448.0, 'age_high': 1484.0})

sites[['age', 'age_sigma', 'age_low', 'age_high']] = sites.apply(_set_age, axis=1)

# du Bray DOI on every row whose age derives from it (SFM host + dike/baked-contact
# rows); the Cambrian conglomerate 6-CT (500 +/- 50 Ma) keeps Meert only
def _cite(row):
    if str(row['site']).strip() == '6-CT':
        return MEERT
    return f'{MEERT}:{DUBRAY}'

sites['citations'] = sites.apply(_cite, axis=1)

# flag rows excluded from the pole as result_quality='b' so that
# pole_tools.load_magic_sites drops them and the pole notebook can select the 18
# pole sites with no component-name filtering: the ca. 1330 Ma dike / baked-contact
# components (3-SG b, 4-GM d/e), the Cambrian conglomerate-test clasts (6-CT), and
# the reversed, poorly grouped site 15-TS. The 'c' component of the multi-component
# sites (3-SG, 4-GM) and the single-component sites remain result_quality='g'. The
# excluded rows are retained in the contribution (they support the conglomerate and
# baked-contact field tests).
def _pole_quality(row):
    site = str(row['site']).strip()
    if site in ('6-CT', '15-TS') or _is_dike(site, row['dir_comp_name']):
        return 'b'
    return 'g'

sites['result_quality'] = sites.apply(_pole_quality, axis=1)

# --- location: bounding box, age, citations ---
locs['lat_s'] = round(float(sites['lat'].min()), 3)
locs['lat_n'] = round(float(sites['lat'].max()), 3)
locs['lon_w'] = round(float(sites['lon'].min()), 3)
locs['lon_e'] = round(float(sites['lon'].max()), 3)
pos = locs.columns.get_loc('age_sigma') + 1
locs.insert(pos, 'age_low', 1448.0)
locs.insert(pos + 1, 'age_high', 1484.0)
locs['age'] = 1466.0          # nominal (matches notebook + Nordic export)
locs['age_sigma'] = np.nan
locs['result_name'] = 'St. Francois Mountains Igneous Province ca. 1466 Ma pole'
locs['citations'] = f'{MEERT}:{DUBRAY}'
# geologic classes / lithologies represented by the pole sites (the 18 good sites:
# granite + Grassy Mountain ignimbrite + ash-flow-tuff rhyolites)
_gd = sites[sites['result_quality'] == 'g']
locs['geologic_classes'] = ':'.join(sorted(_gd['geologic_classes'].dropna().unique()))
locs['lithologies'] = ':'.join(sorted(_gd['lithologies'].dropna().unique()))

# --- location pole = Fisher mean of the 18 pole-site VGPs ---
# This is the value the pole notebook computes (pt.compute_mean_pole, unify_polarity
# =False) and that flows to data/nordic_summaries/1466_St_Francois_Mountains.csv.
# It replaces the source file's Mean-TC dp/dm pole and uses the repo's standard
# pole_alpha95/pole_k columns. The `sites` list is exactly the 18 contributing sites.
pole_rows = sites[(sites['dir_tilt_correction'] == 100) & (sites['result_quality'] == 'g')]
pm = ipmag.fisher_mean(dec=pole_rows['vgp_lon'].tolist(), inc=pole_rows['vgp_lat'].tolist())
locs = locs.drop(columns=[c for c in ('pole_dp', 'pole_dm') if c in locs.columns])
ppos = locs.columns.get_loc('pole_lon') + 1
locs.insert(ppos, 'pole_alpha95', round(pm['alpha95'], 1))
locs.insert(ppos + 1, 'pole_k', round(pm['k'], 1))
locs['pole_lat'] = round(pm['inc'], 1)
locs['pole_lon'] = round(pm['dec'], 1)
locs['pole_n_sites'] = int(pm['n'])
locs['sites'] = ':'.join(pole_rows['site'].tolist())

locs['description'] = (
    'Pole from Meert & Stuckey (2002), recreated at the site level. The location '
    'pole is the Fisher mean of the tilt-corrected site VGPs. For comparison, '
    'Meert & Stuckey also report a Mean-TC pole from the mean tilt-corrected '
    'direction (233.4/+36.9): 13.2 S, 219.0 E (dp 4.7, dm 8.0). The 18 sites in '
    'the `sites` list are the stable-direction '
    'pole sites; excluded from the mean (retained in the contribution as '
    'result_quality=b for the field tests) are the reversed, poorly grouped site '
    '15-TS, the Cambrian boulder-conglomerate clasts (6-CT; conglomerate test), '
    'and the ca. 1.33-1.28 Ga dike / baked-contact component rows of 3-SG and 4-GM. '
    'Per-site VGPs were recomputed (signed) from the listed tilt-corrected '
    'directions, so they are internally consistent with the directions and may '
    'differ slightly from the rounded VGPs in Table 1 (and in hemisphere for the '
    'excluded younger-pulse rows). The pole is supported by positive fold, conglomerate, '
    'and inverse baked-contact tests. Age 1466 Ma (1448-1484) adopted from du Bray '
    'et al. (2021); originally dated 1476 +/- 16 Ma (Van Schmus et al., 1993).'
)

# --- structured field-test CV columns (Meert & Stuckey 2002): positive
#     conglomerate (G+, >500 Ma), positive inverse baked-contact (IC+, >1330 Ma),
#     and positive fold test (F+). Recorded on the pole location the tests support. ---
locs['conglomerate_test'] = 'G+'
locs['contact_test'] = 'IC+'
locs['fold_test'] = 'F+'
# fold test passes above the 99% confidence level (McFadden 1990; McElhinny 1964
# k2/k1 = 3.46 vs critical 1.78) -- Meert & Stuckey (2002)
locs['fold_test_significance'] = 99
# single-polarity pole (the one reversed site, 15-TS, is excluded from the mean)
locs['pole_reversed_perc'] = 0

# --- geographic / tectonic metadata (restored from the 2006 MagIC-team pole-only
#     contribution 14744, which carried richer location descriptors than the 2026
#     student contribution 20670; 'U.S.A.' corrected to the CV country name) ---
locs['continent_ocean'] = 'North America'
locs['country'] = 'United States of America'
locs['state_province'] = 'Missouri'
locs['terranes'] = 'Laurentia'
locs['geological_province_sections'] = 'St. Francois Mountains Igneous Province'

# --- location mean direction (tilt-corrected) underlying the pole ---
dm = ipmag.fisher_mean(dec=pole_rows['dir_dec'].tolist(), inc=pole_rows['dir_inc'].tolist())
locs['dir_dec'] = round(dm['dec'], 1)
locs['dir_inc'] = round(dm['inc'], 1)
locs['dir_alpha95'] = round(dm['alpha95'], 1)
locs['dir_k'] = round(dm['k'], 1)
locs['dir_n_sites'] = int(dm['n'])
locs['dir_tilt_correction'] = 100


# NOTE: no separate `ages` table for now. The location + site rows carry the
# age bracket (age_low/age_high); building a proper ages table around the newer
# du Bray et al. (2021) dates is deferred (a larger geochronology task).


def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


def write_combined(path):
    with open(path, 'w', encoding='utf-8') as f:
        for kind, df in [('locations', locs), ('sites', sites)]:
            if kind != 'locations':
                f.write('>>>>>>>>>>\n')
            f.write(f'tab delimited\t{kind}\n')
            df.to_csv(f, sep='\t', index=False)


def validate(path):
    import sys
    root = os.path.dirname(os.path.dirname(OUT))
    sys.path.insert(0, os.path.join(root, 'scripts'))
    from validate_magic_contribution import validate_upload_file
    return validate_upload_file(path, tables=['locations', 'sites'])


from datetime import date
write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'St.-Francois-Mountains_{stamp}.txt')
write_combined(combined)
print(f'-I- wrote sites.txt ({len(sites)} rows), locations.txt ({len(locs)} pole), '
      f'and combined {os.path.basename(combined)}')
validate(combined)
