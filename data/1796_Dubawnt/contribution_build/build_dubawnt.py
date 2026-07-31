"""Build the Dubawnt Group MagIC contribution (sites.txt + locations.txt + ages.txt).

Source
------
``Park1973_Dubawnt_Table1.csv`` -- site-level cleaned, bedding-corrected
directions of the Dubawnt Group transcribed from Table 1 of Park, Irving &
Donaldson (1973), "Paleomagnetism of the Precambrian Dubawnt Group," Geol. Soc.
Am. Bull. 84, 859-870 (doi:10.1130/0016-7606(1973)84<859:POTPDG>2.0.CO;2). This
is a pre-MagIC study with no specimen- or measurement-level data, so the
contribution is built at the site + location level (cf. the Cleaver Dykes and
Lake Shore Traps legacy contributions).

Notes on the source data
------------------------
- Park et al. (1973) do NOT tabulate per-site alpha95/k or per-site coordinates;
  only group-mean statistics (their Table 2) and a single mean sampling location
  (64.1 N / 94.7 W = 265.3 E). All site rows therefore share that location and
  carry no dir_alpha95/dir_k.
- Per-site VGPs are computed from the cleaned bedding-corrected directions with
  pmag.dia_vgp at the mean location. The Fisher mean of the 30 accepted-site VGPs
  reproduces the published pole (7 N / 083 W = 277 E, K 12, A95 8, N 30) -- the
  recomputed value is 7.7 N / 276.5 E, K 12.0, A95 7.9. (Park's own Table-1 "Pole"
  column uses an inconsistent longitude convention and is not used.)
- Directions are bedding-corrected (dir_tilt_correction = 100); the fold test is
  only marginally positive (k 14 -> 17, not significant), and the authors prefer
  the corrected directions.

Field tests (Park et al., 1973): positive baked-contact test (dikes at sites 2
and 9 bake Hudsonian basement / Kazan sandstone; the baked rock carries the dike
direction) -> contact_test C+. Positive reversal test (normal and reversed means
4 deg apart, "reversal is exact within error") -> reversal_test +.

Age
---
Nominal 1796 Ma, bracket 1758-1833 (Baker Lake Group). The sampled Kazan and
Christopher Island Formations and the Martell Syenite are bracketed between the
U-Pb zircon age of the basal Christopher Island felsic minette flow (1833 +/- 3
Ma, 2sigma; Rainbird et al., 2006) and a U-Pb zircon date from the overlying
Wharton Group (1758 Ma). A syenite plug that cross-cuts the Kazan and Christopher
Island Formations gives an Ar-Ar age of 1811 +/- 12 Ma (Rainbird et al., 2006) --
a minimum age for the sampled units, so the pole most likely sits in the older
half of the bracket. Nominal 1796 Ma is the bracket midpoint. Park et al.'s
original K-Ar / Rb-Sr age (~1725 Ma, old decay constants) is superseded.
"""
import os
import sys
from datetime import date

import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
SRC = os.path.join(HERE, 'Park1973_Dubawnt_Table1.csv')

PARK = '10.1130/0016-7606(1973)84<859:POTPDG>2.0.CO;2'
RAINBIRD = '10.1086/498097'                     # Rainbird et al. (2006), J. Geol. 114, 1-17
RAINBIRD_DAVIS = '10.1130/B25989.1'             # Rainbird & Davis (2007), GSA Bull. 119, 314-328
SLAT, SLON = 64.1, 265.3           # mean sampling location (64.1 N / 94.7 W)
AGE, AGE_LOW, AGE_HIGH = 1796, 1758, 1833

df = pd.read_csv(SRC, comment='#')

POL = {'N': 'n', 'R': 'r', 'M': ''}     # 'm' (mixed) is not a MagIC polarity CV value -> blank
DEMAG = {'thermal': 'LP-DIR-T:DE-BLANKET:DE-FM', 'AF': 'LP-DIR-AF:DE-BLANKET:DE-FM'}

site_rows = []
for _, r in df.iterrows():
    plon, plat, dp, dm = pmag.dia_vgp(r['dec'], r['inc'], 5.0, SLAT, SLON)  # a95 placeholder; dp/dm not reported
    # sandstone/mudstone/breccia are lithologies, not geologic_types -> use the CV type 'Sediment Layer'
    gtype = 'Sediment Layer' if r['gclass'] == 'Sedimentary' else r['gtype']
    site_rows.append({
        'site': str(r['site']), 'location': 'Dubawnt Group',
        'result_type': 'i', 'result_quality': r['quality'],
        'method_codes': DEMAG[r['demag']], 'citations': PARK,
        'geologic_classes': r['gclass'], 'geologic_types': gtype,
        'lithologies': r['lith'], 'lat': SLAT, 'lon': SLON,
        'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': 100, 'dir_dec': r['dec'], 'dir_inc': r['inc'],
        'dir_n_samples': int(r['n']),
        'dir_polarity': POL.get(str(r['polarity']), ''),
        'description': f"{r['formation']}: {r['rock']}",
        'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
    })
sites = pd.DataFrame(site_rows)

good = sites[sites['result_quality'] == 'g']
# pole = Fisher mean of the accepted-site VGPs (polarity unified)
pblk = pmag.flip(ipmag.make_di_block(good['vgp_lon'].tolist(), good['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(pblk)
# mean direction (polarity unified), reported in the reversed sense Park used (D~347/I~-50)
dblk = pmag.flip(ipmag.make_di_block(good['dir_dec'].tolist(), good['dir_inc'].tolist()), combine=True)
d = pmag.fisher_mean(dblk)
if d['inc'] > 0:                          # report the up/NW (reversed) mean, as Park (347/-50)
    d['dec'] = (d['dec'] + 180) % 360
    d['inc'] = -d['inc']
# reversed fraction, matching Park's grouping (20 of 30 = the reversed + mixed sites)
good_src = df[df['quality'] == 'g']
pct_rev = round(100 * good_src['polarity'].isin(['R', 'M']).sum() / len(good_src))

locs = pd.DataFrame([{
    'location': 'Dubawnt Group', 'location_type': 'Region',
    'result_name': 'Dubawnt Group ca. 1796 Ma pole', 'result_type': 'a',
    'result_quality': 'g', 'sites': ':'.join(good['site'].tolist()),
    'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BLANKET:DE-FM:DE-VGP:ST-C', 'citations': PARK,
    'geologic_classes': 'Sedimentary:Extrusive:Intrusive',
    'lithologies': 'Sandstone:Trachyte:Syenite',
    'lat_s': SLAT, 'lat_n': SLAT, 'lon_w': SLON, 'lon_e': SLON,
    'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
    'dir_tilt_correction': 100,
    'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
    'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1), 'dir_n_sites': int(d['n']),
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
    'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'pole_reversed_perc': pct_rev, 'contact_test': 'C+', 'reversal_test': '+',
    'continent_ocean': 'North America', 'country': 'Canada',
    'state_province': 'Nunavut', 'terranes': 'Laurentia',
    'description': (
        "Paleomagnetic pole for the Dubawnt Group (Baker Lake Basin, western "
        "Churchill Province, Nunavut), from Park, Irving & Donaldson (1973). The "
        "pole is the Fisher mean of the bedding-corrected site virtual geomagnetic "
        "poles from 30 sites spanning the Kazan Formation red beds (hematite-borne "
        "remanence isolated by thermal demagnetization), the Christopher Island "
        "Formation trachyte-andesite lavas, dikes and interbedded red sandstone "
        "(specular hematite and magnetite, isolated by alternating-field "
        "demagnetization) and the Martell Syenite. The magnetization is dual "
        "polarity with normal and reversed means indistinguishable within error (a "
        "positive reversal test), and is supported as primary by a positive "
        "baked-contact test where dikes at sites 2 and 9 bake the Hudsonian "
        "basement and Kazan sandstone. Per-site alpha95/k and coordinates were not "
        "tabulated by Park et al. (1973); all sites carry the mean sampling location "
        "and no site-level precision. Age: the sampled units are bracketed between "
        "the U-Pb zircon age of the basal Christopher Island volcanic flow "
        "(1833 +/- 3 Ma; Rainbird et al., 2006) and U-Pb zircon dating of Pitz "
        "Formation rhyolite in the overlying Wharton Group (1757.6 +/- 3.4 Ma; "
        "Rainbird & Davis, 2007); a syenite plug that cross-cuts the sampled "
        "formations (Ar-Ar 1811 +/- 12 Ma; Rainbird et al., 2006) is a minimum age "
        "for the sampled units, so the pole most likely lies in the older half of "
        "the 1758-1833 Ma bracket. Nominal 1796 Ma is the bracket midpoint.")}])

# ages table: the radiometric ages that bracket the pole (Rainbird et al., 2006).
# Published uncertainties are 2 sigma; age_sigma is the 1 sigma (half).
ages = pd.DataFrame([
    {'location': 'Dubawnt Group', 'age': 1833, 'age_sigma': 1.5, 'age_unit': 'Ma',
     'timescale_eon': 'Proterozoic', 'timescale_era': 'Paleoproterozoic',
     'timescale_period': 'Orosirian', 'method_codes': 'GM-UPB',
     'timescale_citations': 'Gradstein et al. 2004', 'citations': RAINBIRD,
     'description': ('U-Pb zircon upper-intercept age of 1833 +/- 3 Ma (2 sigma) on a '
                     'basal Christopher Island Formation felsic minette flow (sample '
                     '89PHA-81; Rainbird et al., 2006) -- the oldest sampled unit and '
                     'the maximum age of the pole.')},
    {'location': 'Dubawnt Group', 'age': 1811, 'age_sigma': 5.8, 'age_unit': 'Ma',
     'timescale_eon': 'Proterozoic', 'timescale_era': 'Paleoproterozoic',
     'timescale_period': 'Orosirian', 'method_codes': 'GM-ARAR',
     'timescale_citations': 'Gradstein et al. 2004', 'citations': RAINBIRD,
     'description': ('Ar-Ar phlogopite age of 1811 +/- 12 Ma (sample 98TX-R40; Rainbird '
                     'et al., 2006) on a later syenite plug that cross-cuts the lower '
                     'Baker Sequence -- the alluvial sandstones of the Kazan Formation '
                     'and the intercalated volcanic and pyroclastic rocks of the '
                     'Christopher Island Formation, i.e. the sedimentary and volcanic '
                     'units sampled for this pole. It is distinct from the comagmatic '
                     'Martell Syenite of Park et al. (1973) and provides a minimum age '
                     'for the sampled units.')},
    {'location': 'Dubawnt Group', 'age': 1758, 'age_sigma': 1.7, 'age_unit': 'Ma',
     'timescale_eon': 'Proterozoic', 'timescale_era': 'Paleoproterozoic',
     'timescale_period': 'Statherian', 'method_codes': 'GM-UPB',
     'timescale_citations': 'Gradstein et al. 2004', 'citations': RAINBIRD_DAVIS,
     'description': ('U-Pb zircon age of 1757.6 +/- 3.4 Ma (2 sigma) on a Pitz Formation '
                     'rhyolite flow (Rainbird & Davis, 2007; a second flow gives 1753.0 '
                     '+/- 1.7 Ma). The Pitz Formation (Wharton Group) unconformably '
                     'overlies the sampled Baker Lake Group units and sets the '
                     'conservative younger bound of the pole age.')},
])


def write_magic(frame, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    frame.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
write_magic(ages, 'ages', os.path.join(OUT, 'ages.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Park1973_Dubawnt_Group_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    for frame, kind in [(locs, 'locations'), (sites, 'sites'), (ages, 'ages')]:
        f.write(f'tab delimited\t{kind}\n'); frame.to_csv(f, sep='\t', index=False)
        f.write('>>>>>>>>>>\n')

print(f'-I- Dubawnt: {len(good)} pole sites (+{len(sites) - len(good)} excluded); '
      f'pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} K {p["k"]:.1f} N {int(p["n"])}; '
      f'dir {d["dec"]:.1f}/{d["inc"]:.1f} a95 {d["alpha95"]:.1f} k {d["k"]:.1f}; %rev {pct_rev}')
print('    Park (1973) published: pole 7 N / 277 E, K 12, A95 8, N 30; dir 347/-50 a95 7')

sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from validate_magic_contribution import validate_upload_file
validate_upload_file(combined, tables=['locations', 'sites', 'ages'])
