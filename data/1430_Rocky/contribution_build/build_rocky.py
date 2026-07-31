"""Build the ca. 1430 Ma Rocky Mountain intrusions MagIC contribution (1430_Rocky).

Three units, each dated ca. 1435-1415 Ma (U-Pb zircon + Ar-Ar):
  - Laramie anorthosite complex (Harlan et al., 1994; doi:10.1029/94JB00580; GPMDB 7493)
  - Sherman Granite               (Harlan et al., 1994; GPMDB 7494)
  - Electra Lake Gabbro           (Harlan & Geissman, 1998; doi:10.1029/98JB01350; GPMDB 8342)

Structure (per user direction): one location per unit reporting that unit's
Fisher-mean-of-site-VGP pole, plus an overall "Mean Rocky Mountain intrusions"
location whose pole is the Fisher mean of ALL site VGPs -- the treatment used in
the Nordic / Evans et al. (2021) compilation. All directions are in geographic
(in-situ) coordinates, reproducing the notebook / Nordic export:
    Laramie      -8.1 / 217.1 A95 3.3 N28
    Sherman      -1.1 / 207.6 A95 9.9 N11
    Electra Lake -20.2 / 225.7 A95 3.7 N19   (in-situ; see note)
    Overall      -10.4 / 217.9 A95 5.9 N62
The Laramie complex has a NEGATIVE fold test (99%; Harlan et al., 1994) so its
directions are retained in situ; Laramie and Sherman carry dual polarity; Electra
Lake is single polarity. Harlan & Geissman (1998) preferred a small (~5 deg)
tilt correction for Electra Lake (overlying Cambrian Ignacio Fm), which the
compilation does not apply -- the in-situ Electra Lake pole (-20.2/225.7) differs
from their tilt-corrected pole (-21.1/221.1) by ~1-5 deg. This build follows the
compilation (in situ) for internal consistency of the grand mean.

Audited fixes to the source site table:
  - the source Harlan1998a.csv labels two sites 'ELA14'; the first (37.5447 N,
    between ELlBlk and ELA3) is the mislabeled ELA1 and is renamed accordingly
    (both are diabase-dike sites in the Electra Lake pole, so the pole is
    unchanged; the fix removes the duplicate key);
  - method_codes replaced from the placeholder 'DE-FM': Harlan used AF + thermal
    demagnetization with principal-component line fits (Kirschvink, 1980), plus
    combined line/great-circle analysis for Electra Lake;
  - geologic_classes set to Intrusive (Metamorphic for the melt-zone gneiss);
    geologic_types/lithologies made parallel;
  - location bounding boxes recomputed from the site coordinates (the source
    locations.txt carried a lat_n 51.78 typo and an anomalous lon_w).

No structured field-test CV column is set: the Laramie fold test is negative, the
dual polarity fails a formal reversal test, and no positive baked-contact test
exists (the Electra Lake contact test was inconclusive).

Age: 1430 Ma, bracket 1415-1445, spanning the U-Pb zircon and Ar-Ar ages of the
three units (Laramie ca. 1435-1440; Sherman U-Pb 1431 +/- 6, Ar-Ar 1421 +/- 6;
Electra Lake U-Pb 1435 +/- 2, Ar-Ar 1431 +/- 3).
"""
import os, io, sys
from datetime import date
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

SRC = os.path.join(HERE, 'Rocky_sites_source.txt')
AGE, AGE_LOW, AGE_HIGH = 1430, 1415, 1445
CIT94, CIT98 = '10.1029/94JB00580', '10.1029/98JB01350'

# per-lithology-token CV mapping (keeps geologic_classes/types parallel to lithologies)
LITH_MAP = {
    'Anorthosite': ('Intrusive', 'Pluton'), 'Syenite': ('Intrusive', 'Pluton'),
    'Troctolite': ('Intrusive', 'Pluton'), 'Granite': ('Intrusive', 'Pluton'),
    'Quartz Monzonite': ('Intrusive', 'Pluton'), 'Gabbro': ('Intrusive', 'Pluton'),
    'Diabase': ('Intrusive', 'Volcanic Dike'), 'Gneiss': ('Metamorphic', 'Baked Rock'),
}

# rock-type code -> (unit location, in per-unit pole?, readable description)
UNIT = {
    'An': ('Laramie Anorthosite Complex', True, 'anorthosite'),
    'Sy': ('Laramie Anorthosite Complex', True, 'monzosyenite'),
    'Tr': ('Laramie Anorthosite Complex', True, 'troctolite'),
    'An xen': ('Laramie Anorthosite Complex', False, 'anorthosite-hosted xenolith'),
    'Sy xen': ('Laramie Anorthosite Complex', False, 'syenite-hosted xenolith'),
    'ShGr': ('Sherman Granite', True, 'Sherman granite / quartz monzonite'),
    'Pgb': ('Electra Lake Gabbro', True, 'Electra Lake gabbro'),
    'Pdb': ('Electra Lake Gabbro', True, 'diabase dike cutting the Electra Lake gabbro'),
    'Pgn/Pdb': ('Electra Lake Gabbro', False, 'contact melt-zone gneiss and diabase'),
    'Pgnmz': ('Electra Lake Gabbro', False, 'contact melt-zone gneiss'),
}

# ---- read source, keep in-situ (dir_tilt_correction == 0) rows ---------------
raw = open(SRC, encoding='utf-8').read().splitlines()
d = pd.read_csv(io.StringIO('\n'.join(raw[1:])), sep='\t')
d = d[d['dir_tilt_correction'] == 0].copy()

# fix the duplicate ELA14: the first one (lat ~37.5447) is the mislabeled ELA1
mask = (d['site'] == 'ELA14')
first_idx = d[mask].sort_values('lat').index[0]   # 37.5447 < 37.5475
d.loc[first_idx, 'site'] = 'ELA1'

d['code'] = d['description']
d['location'] = d['code'].map(lambda c: UNIT[c][0])
d['in_pole'] = d['code'].map(lambda c: UNIT[c][1])
d['description'] = d['code'].map(lambda c: UNIT[c][2])
d['method_codes'] = d['citations'].map(
    lambda c: 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM' if c == CIT94
    else 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-BFP:DE-FM')
d['result_type'] = 'i'
d['result_quality'] = 'g'
d['age'] = AGE
d['age_low'] = AGE_LOW
d['age_high'] = AGE_HIGH
d['age_unit'] = 'Ma'


def gclass_gtype(lith):
    toks = str(lith).split(':')
    cls = ':'.join(LITH_MAP.get(t, ('Intrusive', 'Pluton'))[0] for t in toks)
    typ = ':'.join(LITH_MAP.get(t, ('Intrusive', 'Pluton'))[1] for t in toks)
    return pd.Series({'geologic_classes': cls, 'geologic_types': typ})


d[['geologic_classes', 'geologic_types']] = d['lithologies'].apply(gclass_gtype)

# recompute in-situ VGPs from the directions (matches the notebook)
d['vgp_lon'], d['vgp_lat'], d['vgp_dp'], d['vgp_dm'] = pmag.dia_vgp(
    d['dir_dec'], d['dir_inc'], d['dir_alpha95'], d['lat'], d['lon'])

SITE_COLS = ['site', 'location', 'result_type', 'result_quality', 'method_codes',
             'citations', 'geologic_classes', 'geologic_types', 'lithologies',
             'lat', 'lon', 'age', 'age_low', 'age_high', 'age_unit',
             'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k',
             'dir_n_samples', 'description', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm']
sites = d[SITE_COLS].copy()

GEO = dict(continent_ocean='North America', country='United States of America',
           terranes='Laurentia')

UNITS = [
    dict(loc='Laramie Anorthosite Complex', codes=['An', 'Sy', 'Tr'],
         rname='Laramie anorthosite complex ca. 1435 Ma pole', gpmdb='7493',
         cit=CIT94, mc='LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM',
         gclasses='Intrusive', liths='Anorthosite:Syenite:Troctolite',
         state='Wyoming', unify=True,
         desc=("Paleomagnetic pole for the Laramie anorthosite complex, a "
               "Mesoproterozoic massif-type intrusion (anorthosite, mafic "
               "anorthosite, troctolite and associated monzosyenite) of the "
               "southern Laramie Range, Wyoming and Colorado (Harlan et al., 1994). "
               "The pole is the Fisher mean of the in-situ site virtual geomagnetic "
               "poles. A fold test using the igneous layering is negative at the 99% "
               "level, indicating the magnetization postdates doming, so the "
               "directions are retained in geographic coordinates; three sites carry "
               "dual (antipodal) polarity. The characteristic remanence resides in "
               "single- to pseudo-single-domain low-Ti titanomagnetite and was "
               "isolated by alternating-field and thermal demagnetization with "
               "principal-component analysis. The complex is bracketed by U-Pb "
               "zircon ages of ca. 1435-1440 Ma on associated plutons.")),
    dict(loc='Sherman Granite', codes=['ShGr'],
         rname='Sherman Granite ca. 1420 Ma pole', gpmdb='7494',
         cit=CIT94, mc='LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM',
         gclasses='Intrusive', liths='Granite:Quartz Monzonite',
         state='Wyoming', unify=True,
         desc=("Paleomagnetic pole for the Sherman Granite, a composite "
               "Mesoproterozoic batholith (granite to quartz monzonite) of the "
               "southern Laramie Range, Wyoming and Colorado (Harlan et al., 1994). "
               "The pole is the Fisher mean of the in-situ site virtual geomagnetic "
               "poles; the sites are in largely undeformed intrusions and carry dual "
               "(antipodal) polarity. The magnetite-borne characteristic remanence "
               "was isolated by alternating-field and thermal demagnetization with "
               "principal-component analysis. U-Pb zircon dates the granite at "
               "1431 +/- 6 Ma; an Ar-Ar hornblende plateau age of 1421 +/- 6 Ma (2 "
               "sigma) records cooling, giving a magnetization age of about "
               "1415-1430 Ma.")),
    dict(loc='Electra Lake Gabbro', codes=['Pgb', 'Pdb'],
         rname='Electra Lake Gabbro ca. 1433 Ma pole', gpmdb='8342',
         cit=CIT98, mc='LP-DIR-AF:LP-DIR-T:DE-BFL:DE-BFP:DE-FM',
         gclasses='Intrusive', liths='Gabbro:Diabase',
         state='Colorado', unify=True,
         desc=("Paleomagnetic pole for the Electra Lake Gabbro, a small "
               "Mesoproterozoic ophitic-gabbro pluton (with associated diabase "
               "dikes) intruding Paleoproterozoic gneisses of the Needle Mountains, "
               "southwestern Colorado (Harlan & Geissman, 1998). The pole is the "
               "Fisher mean of the in-situ site virtual geomagnetic poles; the unit "
               "carries a single polarity. The characteristic remanence resides in "
               "single- to pseudo-single-domain magnetite and was isolated by "
               "alternating-field and thermal demagnetization with principal-"
               "component and great-circle analysis. A baked-contact test was "
               "attempted but proved inconclusive. U-Pb zircon dates the gabbro at "
               "1435 +/- 2 Ma, with an indistinguishable Ar-Ar biotite plateau age "
               "of 1431 +/- 3 Ma (2 sigma). Directions are reported in geographic "
               "coordinates to match the compilation grand mean; Harlan & Geissman "
               "(1998) preferred a small tilt correction for post-Cambrian tilting.")),
]


def pole_of(sub, unify):
    """Replicates pole_tools.compute_mean_pole(sub, unify_polarity=unify, flip=True)."""
    blk = ipmag.make_di_block(sub['vgp_lon'].tolist(), sub['vgp_lat'].tolist())
    if unify:
        blk = pmag.flip(blk, combine=True)
    blk = ipmag.do_flip(di_block=blk)
    return ipmag.fisher_mean(di_block=blk)


loc_rows = []
for u in UNITS:
    sub = sites[sites['code'].isin(u['codes'])] if 'code' in sites.columns else \
        sites[d['code'].isin(u['codes'])]
    sub = sites.loc[d[d['code'].isin(u['codes'])].index]
    p = pole_of(sub, u['unify'])
    r = {'location': u['loc'], 'location_type': 'Region', 'result_name': u['rname'],
         'result_type': 'a', 'result_quality': 'g',
         'sites': ':'.join(sub['site'].tolist()),
         'method_codes': u['mc'] + ':DE-VGP', 'citations': u['cit'],
         'geologic_classes': u['gclasses'], 'lithologies': u['liths'],
         'lat_s': round(sub['lat'].min(), 3), 'lat_n': round(sub['lat'].max(), 3),
         'lon_w': round(sub['lon'].min(), 3), 'lon_e': round(sub['lon'].max(), 3),
         'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
         'dir_tilt_correction': 0,
         'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
         'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
         'pole_n_sites': int(p['n']),
         'state_province': u['state'], 'description': u['desc']}
    r.update(GEO)
    loc_rows.append(r)

# overall grand-mean pole: ALL site VGPs, flip, no polarity unification (Nordic treatment)
p_all = pole_of(sites, unify=False)
loc_rows.append({
    'location': 'Mean Rocky Mountain intrusions', 'location_type': 'Region',
    'result_name': 'Mean Rocky Mountain intrusions ca. 1430 Ma pole',
    'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(sites['site'].tolist()),
    'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-BFP:DE-FM:DE-VGP',
    'citations': f'{CIT94}:{CIT98}', 'geologic_classes': 'Intrusive',
    'lithologies': 'Anorthosite:Syenite:Troctolite:Granite:Quartz Monzonite:Gabbro:Diabase:Gneiss',
    'lat_s': round(sites['lat'].min(), 3), 'lat_n': round(sites['lat'].max(), 3),
    'lon_w': round(sites['lon'].min(), 3), 'lon_e': round(sites['lon'].max(), 3),
    'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
    'dir_tilt_correction': 0,
    'pole_lat': round(p_all['inc'], 1), 'pole_lon': round(p_all['dec'], 1),
    'pole_alpha95': round(p_all['alpha95'], 1), 'pole_k': round(p_all['k'], 1),
    'pole_n_sites': int(p_all['n']),
    'description': (
        "Grand-mean paleomagnetic pole for the ca. 1430 Ma Rocky Mountain "
        "intrusions of Colorado and Wyoming, combining all site virtual geomagnetic "
        "poles of the Laramie anorthosite complex, the Sherman Granite (Harlan et "
        "al., 1994) and the Electra Lake Gabbro (Harlan & Geissman, 1998). The pole "
        "is the Fisher mean of all site VGPs in geographic (in-situ) coordinates, "
        "following the treatment of these units in the Laurentia apparent-polar-"
        "wander compilation. The three units yield indistinguishable ca. 1435-1415 "
        "Ma U-Pb zircon and Ar-Ar ages."),
    **GEO})

sites = sites.drop(columns=[c for c in ['code'] if c in sites.columns])
locs = pd.DataFrame(loc_rows)


def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Harlan1994_1998_Rocky_Mountain_intrusions_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    f.write('tab delimited\tlocations\n'); locs.to_csv(f, sep='\t', index=False)
    f.write('>>>>>>>>>>\n')
    f.write('tab delimited\tsites\n'); sites.to_csv(f, sep='\t', index=False)

for r in loc_rows:
    print(f"  {r['location']:32s}: {r['pole_lat']}/{r['pole_lon']} A95 {r['pole_alpha95']} "
          f"K {r['pole_k']} N {r['pole_n_sites']}")

from validate_magic_contribution import validate_upload_file
validate_upload_file(combined, tables=['locations', 'sites'])
