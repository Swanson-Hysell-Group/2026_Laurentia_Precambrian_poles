"""Build the MagIC contribution for the Catoctin volcanic province of Meert et al. (1994).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Meert, J. G., Van der Voo, R., & Payne, T. W. (1994). Paleomagnetism of the Catoctin volcanic
province: A new Vendian-Cambrian apparent polar wander path for North America. Journal of
Geophysical Research 99(B3), 4625-4641. doi:10.1029/93JB01723

Dates: Aleinikoff, J. N., Zartman, R. E., Walter, M., Rankin, D. W., Lyttle, P. T., & Burton,
W. C. (1995). U-Pb ages of metarhyolites of the Catoctin and Mount Rogers Formations, central
and southern Appalachians: evidence for two pulses of Iapetan rifting. American Journal of
Science 295, 428-454. doi:10.2475/ajs.295.4.428 (pp. 442-444; the intercept uncertainty level
is not stated, so the ages table gives the intercept ranges as age_low/age_high).

Sources and checks
------------------
- Table 1 (p. 4630, printed rotated; B and C components) and Table 2 (p. 4633; A component)
  were transcribed from the page images: site, n/N (samples used / collected), in-situ and
  tilt-corrected directions, k and alpha95. Tilt-corrected directions are not given for the C
  component, which fails the fold test.
- Site VGPs are calculated from the site directions. The VGP latitudes in Tables 1 and 2 agree
  with these, but the tabulated VGP longitudes, and those of the three mean poles, are ca. 10
  deg east of the values calculated from the directions at the study area (median 10.2 deg
  over the 47 tabulated VGPs; as if calculated at ca. 68 W rather than ca. 78 W). For example,
  the A-component mean direction 68/+84 gives 42.0S 116.7E (Table 2: 43S 128E). The poles here
  are therefore ca. 10 deg west of the published ones, in agreement with the corrected
  longitudes of the Global Paleomagnetic Database (results 7474-7476) and with McCausland et
  al. (2007, 2011; Catoctin A 42.0N 296.7E).
- Table 1 prints the C-component VGP of site 22 as 21N 168E, as for site 29; the VGP is
  calculated from the direction (143/+58).
- The three published mean directions are reproduced from the site directions (printed by
  the build); the A-component pole is asserted against the compilation value below.
- The paper gives no site coordinates. Approximate positions of the labelled sampling groups
  are read from Fig. 1 (79W/78W and 39N/38N ticks). The map is schematic (Charlottesville
  plots within 0.06 deg, the Pennsylvania-Maryland border ca. 0.3 deg too far south), so the
  coordinates are good to a few tenths of a degree; sites within a group share one coordinate.
  A single study-area coordinate (38.5N, 78.2W) instead gives 42.5N 297.0E for the A pole.

Structure and interpretation
----------------------------
- One location with the three component poles reported by the authors: the A component
  (dir_comp_name 'characteristic'; magnetite, dual polarity, tilt corrected, the pole used in
  the Swanson-Hysell et al. (2026) compilation), the B component (hematite, Late Cambrian
  remagnetization, tilt corrected) and the C component (Taconic remagnetization, in situ).
  The B and C poles are published results of the paper and are kept as pole rows.
- Site rows: A and B components in situ and tilt corrected, C component in situ. A sites 6
  and 16 have one sample each and are 'b' rows (the authors' preferred Mean-B omits them).
- Each pole is the Fisher mean of the site VGPs after polarity unification: A in the
  northern hemisphere (as in the compilation), B in the southern and C in the northern
  hemisphere, as in the paper. Mean directions are in the majority (or the paper's) polarity.
- Lithology follows the paper: basalt flows (site 2 porphyritic basalt, where the A
  component is carried by hematite); site 8, a 4 m feeder dike, and site 28, an 8 m sill,
  intrude the Pedlar Granite and their A means include baked granite samples (positive
  contact test); site 10 is a feeder dike or sill (Fig. 1 caption).
- Field tests recorded: A fold_test F+ (McFadden, 1990 outlier test), reversal_test Rc
  (McFadden & McElhinny, 1990, tilt corrected) and contact_test C+; B fold_test F+. The C
  fold test is negative (in-situ hypothesis preferred) and is not recorded.
- Ages: the A rows carry the compilation age, 568 Ma (555-577 Ma), spanning the two U-Pb
  zircon dates on Catoctin felsic rocks; B and C carry the authors' interpreted
  remagnetization brackets (505-515 Ma and 450-470 Ma), which are not measured ages. The
  ages table holds the two U-Pb dates.
"""
import os
import sys
from datetime import date
import numpy as np
import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

CIT = '10.1029/93JB01723'
CIT_AGE = '10.2475/ajs.295.4.428'
LOC = 'Catoctin volcanic province'
CA, CB, CC = 'characteristic', 'B component', 'C component'
AGES = {CA: (568, 555, 577), CB: (None, 505, 515), CC: (None, 450, 470)}

# ---- sampling-group coordinates read from Fig. 1 (approximate) -----------------
# pixel positions on a 200 dpi render of p. 4626 (map crop starting at 40% of the
# page height); ticks: 79W x = 653, 78W x = 979, 39N y = 250, 38N y = 763


def _fig1(x, y):
    return round(39 - (y - 250) / 513, 2), round(360 - (79 - (x - 653) / 326), 2)


GROUP_XY = {'1-6': (998, 295), '7-12': (960, 338), '13-14': (923, 373), '15': (904, 395),
            '16': (856, 445), '17-18': (845, 510), '19-21': (799, 587),
            '22-24': (745, 664), '25': (1064, 316), '26': (1051, 360), '28': (1023, 469),
            '29': (929, 636), '30': (890, 676), '31': (795, 768)}


def group_of(site):
    s = int(site)
    for g in GROUP_XY:
        lo, hi = (g.split('-') + [g])[:2]
        if int(lo) <= s <= int(hi):
            return g
    raise KeyError(site)


def coords(site):
    return _fig1(*GROUP_XY[group_of(site)])


# ---- Table 2 (p. 4633): A component ---------------------------------------------
# site, n, N, Ds, Is, Dc, Ic, k, a95   (published VGPs not used, see docstring)
TAB_A = [
    (2, 3, 5, 191, -81, 256, -74, 504, 6),
    (3, 6, 6, 289, 52, 281, 80, 154, 6),
    (6, 1, 6, 41, -79, 133, -80, None, None),
    (8, 6, 9, 164, 89, 164, 89, 43, 8),
    (9, 5, 6, 55, 80, 86, 76, 334, 5),
    (15, 4, 6, 194, -76, 203, -86, 97, 9),
    (16, 1, 5, 297, 82, 159, 82, None, None),
    (28, 4, 5, 189, -75, 248, -75, 730, 4),
]
# ---- Table 1 (p. 4630): B component ----------------------------------------------
TAB_B = [
    (1, 5, 6, 100, 43, 102, 28, 87, 8),
    (4, 6, 6, 103, 29, 94, 31, 13, 19),
    (12, 5, 5, 70, 29, 79, 18, 65, 12),
    (17, 4, 5, 107, -4, 106, -16, 20, 25),
    (18, 6, 6, 251, -8, 252, -6, 100, 8),
    (20, 4, 5, 102, -3, 101, 16, 20, 22),
    (25, 7, 7, 51, 49, 77, 20, 27, 10),
    (26, 5, 5, 117, 64, 104, 35, 26, 15),
    (30, 8, 8, 281, -30, 279, -5, 21, 12),
]
# ---- Table 1 (p. 4630): C component (in situ only) ------------------------------
TAB_C = [
    (1, 3, 6, 150, 72, 55, 13), (4, 3, 6, 152, 34, 24, 25), (6, 4, 6, 111, 55, 33, 17),
    (10, 3, 5, 150, 31, 70, 10), (11, 8, 8, 164, 28, 51, 7), (13, 4, 8, 131, 59, 53, 24),
    (16, 3, 5, 137, 55, 36, 21), (19, 4, 5, 157, 29, 154, 9), (21, 3, 3, 146, 30, 21, 26),
    (22, 17, 17, 143, 58, 27, 7), (23, 8, 8, 155, 44, 50, 7), (24, 6, 6, 159, 45, 30, 14),
    (29, 6, 6, 129, 21, 58, 7), (30, 6, 8, 144, 53, 17, 16), (31, 6, 9, 157, 37, 24, 14),
]
assert sum(r[1] for r in TAB_A if r[1] > 1) == 28 and sum(r[2] for r in TAB_A if r[1] > 1) == 37
assert sum(r[1] for r in TAB_B) == 50 and sum(r[1] for r in TAB_C) == 84

# lithology by site: Fig. 1 caption 'Sites 5, 8, 10 and 28 were Catoctin feeder dikes or
# sills'; p. 4629-4631 (site 2 porphyritic basalt; sites 8 and 28 intrude Pedlar Granite)
GEO_BASALT = ('Extrusive', 'Lava Flow', 'Basalt', 'Catoctin Formation basalt flow.')
GEO = {
    2: ('Extrusive', 'Lava Flow', 'Porphyritic Basalt', 'Catoctin Formation porphyritic '
        'basalt flow; the A component here is carried by hematite.'),
    8: ('Intrusive', 'Volcanic Dike:Baked Contact', 'Diabase:Granite',
        'Catoctin feeder dike, 4 m wide, intruding the Mesoproterozoic Pedlar Granite; the '
        'site mean includes samples of the dike and of baked granite at the contact (6 cm '
        'from the contact). Granite 2 m and 5 km from the dike has shallow to intermediate, '
        'poorly grouped directions (positive contact test).'),
    10: ('Intrusive', 'Volcanic Dike', 'Diabase', 'Catoctin feeder dike or sill.'),
    28: ('Intrusive', 'Sill:Baked Contact', 'Diabase:Granite',
         'Catoctin sill, 8 m wide, intruding the Mesoproterozoic Pedlar Granite; the site '
         'mean includes samples of the sill and of baked granite (positive contact test).'),
}
MC_SITE = 'SO-MAG:LP-DIR-T:LP-DIR-AF:DE-BFL:DE-FM'

rows = []


def vgp(dec, inc, lat, lon):
    plon, plat, _, _ = pmag.dia_vgp(dec, inc, 0, lat, lon)
    return round(float(plat), 1), round(float(plon), 1)


def add(site, comp, tc, n, N, dec, inc, k, a95, extra='', quality='g'):
    lat, lon = coords(site)
    gclass, gtype, lith, desc = GEO.get(site, GEO_BASALT)
    age, lo, hi = AGES[comp]
    plat, plon = vgp(dec, inc, lat, lon)
    if tc == 100:
        extra += ' Tilt-corrected direction.'
    rows.append({
        'site': f'CT{site}', 'location': LOC, 'result_type': 'i', 'result_quality': quality,
        'method_codes': MC_SITE, 'citations': CIT,
        'geologic_classes': gclass, 'geologic_types': gtype, 'lithologies': lith,
        'lat': lat, 'lon': lon, 'age': age, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
        'dir_comp_name': comp, 'dir_tilt_correction': tc,
        'dir_dec': float(dec), 'dir_inc': float(inc),
        'dir_alpha95': None if a95 is None else float(a95),
        'dir_k': None if k is None else float(k),
        'dir_n_samples': n, 'dir_n_total_samples': N,
        'vgp_lat': plat, 'vgp_lon': plon,
        'description': (f'{desc}{extra} Site coordinates are approximate (sampling group '
                        f'{group_of(site)}), read from Fig. 1 of the paper; the VGP is '
                        'calculated from the tabulated site direction rather than taken '
                        'from the VGP tabulated in the paper.')})


COMP_TXT = {CA: ' A component of Meert et al. (1994), magnetite (unblocking 500-585 C).',
            CB: ' B component of Meert et al. (1994), hematite (unblocking 640-690 C); a '
                'remagnetization interpreted by the authors as Late Cambrian.',
            CC: ' C component of Meert et al. (1994), magnetite (unblocking 400-580 C); a '
                'remagnetization interpreted by the authors as Taconic.'}

for site, n, N, ds, is_, dc, ic, k, a95 in TAB_A:
    q, ex = ('b', ' [excluded: one sample only; omitted from the authors\' preferred mean]') \
        if n == 1 else ('g', '')
    add(site, CA, 0, n, N, ds, is_, k, a95, COMP_TXT[CA] + ex, q)
    add(site, CA, 100, n, N, dc, ic, k, a95, COMP_TXT[CA] + ex, q)
for site, n, N, ds, is_, dc, ic, k, a95 in TAB_B:
    add(site, CB, 0, n, N, ds, is_, k, a95, COMP_TXT[CB])
    add(site, CB, 100, n, N, dc, ic, k, a95, COMP_TXT[CB])
for site, n, N, ds, is_, k, a95 in TAB_C:
    add(site, CC, 0, n, N, ds, is_, k, a95, COMP_TXT[CC])

sites = pd.DataFrame(rows)
sites['age'] = sites['age'].astype('Int64')


def sel(comp, tc):
    return sites[(sites['dir_comp_name'] == comp) & (sites['dir_tilt_correction'] == tc) &
                 (sites['result_quality'] == 'g')]


A_tc, B_tc, C_geo = sel(CA, 100), sel(CB, 100), sel(CC, 0)
assert (len(A_tc), len(B_tc), len(C_geo)) == (6, 9, 15)


def cart(d, i):
    return np.array(pmag.dir2cart([d, i, 1.0]))


def minority_mask(sub, ref):
    """True where a site direction is > 90 deg from the reference (the other polarity)."""
    r = cart(*ref)
    return np.array([np.dot(cart(x.dir_dec, x.dir_inc), r) < 0 for x in sub.itertuples()])


def dir_mean(sub, ref):
    """Fisher mean of the site directions in the polarity of the reference direction."""
    m = minority_mask(sub, ref)
    decs = np.where(m, (sub['dir_dec'] + 180) % 360, sub['dir_dec'])
    incs = np.where(m, -sub['dir_inc'], sub['dir_inc'])
    return ipmag.fisher_mean(di_block=ipmag.make_di_block(list(decs), list(incs))), m


def pole_mean(sub, m, hemisphere):
    """Fisher mean of the site VGPs after inverting the minority-polarity sites, reported in
    the requested hemisphere (+1 north, -1 south)."""
    lons = np.where(m, (sub['vgp_lon'] + 180) % 360, sub['vgp_lon'])
    lats = np.where(m, -sub['vgp_lat'], sub['vgp_lat'])
    p = ipmag.fisher_mean(di_block=ipmag.make_di_block(list(lons), list(lats)))
    if np.sign(p['inc']) != hemisphere:
        p['dec'], p['inc'] = (p['dec'] + 180) % 360, -p['inc']
    return p


GEOG = dict(continent_ocean='North America', country='United States of America',
            state_province='Virginia', terranes='Laurentia')
LOC_MC = 'SO-MAG:LP-DIR-T:LP-DIR-AF:DE-BFL:DE-BFP:DE-FM-LP:DE-FM:DE-VGP'
LONG_NOTE = (' The site virtual geomagnetic poles (VGPs) in this contribution are calculated '
             'from the site-mean directions tabulated in Tables 1 and 2 of Meert et al. '
             '(1994). The VGP longitudes tabulated in those tables, and the longitudes of '
             'the tabulated mean poles, are ca. 10 degrees east of the values calculated '
             'from the tabulated directions at the study area; the tabulated VGP latitudes '
             'agree with the calculated values.')
DESC = {
    CA: ("Paleomagnetic pole from the A component of the Catoctin volcanic province, Blue "
         "Ridge anticlinorium, central Virginia (Meert et al., 1994): subaerial basalt flows "
         "and feeder dikes and sills of the Catoctin Formation, which records rifting of the "
         "eastern Laurentian margin during opening of the Iapetus Ocean. The A component is "
         "carried by magnetite (unblocking 500-585 C), is steep and of dual polarity, and "
         "was isolated by thermal (and pilot alternating-field) demagnetization with "
         "principal-component analysis. The pole is the Fisher mean of the tilt-corrected "
         "site virtual geomagnetic poles of the six sites with more than one sample (the "
         "authors' preferred Mean-B), after polarity unification, reported in the northern "
         "hemisphere. The A component passes an outlier fold test (McFadden, 1990) and a "
         "class C reversal test in tilt-corrected coordinates, and a baked-contact test is "
         "positive where a feeder dike (site 8) and a sill (site 28) intrude the Pedlar "
         "Granite; the authors consider that it 'may well be a primary remanence'. Two "
         "younger components (B, Late Cambrian; C, Taconic) are reported as separate "
         "results. The U-Pb zircon dates are from a felsic dike at Point of Rocks and a "
         "metarhyolite at South Mountain, Pennsylvania, north of the sampled area."
         + LONG_NOTE),
    CB: ("Pole from the B component of the Catoctin volcanic province (Meert et al., 1994): "
         "a dual-polarity, hematite-carried magnetization (unblocking 640-690 C) that "
         "passes a fold test at 99% confidence and is interpreted by the authors as a "
         "remagnetization during a Late Cambrian (505-515 Ma) orogenic event in the central "
         "Appalachians. The pole is the Fisher mean of the tilt-corrected site virtual "
         "geomagnetic poles after polarity unification, reported in the southern hemisphere "
         "as in the paper; the mean direction is in the majority (easterly, downward) "
         "polarity. The age bracket is the authors' interpretation, not a measured age."
         + LONG_NOTE),
    CC: ("Pole from the C component of the Catoctin volcanic province (Meert et al., 1994): "
         "a single-polarity, magnetite-carried magnetization (unblocking 400-580 C) in "
         "basalts of higher metamorphic grade and stronger deformation. It fails the fold "
         "test (McFadden, 1990; in-situ hypothesis preferred), so the pole is the Fisher "
         "mean of the in-situ site virtual geomagnetic poles, reported in the northern "
         "hemisphere as in the paper. The authors interpret it as a Taconic "
         "remagnetization, correlated with Rb-Sr resetting ages of 450-470 Ma; the age "
         "bracket is not a measured age of the magnetization." + LONG_NOTE),
}
LOCDEF = [
    (CA, A_tc, (68, 84), +1, 'Catoctin volcanic province A component pole', 100,
     dict(fold_test='F+', reversal_test='Rc', contact_test='C+'), ':ST-F:ST-R:ST-C',
     'Basalt:Porphyritic Basalt:Diabase:Granite', 'Extrusive:Intrusive'),
    (CB, B_tc, (92, 17), -1, 'Catoctin volcanic province B component (Late Cambrian '
     'remagnetization) pole', 100, dict(fold_test='F+'), ':ST-F', 'Basalt', 'Extrusive'),
    (CC, C_geo, (147, 44), +1, 'Catoctin volcanic province C component (Taconic '
     'remagnetization) pole', 0, {}, '', 'Basalt:Diabase', 'Extrusive:Intrusive'),
]
GPMDB = {CA: 7474, CB: 7475, CC: 7476}
loc_rows = []
poles = {}
for comp, sub, ref, hemi, rname, tc, tests, st, liths, classes in LOCDEF:
    d, m = dir_mean(sub, ref)
    p = pole_mean(sub, m, hemi)
    poles[comp] = (p, d, sub)
    age, lo, hi = AGES[comp]
    r = {'location': LOC, 'location_type': 'Region', 'result_name': rname,
         'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(sub['site']),
         'method_codes': LOC_MC + st, 'citations': f'{CIT}:{CIT_AGE}' if comp == CA else CIT,
         'geologic_classes': classes, 'lithologies': liths,
         'lat_s': sub['lat'].min(), 'lat_n': sub['lat'].max(),
         'lon_w': sub['lon'].min(), 'lon_e': sub['lon'].max(),
         'age': age, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
         'dir_tilt_correction': tc,
         'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
         'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1),
         'dir_n_sites': int(d['n']), 'pole_comp_name': comp,
         'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
         'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
         'pole_n_sites': int(p['n']), 'pole_reversed_perc': round(100 * m.sum() / len(m), 1),
         'external_database_ids': f'GPMDB[{GPMDB[comp]}]', 'description': DESC[comp]}
    r.update(tests)
    r.update(GEOG)
    loc_rows.append(r)
locs = pd.DataFrame(loc_rows)
locs['age'] = locs['age'].astype('Int64')

# ---- ages table: measured dates (Aleinikoff et al., 1995) -------------------------
ages = pd.DataFrame([
    {'location': LOC, 'age': 572, 'age_low': 567, 'age_high': 577, 'age_unit': 'Ma',
     'method_codes': 'GM-UPB-CC-T0', 'citations': CIT_AGE,
     'timescale_eon': 'Proterozoic', 'timescale_era': 'Neoproterozoic',
     'timescale_period': 'Ediacaran',
     'description': ('U-Pb zircon discordia upper-intercept age of five size fractions, '
                     '571.5 +/- 4.7 Ma (MSWD 0.41), quoted as 572 +/- 5 Ma, from a 2 m wide '
                     'metarhyolite dike cutting Mesoproterozoic granitic gneiss and Late '
                     'Proterozoic metadiabase near Point of Rocks, south side of the Potomac '
                     'River (sample BR-38-90), thought to be a feeder to felsic volcanic rocks '
                     'of the Catoctin Formation (Aleinikoff et al., 1995; isotope-ratio errors '
                     '2 sigma, intercept uncertainty level not stated)')},
    {'location': LOC, 'age': 564, 'age_low': 555, 'age_high': 573, 'age_unit': 'Ma',
     'method_codes': 'GM-UPB-CC-T0', 'citations': CIT_AGE,
     'timescale_eon': 'Proterozoic', 'timescale_era': 'Neoproterozoic',
     'timescale_period': 'Ediacaran',
     'description': ('U-Pb zircon discordia upper-intercept age of four fractions, 564.3 +/- '
                     '9.3 Ma (MSWD 0.047), quoted as 564 +/- 9 Ma, from metarhyolite of the '
                     'Catoctin Formation at South Mountain, Pennsylvania (sample SM-10-89); '
                     'the weighted mean 207Pb/206Pb age of the same fractions is 560 +/- 3 Ma '
                     '(Aleinikoff et al., 1995; intercept uncertainty level not stated)')},
])


def fmt(df):
    """Integers without '.0', floats without trailing digits, blanks for missing values."""
    df = df.dropna(axis=1, how='all').copy()
    for c in df.columns:
        if df[c].dtype.kind == 'f':
            df[c] = df[c].map(lambda x: '' if pd.isna(x) else f'{x:.10g}')
    return df.astype(object).where(df.notna(), '')


def table_text(df, kind):
    return f'tab delimited\t{kind}\n' + fmt(df).to_csv(sep='\t', index=False)


TABLES = (('locations', locs), ('sites', sites), ('ages', ages))
for kind, df in TABLES:
    with open(os.path.join(OUT, f'{kind}.txt'), 'w', encoding='utf-8') as f:
        f.write(table_text(df, kind))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Meert1994_Catoctin_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {combined}')

print(f'-I- {len(sites)} site rows')
for r in loc_rows:
    print(f"  {r['pole_comp_name']:15s} tc{r['dir_tilt_correction']:3d}: pole {r['pole_lat']}/"
          f"{r['pole_lon']} A95 {r['pole_alpha95']} K {r['pole_k']} N {r['pole_n_sites']} "
          f"rev% {r['pole_reversed_perc']}; dir {r['dir_dec']}/{r['dir_inc']} "
          f"a95 {r['dir_alpha95']} k {r['dir_k']}")
print('  published: A 68/+84 k 59 a95 9, pole 43S 128E; B 92/+17 k 16 a95 13, pole 4S 193E; '
      'C 147/+44 k 21 a95 9, pole 19N 143E')

# published mean directions (Tables 1 and 2) are reproduced from the site directions
for comp, (D, I) in {CA: (68, 84), CB: (92, 17), CC: (147, 44)}.items():
    d = poles[comp][1]
    assert pmag.angle([d['dec'], d['inc']], [D, I])[0] < 1.0, comp  # integer site values

# pole of the Swanson-Hysell et al. (2026) compilation (lat, lon, A95, B, N)
EXPECT = (42.8, 297.2, 16.7, 6, 28)
p, d, sub = poles[CA]
got = (round(p['inc'], 1), round(p['dec'], 1), round(p['alpha95'], 1), int(p['n']),
       int(sub['dir_n_samples'].sum()))
assert got == EXPECT, got
print('-I- the A-component pole matches the compilation')

from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(combined, tables=['locations', 'sites', 'ages'], server=True)
