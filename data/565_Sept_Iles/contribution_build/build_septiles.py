"""Build the MagIC contribution for the Sept-Iles layered mafic intrusion of Tanczyk et al. (1987).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Tanczyk, E. I., Lapointe, P., Morris, W. A., & Schmidt, P. W. (1987). A paleomagnetic study
of the layered mafic intrusion at Sept-Iles, Quebec. Can. J. Earth Sci. 24, 1431-1438.
doi:10.1139/e87-135

Dates: Higgins, M. D., & van Breemen, O. (1998). The age of the Sept Iles layered mafic
intrusion, Canada: implications for the late Neoproterozoic/Cambrian history of
southeastern Canada. J. Geol. 106, 421-432. doi:10.1086/516033. Legros, H., Kieffer, M. A.,
Dare, S. A. S., Namur, O., Higgins, M. D., Savard, D., & Glorie, S. (2026). Apatite
geochronology of the Sept-Iles Intrusive Suite (Grenville province, Canada) reveals two
magmatic episodes 55 Myrs apart related to plume arrival and subsequent rifting. Contrib.
Mineral. Petrol. 181, 32. doi:10.1007/s00410-026-02320-1. Also cited: Bono, R. K., &
Tarduno, J. A. (2015). A stable Ediacaran Earth recorded by single silicate crystals of the
ca. 565 Ma Sept-Iles intrusion. Geology 43, 131-134. doi:10.1130/G36247.1

Sources and checks
------------------
- Table 1 (p. 1433; site means of Remanences A and B, N = samples) was transcribed from the
  page image; Table 3 (p. 1435) gives the published means and poles.
- The 10-site Remanence A mean (R 9.732), the 12-site host-rock Remanence B mean (R 11.758)
  and the 4-dyke Remanence B mean are reproduced from the site directions (asserted below).
- Table 1 prints alpha95 = 8 for site 11 (N 6, k 100), which implies 6.7; values are stored
  as printed.
- Sample remanences were calculated 'giving unit weight to treatment steps' (Table 2 notes),
  i.e. as Fisher means of the stable steps (DE-FM). LP-DIR-AF is given where Table 1 lists an
  alternating-field stability range, and to the site-2 dyke, where AF isolated Remanence B in
  some samples; all sites were thermally demagnetized.
- The paper gives no site coordinates. Approximate positions of the seven sampling groups are
  read from Fig. 1b. Its longitude ticks are printed 50 deg 30'W and 50 deg 15'W, which are
  66 deg 30'W and 66 deg 15'W (their spacing equals 15' of longitude at 50 N and matches the
  10 km scale bar); the georeference reproduces the Natural Resources Canada gazetteer
  positions of Ile Grande Basque (50.162N, 66.364W) and Ile Petite Boule (50.166N, 66.300W)
  to within 0.005 deg. Site VGPs are calculated from the site directions at these positions.

Structure and interpretation
----------------------------
- Location 'Sept-Iles Layered Intrusion': the Remanence A pole (component 'characteristic'),
  the pole used in the Swanson-Hysell et al. (2026) compilation, from the 10 gabbro and
  anorthosite sites. Remanence A is the primary thermochemical remanence of the intrusion
  (Tanczyk et al., 1987); single polarity, in situ (subhorizontal layering, no tilt
  correction). It is reported at the southern-hemisphere antipode of the authors' 20N/141E.
- Location 'Sept-Iles diabase dykes': the Remanence B pole of the four diabase dykes that
  cut the intrusion (the dykes' own magnetization; the authors' pole B_D), with the positive
  contact test (C+). Its mean direction is given in the majority (downward) polarity.
- Remanence B of the host rock near the dykes (12 sites; the authors' pole B_C) is the
  remagnetized side of the contact test and is not entered as a pole: its site means are
  rows of component 'Remanence B' that are not included in a pole. Seven localities carry both
  remanences and have one row per component. The transition from A to B towards the dykes
  is an inverse contact test (IC+) on the Remanence A pole.
- Ages are ages of the rocks. Rows of the intrusion carry 565 Ma (561-569 Ma) from the U-Pb
  zircon date of 565 +/- 4 Ma on a granophyric segregation in the Layered Series on the
  mainland (Higgins and van Breemen, 1998), redated by Legros et al. (2026) at 563.9 +/- 1.9
  Ma (concordia) and 564.7 +/- 2.3 Ma (206Pb/238U weighted mean). Legros et al. (2026) show
  that the Layered Series comprises a ca. 620 Ma unit and the ca. 565 Ma Lac au Castor unit
  that holds the dated granophyre; the sampled anorthosite and gabbro of Presqu'ile Marconi
  and the islands are not directly dated. The diabase dykes are undated: they cut the
  intrusion and are absent from the locally exposed Ordovician limestone (Tanczyk et al., 1987),
  so they carry 487-569 Ma (start of the Ordovician to the intrusion; not a measured age). The
  ages table holds the three measured dates.
"""
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))

CIT = '10.1139/e87-135'
CIT_AGE_HvB = '10.1086/516033'
CIT_AGE_L = '10.1007/s00410-026-02320-1'
CIT_BONO = '10.1130/G36247.1'
LOC_INT = 'Sept-Iles Layered Intrusion'
LOC_DYK = 'Sept-Iles diabase dykes'
AGE_INT = (565, 561, 569)        # age of the intrusion (compilation nominal age and bracket)
AGE_DYK = (None, 487, 569)       # diabase dykes: cut the intrusion, absent from the Ordovician limestone
COMP_A, COMP_B = 'characteristic', 'Remanence B'

# ---- sampling-group coordinates read from Fig. 1b (approximate) ---------------
# pixel positions on a 250 dpi render of p. 1432, georeferenced with the (misprinted)
# 66 deg 30'W / 66 deg 15'W ticks at x = 522 / 959 and 50 deg 10'N / 50 deg 05'N at
# y = 1087 / 1310


def _fig1b(x, y):
    lat = 50 + 10 / 60 - (y - 1087) / 44.6 / 60
    lon = 360 - (66.25 + (959 - x) / 29.13 / 60)
    return round(lat, 2), round(lon, 2)


GROUPS = {
    'GB': (_fig1b(765, 1090), 'east shore of Ile Grande Basque'),
    'PB': (_fig1b(862, 1076), 'Ile Petite Boule'),
    'M10': (_fig1b(642, 1092), 'north-east point of Presqu\'ile Marconi'),
    'M11': (_fig1b(604, 1145), 'east shore of Presqu\'ile Marconi'),
    'M12': (_fig1b(585, 1190), 'southern tip of Presqu\'ile Marconi'),
    'M17': (_fig1b(556, 1100), 'north shore of Presqu\'ile Marconi'),
    'M19': (_fig1b(504, 1110), 'north-west shore of Presqu\'ile Marconi'),
}
SITE_GROUP = {**{s: 'GB' for s in range(1, 10)}, 10: 'M10', 11: 'M11', 12: 'M12', 13: 'M12',
              14: 'M12', 15: 'PB', 16: 'PB', 17: 'M17', 18: 'M17', 19: 'M19'}

# ---- Table 1 (p. 1433), transcribed from the page image ------------------------
# Remanence A: site, D, I, N (samples), k, alpha95, afs (mT), T_UB (C)
TAB_A = [
    (1, 320, -34, 13, 61, 5, '120', '560'),
    (3, 336, -40, 9, 61, 7, '75-150', '570'),
    (4, 331, -41, 3, 16, 31, None, '570'),
    (8, 334, -21, 5, 21, 17, '75-150', '570'),
    (10, 342, -24, 3, 57, 16, '75', '540'),
    (11, 326, -39, 6, 100, 8, '75-120', '550'),
    (13, 323, -30, 4, 28, 18, '50-110', None),
    (14, 330, -10, 5, 173, 6, None, '570'),
    (15, 344, -38, 2, None, None, '40', '580'),
    (17, 343, -10, 1, None, None, None, '550'),
]
# Remanence B: site, D, I, N (samples), k, alpha95, T_UB (C)
TAB_B = [
    (2, 8, 80, 8, 32, 10, '540'),
    (4, 256, -79, 4, 22, 20, '500'),
    (5, 174, -88, 2, None, None, '530'),
    (6, 262, 89, 8, 18, 13, '550'),
    (7, 347, 65, 6, 126, 6, '500'),
    (8, 31, 83, 2, None, None, '500'),
    (9, 38, 71, 6, 10, 22, '525'),
    (10, 23, -78, 4, 30, 17, '550'),
    (11, 254, -88, 2, None, None, '540'),
    (12, 28, -83, 8, 32, 10, '550'),
    (13, 263, -88, 2, None, None, '560'),
    (15, 25, 85, 6, 23, 14, '550'),
    (16, 349, 76, 4, 94, 10, '525'),
    (17, 21, 79, 4, 27, 18, '550'),
    (18, 19, 77, 6, 25, 14, '555'),
    (19, 73, -73, 4, 80, 10, '530'),
]
DYKES = {2, 6, 12, 18}
DYKE_CONTACT = {1: '0.25 m from the dyke of site 2', 5: 'within 0.1 m of the dyke of site 6',
                13: 'within 0.05 m of the dyke of site 12'}
SINGLE_POLARITY_B = {7, 9, 10, 17}      # 'Exclusively one polarity ... at sites 7, 9, 10, and 17'
assert sum(r[3] for r in TAB_A) == 51   # Table 3b: 51 samples
HOST_B_EXCLUDED = (' This dyke-related overprint of the host rock, the remagnetized side of '
                   'the contact test, is not included in a pole.')

rows = []


def vgp(dec, inc, lat, lon):
    plon, plat, _, _ = pmag.dia_vgp(dec, inc, 0, lat, lon)
    return round(float(plat), 1), round(float(plon), 1)


def add(site, comp, quality, dec, inc, n, k, a95, mc, gtype, lith, age, desc):
    (lat, lon), place = GROUPS[SITE_GROUP[site]]
    plat, plon = vgp(dec, inc, lat, lon)
    a, lo, hi = age
    rows.append({
        'site': f'SI{site}', 'location': LOC_DYK if site in DYKES else LOC_INT,
        'result_type': 'i', 'result_quality': quality, 'method_codes': mc, 'citations': CIT,
        'geologic_classes': 'Intrusive', 'geologic_types': gtype, 'lithologies': lith,
        'lat': lat, 'lon': lon, 'age': a, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
        'dir_comp_name': comp, 'dir_tilt_correction': 0,
        'dir_dec': float(dec), 'dir_inc': float(inc),
        'dir_alpha95': None if a95 is None else float(a95),
        'dir_k': None if k is None else float(k), 'dir_n_samples': n,
        'vgp_lat': plat, 'vgp_lon': plon,
        'description': desc + f' Site coordinates are approximate, read from Fig. 1b of the '
                              f'paper ({place}).'})


for site, dec, inc, n, k, a95, afs, tub in TAB_A:
    mc = 'LP-DIR-T:LP-DIR-AF:DE-FM' if afs else 'LP-DIR-T:DE-FM'
    d = 'Gabbro or anorthosite of the Sept-Iles intrusion, Remanence A of Tanczyk et al. (1987).'
    if afs:
        d += f' Alternating-field stability {afs} mT.'
    if tub:
        d += f' Unblocking temperature {tub} C.'
    if site in DYKE_CONTACT:
        d += f' Nearest sample {DYKE_CONTACT[site]}.'
    if site == 13:
        d = d.replace('Gabbro or anorthosite', 'Gabbroic anorthosite')
    if n < 3:
        d += ' No site statistics (fewer than three samples).'
    lith = 'Anorthosite' if site == 13 else 'Anorthosite:Gabbro'
    add(site, COMP_A, 'g', dec, inc, n, k, a95, mc, 'Layered Intrusion', lith, AGE_INT, d)

for site, dec, inc, n, k, a95, tub in TAB_B:
    pol = (' All samples of one polarity.' if site in SINGLE_POLARITY_B
           else ' Samples of both polarities; the site mean is expressed in the prevailing '
                'polarity.')
    stats = ' No site statistics (fewer than three samples).' if n < 3 else ''
    if site in DYKES:
        mc = 'LP-DIR-T:LP-DIR-AF:DE-FM' if site == 2 else 'LP-DIR-T:DE-FM'
        d = ('Diabase dyke cutting the Sept-Iles intrusion, Remanence B of Tanczyk et al. '
             f'(1987). Unblocking temperature {tub} C.{pol}{stats}')
        if site == 18:
            d += (' The site includes the dyke and a 10 m wide host-rock contact zone, whose '
                  'remanences are indistinguishable.')
        if site == 2:
            d += ' Alternating-field treatment isolated this remanence in some samples.'
        add(site, COMP_B, 'g', dec, inc, n, k, a95, mc, 'Volcanic Dike', 'Diabase', AGE_DYK, d)
    else:
        d = ('Gabbro or anorthosite of the Sept-Iles intrusion carrying the dyke-related '
             f'Remanence B of Tanczyk et al. (1987). Unblocking temperature {tub} C.{pol}'
             f'{stats}')
        if site in DYKE_CONTACT:
            d += f' Nearest sample {DYKE_CONTACT[site]}.'
        if site == 13:
            d = d.replace('Gabbro or anorthosite', 'Gabbroic anorthosite')
        lith = 'Anorthosite' if site == 13 else 'Anorthosite:Gabbro'
        add(site, COMP_B, 'g', dec, inc, n, k, a95, 'LP-DIR-T:DE-FM', 'Layered Intrusion', lith,
            AGE_INT, d + HOST_B_EXCLUDED)

sites = pd.DataFrame(rows)
sites['_n'] = sites['site'].str[2:].astype(int)
sites = sites.sort_values(['dir_comp_name', '_n'], ascending=[False, True]).drop(columns='_n')
sites = sites.reset_index(drop=True)
for c in ['age', 'age_low', 'age_high', 'dir_n_samples']:
    sites[c] = sites[c].astype('Int64')

A = sites[sites['dir_comp_name'] == COMP_A]
BH = sites[(sites['dir_comp_name'] == COMP_B) & (sites['location'] == LOC_INT)]
BD = sites[(sites['dir_comp_name'] == COMP_B) & (sites['location'] == LOC_DYK)]
assert (len(A), len(BH), len(BD)) == (10, 12, 4)


def unify(decs, incs, ref_inc_sign):
    """Put every direction (or VGP) into the hemisphere of a reference inclination sign."""
    blk = []
    for d, i in zip(decs, incs):
        if np.sign(i) != ref_inc_sign and i != 0:
            d, i = (d + 180) % 360, -i
        blk.append([d, i])
    return blk


def pole_south(sub):
    """Fisher mean of the polarity-unified VGPs, reported in the southern hemisphere."""
    return ipmag.fisher_mean(di_block=unify(sub['vgp_lon'], sub['vgp_lat'], -1))


def dir_mean(sub, inc_sign):
    """Fisher mean of the site directions in one polarity (by inclination sign)."""
    # Both remanences separate cleanly into polarity by inclination sign (A is shallow NW-up;
    # B is within 20 deg of vertical).
    return ipmag.fisher_mean(di_block=unify(sub['dir_dec'], sub['dir_inc'], inc_sign))


def minority_perc(sub, inc_sign):
    return round(100 * (np.sign(sub['dir_inc']) != inc_sign).sum() / len(sub), 1)


# published means (Table 3a): N, D, I, R
for sub, sign, (n, D, I, R) in [(A, -1, (10, 333, -29, 9.732)),
                                (BH, -1, (12, 186, -85, 11.758)),
                                (BD, -1, (4, 188, -86, 3.963))]:
    d = dir_mean(sub, sign)
    # within 0.5 deg of the published direction (the steep dyke mean: 184.8/-86.0 vs 188/-86)
    assert int(d['n']) == n and pmag.angle([d['dec'], d['inc']], [D, I])[0] < 0.5, (n, d)
    assert abs(d['r'] - R) < 0.005, (n, d['r'])

GEOG = dict(continent_ocean='North America', country='Canada', state_province='Quebec',
            terranes='Laurentia')

DESC_A = (
    "Paleomagnetic pole from Remanence A of the Sept-Iles layered mafic intrusion, north "
    "shore of the St. Lawrence River, Quebec (Tanczyk et al., 1987): gabbroic and "
    "anorthositic members sampled on Presqu'ile Marconi and the islands of Grande Basque "
    "and Petite Boule. The remanence is north-west and up, unblocks between 550 and 580 C "
    "with an abrupt intensity drop, resists alternating fields of 75 mT or more, and is "
    "carried by magnetite dust in plagioclase; the authors interpret it as a primary "
    "thermochemical remanence acquired during initial cooling of the intrusion and before "
    "emplacement of the diabase dykes. Igneous layering is subhorizontal and no "
    "post-emplacement tilting was recognized, so directions are in situ. Sample remanences "
    "were calculated giving unit weight to treatment steps and accepted if alpha95 < 12 deg. "
    "The pole is the Fisher mean of the site virtual geomagnetic poles, reported at the "
    "southern-hemisphere antipode of the authors' 20 N, 141 E. The dykes and the host rock "
    "near dyke contacts carry the steep, dual-polarity Remanence B (host-rock site means "
    "included but not part of a pole); the transition from A to B towards the dykes "
    "constitutes an inverse contact test indicating that A predates the dykes, which are "
    "absent from the locally exposed Ordovician limestone. Bono and Tarduno (2015) report that "
    "single plagioclase crystals from the same gabbro-anorthosite unit carry a shallow "
    "dual-polarity remanence resembling Remanence A, and interpret the steep Remanence B as a "
    "secondary overprint carried by multidomain magnetite. "
    "The U-Pb zircon dates are from a granophyric segregation in the Layered Series on the "
    "mainland (Higgins and van Breemen, 1998; Legros et al., 2026). Legros et al. (2026) show "
    "that the Layered Series comprises a ca. 620 Ma unit and the ca. 565 Ma Lac au Castor "
    "unit that holds the dated granophyre; the anorthosite and gabbro of the peninsula and "
    "islands are not directly dated.")
DESC_BD = (
    "Pole from Remanence B of four diabase dykes that cut the Sept-Iles intrusion "
    "(Tanczyk et al., 1987, Table 3, 'Remanence B (dykes)'). The dykes are the last igneous "
    "event in the area and are undated; they are absent from the locally exposed Ordovician "
    "limestone, and field evidence suggests they followed the main "
    "intrusive event not long after, and the authors state that the timing of Remanence B "
    "must remain in question until the dykes are dated. The remanence is steep, of both "
    "polarities, with distributed unblocking temperatures, and is largely removed by "
    "alternating fields. The contact test is positive: the host intrusion is remagnetized to "
    "the dyke direction near the contacts (a remanence statistically identical to that of "
    "the dykes, with pole 61 S, 115 E in the paper) and retains Remanence A away from them. "
    "The pole is the Fisher mean of the polarity-unified site virtual geomagnetic poles, "
    "reported in the southern hemisphere as in the paper; the mean direction is given in the "
    "majority (downward) polarity, whereas the paper lists it with upward inclination. The "
    "age bracket (start of the Ordovician, ca. 487 Ma, to the age of the intrusion) is not a "
    "measured age. "
    "Bono and Tarduno (2015) interpret the steep "
    "direction as a secondary magnetization carried by multidomain magnetite.")


def loc_row(loc, rname, sub, comp, inc_sign, age, mc, liths, desc, cits, **extra):
    p, d = pole_south(sub), dir_mean(sub, inc_sign)
    a, lo, hi = age
    r = {'location': loc, 'location_type': 'Region', 'result_name': rname,
         'result_type': 'a', 'result_quality': 'g',
         'sites': ':'.join(sub['site']), 'method_codes': mc, 'citations': cits,
         'geologic_classes': 'Intrusive', 'lithologies': liths,
         'lat_s': sub['lat'].min(), 'lat_n': sub['lat'].max(),
         'lon_w': sub['lon'].min(), 'lon_e': sub['lon'].max(),
         'age': a, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
         'dir_tilt_correction': 0,
         'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
         'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1),
         'dir_n_sites': int(d['n']), 'pole_comp_name': comp,
         'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
         'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
         'pole_n_sites': int(p['n']), 'pole_reversed_perc': minority_perc(sub, inc_sign),
         'description': desc}
    r.update(extra)
    r.update(GEOG)
    return r


locs = pd.DataFrame([
    loc_row(LOC_INT, 'Sept-Iles Layered Intrusion Remanence A pole', A, COMP_A, -1, AGE_INT,
            'LP-DIR-T:LP-DIR-AF:DE-FM:DE-VGP:ST-C-I', 'Anorthosite:Gabbro', DESC_A,
            f'{CIT}:{CIT_AGE_HvB}:{CIT_AGE_L}:{CIT_BONO}', contact_test='IC+',
            external_database_ids='GPMDB[1752]'),
    loc_row(LOC_DYK, 'Sept-Iles diabase dykes Remanence B pole', BD, COMP_B, +1, AGE_DYK,
            'LP-DIR-T:LP-DIR-AF:DE-FM:DE-VGP:ST-C', 'Diabase', DESC_BD, f'{CIT}:{CIT_BONO}',
            contact_test='C+'),
])
for c in ['age', 'age_low', 'age_high']:
    locs[c] = locs[c].astype('Int64')

# ---- ages table: measured dates only --------------------------------------------
AGE_ROWS = [
    (565, None, 561, 569, 'GM-UPB', CIT_AGE_HvB,
     'U-Pb ID-TIMS zircon, granophyric segregation in gabbro of the Layered Series '
     '(mainland), weighted mean of the 207Pb/206Pb model ages of two concordant single '
     'grains (A 565.6 +/- 1.0, D 564.3 +/- 1.0 Ma, 2 sigma; a third grain was not used), '
     '565 +/- 4 Ma (Higgins and van Breemen, 1998; uncertainty level of the weighted mean '
     'not stated)'),
    (563.9, 0.95, None, None, 'GM-UPB-CC', CIT_AGE_L,
     'U-Pb LA-ICP-MS zircon concordia date (eight analyses) from the same granophyre '
     'zircon concentrate dated by Higgins and van Breemen (1998), 563.9 +/- 1.9 Ma '
     '(2 sigma; MSWD 1.6) (Legros et al., 2026)'),
    (564.7, 1.15, None, None, 'GM-UPB', CIT_AGE_L,
     'U-Pb LA-ICP-MS zircon 206Pb/238U weighted-mean date from the same granophyre zircon '
     'concentrate, 564.7 +/- 2.3 Ma (2 sigma; MSWD 0.79), retained by Legros et al. (2026) '
     'as the crystallization age of the unit (former Mega Cyclic Unit II, their Lac au '
     'Castor unit)'),
]
ages = pd.DataFrame([{
    'location': LOC_INT, 'age': a, 'age_sigma': s, 'age_low': lo, 'age_high': hi,
    'age_unit': 'Ma', 'method_codes': mc, 'citations': cit,
    'timescale_eon': 'Proterozoic', 'timescale_era': 'Neoproterozoic',
    'timescale_period': 'Ediacaran', 'description': desc}
    for a, s, lo, hi, mc, cit, desc in AGE_ROWS])
for c in ['age_low', 'age_high']:
    ages[c] = ages[c].astype('Int64')


def fmt(df):
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
upload = os.path.join(OUT, f'Tanczyk1987_Sept_Iles_{stamp}.txt')
with open(upload, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {upload}')
print(f'    {len(sites)} site rows (A {len(A)}, host-rock B {len(BH)}, dyke B {len(BD)}), '
      f'{len(locs)} pole rows, {len(ages)} ages')
for r in locs.itertuples():
    print(f'  {r.result_name:45s} pole {r.pole_lat}/{r.pole_lon} A95 {r.pole_alpha95} '
          f'K {r.pole_k} N {r.pole_n_sites}; dir {r.dir_dec}/{r.dir_inc} a95 {r.dir_alpha95} '
          f'k {r.dir_k} minority% {r.pole_reversed_perc}')
print('  Table 3a (published): A 333/-29 k 34 a95 8, pole 20N 141E; B dykes 188/-86 R 3.963 '
      'k 81 a95 10, pole 59S 116E')

# poles (lat, lon, A95, B, N samples); the Remanence A pole is the pole of the
# Swanson-Hysell et al. (2026) compilation
EXPECT = {'Sept-Iles Layered Intrusion Remanence A pole': (-19.7, 321.3, 6.5, 10, 51),
          'Sept-Iles diabase dykes Remanence B pole': (-57.8, 114.6, 20.7, 4, 30)}
for r, sub in zip(locs.itertuples(), (A, BD)):
    got = (r.pole_lat, r.pole_lon, r.pole_alpha95, r.pole_n_sites,
           int(sub['dir_n_samples'].sum()))
    assert got == EXPECT[r.result_name], (r.result_name, got)
print('-I- poles match the expected values')

sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True)
