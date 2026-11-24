"""Build the MagIC contribution for the Callander alkaline complex of Symons & Chiasson (1991).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Symons, D. T. A., & Chiasson, A. D. (1991). Paleomagnetism of the Callander Complex and the
Cambrian apparent polar wander path for North America. Can. J. Earth Sci. 28, 355-363.
doi:10.1139/e91-032

Dates: Kamo, S. L., Krogh, T. E., & Kumarapeli, P. S. (1995). Age of the Grenville dyke
swarm, Ontario-Quebec: implications for the timing of Iapetan rifting. Can. J. Earth Sci.
32, 273-280. doi:10.1139/e95-022 (p. 278: concordant U-Pb zircon age of 577 +/- 1 Ma on
nepheline syenite of the complex, S. L. Kamo, unpublished data). The K-Ar and Pb-Pb dates
are as reported by Symons & Chiasson (1991). No uncertainty level is stated for any of these
dates, so each is given as age_low/age_high.

Sources and checks
------------------
- Table 1 (p. 357; site means, specimen counts and demagnetization ranges) was transcribed
  from the page image; Table 2 (p. 359) gives the group means and poles, and its notes the
  locality, 46.2N 79.4W.
- The published A mean (82.2/82.7, alpha95 3.1, k 83, N 26) and the syenite, lamprophyre and
  fenite group means are reproduced from the site directions (asserted below).
- Table 1 counts are specimens (Nm measured, Ne end points, Nr remagnetization circles; one
  or two specimens were sliced from each of four or five cores per site). They are stored in
  dir_n_specimens (Ne + Nr), dir_n_specimens_lines, dir_n_specimens_planes and
  dir_n_total_specimens; dir_n_samples is left blank.
- Site 4 is printed with rock code 'SF', which the Table 1 notes do not define. It is one of
  the five syenitic-core sites (with 6, 7, 23 and 24: the only grouping that gives the five
  syenitic and five fenite sites of Table 2) and is taken to be leucocratic syenite (SL).
- Cores were oriented by sun compass (SO-SUN). Two specimens per site were AF demagnetized
  in 13 steps to 160 mT and two thermally in 13 steps to 680 C, and the rest AF demagnetized
  in four steps (LP-DIR-AF:LP-DIR-T); directions are line fits with MAD <= 10 deg (DE-BFL),
  combined with remagnetization circles at the six sites with Nr > 0 (DE-BFP:DE-FM-LP).
- The complex is ca. 3 km across and the paper's Fig. 2 has no graticule, so all sites carry
  the authors' locality coordinate; site VGPs are calculated from the site directions there.

Structure and interpretation
----------------------------
- One location, 'Callander Alkaline Complex', with the pole of the characteristic A
  component (the pole used in the Swanson-Hysell et al. (2026) compilation): the 26 sites
  accepted by the authors (16 lamprophyre dikes, five syenitic-core sites, five fenite
  sites), polarity unified and reported in the northern hemisphere at the antipode of the
  authors' 46.3S 121.4E. Directions are in situ (the complex is untilted).
- Sites 3 (alpha95 25 deg) and 22 (fewer than four reliable specimen directions) were not
  accepted by the authors and have result_quality 'b'.
- Contact test (C+): the 0.7 m wide lamprophyre dike of site 19 intrudes unfenitized
  Grenville granite gneiss (site 29). Block 29A, 0.1 m from the contact, is remagnetized with
  the dike direction; blocks 29B-D, 0.6-4.1 m away, retain the hematite-borne B component of
  the gneiss. The four blocks are host-rock rows ('g', not included in a pole). The authors' B pole
  from the three unbaked blocks (22S 134E), which they correlate with the ca. 980 Ma
  Grenville path, is a single-outcrop result and is not entered as a pole.
- No reversal_test is set: the authors describe their reversal test as 'not robust' with
  three reversed sites.
- Ages are ages of the rocks. Rows of the complex carry 575 Ma (570-580 Ma), the authors'
  estimate from the K-Ar and Pb-Pb dates. The granite gneiss rows carry 975-1590 Ma: younger
  than its 1590 Ma Sm-Nd model age (Dickin and McNutt, 1989) and metamorphosed before the ca.
  975 Ma regional K-Ar biotite cooling (both as reported by Symons & Chiasson, 1991). The
  ages table holds the measured dates on the complex with stated uncertainties.
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

CIT = '10.1139/e91-032'
CIT_KAMO = '10.1139/e95-022'
LOC = 'Callander Alkaline Complex'
LAT, LON = 46.2, 280.6          # 46.2N, 79.4W (Table 2 notes)
CA, CB = 'characteristic', 'B component'
AGE_A = (575, 570, 580)          # Callander Complex (compilation nominal age and bracket)
AGE_GNEISS = (None, 975, 1590)   # Grenville granite gneiss host

# ---- Table 1 (p. 357), transcribed from the page image -------------------------
# site, rock code, Nm, Ne, Nr, AF range (mT), thermal range (C), D, I, a95, k
TAB = [
    ('1', 'L', 9, 8, 0, '20-160', '450-580', 106, 74, 13, 18),
    ('2', 'FI', 9, 8, 0, '10-60', '525-580', 89, 69, 12, 24),
    ('3', 'C', 6, 5, 0, '0-160', '525-630', 65, 55, 25, 10),
    ('4', 'SF', 12, 12, 0, '15-60', '300-580', 6, 85, 10, 21),
    ('5', 'L', 9, 4, 3, '20-60', '525-580', 332, -78, 8, 54),
    ('6', 'SM', 9, 6, 3, '60-120', '500-580', 335, -82, 13, 17),
    ('7', 'T', 6, 5, 0, '10-160', '600-680', 68, 75, 14, 32),
    ('8', 'L', 9, 8, 0, '10-60', '100-580', 119, 83, 12, 21),
    ('9', 'L', 9, 9, 0, '10-80', '200-580', 5, 81, 10, 30),
    ('10', 'L', 11, 8, 3, '10-30', '300-580', 28, -88, 9, 29),
    ('11', 'L', 8, 6, 0, '10-60', '300-580', 36, 83, 13, 29),
    ('12', 'L', 7, 7, 0, '0-160', '200-525', 137, 87, 14, 18),
    ('13', 'L', 8, 8, 0, '20-50', '450-580', 84, 74, 8, 46),
    ('14', 'L', 8, 6, 0, '10-30', '100-580', 125, 78, 6, 117),
    ('15', 'L', 10, 9, 0, '10-80', '500-580', 73, 79, 6, 86),
    ('16', 'L', 10, 9, 0, '10-60', '400-580', 58, 85, 4, 136),
    ('17', 'L', 8, 8, 0, '10-60', '525-580', 95, 80, 10, 31),
    ('18', 'L', 8, 8, 0, '10-50', '300-580', 89, 83, 8, 44),
    ('19', 'L', 12, 12, 0, '10-30', '400-580', 350, 83, 5, 86),
    ('20', 'L', 12, 11, 1, '10-30', '200-570', 116, 65, 7, 36),
    ('21', 'L', 10, 9, 0, '10-60', '400-580', 39, 88, 12, 20),
    ('22', 'FM', 4, 1, 2, '15-25', '400-500', 328, -73, 9, 191),
    ('23', 'SM', 6, 5, 0, '10-30', '450-580', 49, 76, 7, 119),
    ('24', 'T', 6, 4, 1, '15-40', '525-580', 77, 83, 9, 79),
    ('25', 'FM', 8, 6, 0, '10-120', '525-600', 112, 87, 15, 20),
    ('26', 'FO', 8, 6, 0, '20-120', '500-600', 58, 82, 10, 44),
    ('27', 'FO', 7, 6, 0, '20-60', '450-570', 65, 77, 12, 34),
    ('28', 'FO', 6, 6, 0, '10-50', '450-580', 350, 77, 14, 25),
    ('29A', 'A', 4, 4, 0, '5-15', '400-500', 159, 72, 31, 10),
    ('29B', 'B', 3, 3, 0, '10-160', '400-660', 304, -71, 14, 81),
    ('29C', 'C', 5, 5, 0, '15-80', '450-600', 293, -67, 8, 89),
    ('29D', 'D', 3, 3, 0, '50-140', '0-660', 309, -71, 14, 83),
]
NOT_ACCEPTED = {'3': 'alpha95 > 15 deg', '22': 'fewer than four reliable specimen directions'}
BLOCK_DIST = {'A': 0.1, 'B': 0.6, 'C': 1.6, 'D': 4.1}

GEO = {  # code -> (class, type, lithology, description)
    'L': ('Intrusive', 'Volcanic Dike', 'Lamprophyre', 'Lamprophyre dike radiating from the '
          'Callander Complex.'),
    'C': ('Intrusive', 'Volcanic Dike', 'Carbonatite', 'Carbonatite dike of the Callander '
          'Complex.'),
    'T': ('Intrusive', 'Volcanic Dike', 'Trachyte', 'Potassium trachyte dike (syenitic-core '
          'group) of the Callander Complex.'),
    'SM': ('Intrusive', 'Pluton', 'Nepheline syenite', 'Mesocratic nepheline syenite of the '
           'outer shell of the Callander Complex core.'),
    'SF': ('Intrusive', 'Pluton', 'Nepheline syenite', 'Nepheline syenite of the Callander '
           'Complex core (rock code printed SF in Table 1, taken to be leucocratic syenite, '
           'SL; one of the five syenitic-core sites).'),
    'FI': ('Metamorphic', 'Baked Contact', 'Fenite', 'Inner-zone fenite of the aureole '
           '(metasomatized host gneiss with hematite).'),
    'FM': ('Metamorphic', 'Baked Contact', 'Fenite', 'Middle-zone fenite of the aureole.'),
    'FO': ('Metamorphic', 'Baked Contact', 'Fenite', 'Outer-zone fenite of the aureole.'),
}
GNEISS_AGE_TEXT = (' Age of the rock: Grenville gneiss younger than its 1590 Ma Sm-Nd model '
                   'age and metamorphosed before the ca. 975 Ma regional K-Ar biotite cooling '
                   '(as reported by the authors).')


def vgp(dec, inc):
    plon, plat, _, _ = pmag.dia_vgp(dec, inc, 0, LAT, LON)
    return round(float(plat), 1), round(float(plon), 1)


rows = []
for site, code, nm, ne, nr, af, th, dec, inc, a95, k in TAB:
    mc = 'SO-SUN:LP-DIR-AF:LP-DIR-T:DE-BFL:' + ('DE-BFP:DE-FM-LP' if nr else 'DE-FM')
    demag = f' Demagnetization ranges used: {af} mT (AF) and {th} C (thermal).'
    quality, comp, age = 'g', CA, AGE_A
    if site.startswith('29'):
        blk = site[-1]
        age = AGE_GNEISS
        gclass, lith = 'Metamorphic', 'Gneiss'
        gtype = 'Baked Contact' if blk == 'A' else 'Not Specified'
        desc = (f'Unfenitized Grenville granite gneiss, block {blk} of contact-test site 29, '
                f'{BLOCK_DIST[blk]} m from the 0.7 m wide lamprophyre dike of site 19.')
        if blk == 'A':
            desc += (' In the baked zone; remagnetized with the dike direction. This baked '
                     'host rock of the contact test is not included in a pole.')
        else:
            comp = CB
            desc += (' B component of the gneiss, carried by hematite (unblocking above 620 C). '
                     'This unbaked host rock of the contact test is not included in a pole.')
            if blk in 'BC':
                desc += (' A lower-temperature normal component (Holocene viscous remanence or a '
                         'partial overprint from the dike) was removed first.')
        desc += GNEISS_AGE_TEXT
    else:
        gclass, gtype, lith, desc = GEO[code]
        if site in NOT_ACCEPTED:
            quality = 'b'
            desc += f' Not accepted by the authors ({NOT_ACCEPTED[site]}).'
        if code == 'T' and site == '7':
            desc += ' The A component here is carried by hematite (unblocking 600-680 C).'
    plat, plon = vgp(dec, inc)
    rows.append({
        'site': site, 'location': LOC, 'result_type': 'i', 'result_quality': quality,
        'method_codes': mc, 'citations': CIT,
        'geologic_classes': gclass, 'geologic_types': gtype, 'lithologies': lith,
        'lat': LAT, 'lon': LON, 'age': age[0], 'age_low': age[1], 'age_high': age[2],
        'age_unit': 'Ma', 'dir_comp_name': comp, 'dir_tilt_correction': 0,
        'dir_dec': float(dec), 'dir_inc': float(inc), 'dir_alpha95': float(a95),
        'dir_k': float(k), 'dir_n_specimens': ne + nr, 'dir_n_specimens_lines': ne,
        'dir_n_specimens_planes': nr, 'dir_n_total_specimens': nm,
        'vgp_lat': plat, 'vgp_lon': plon,
        'description': desc + demag + (' Remagnetization circles combined with end points.'
                                       if nr else '')})
sites = pd.DataFrame(rows)
for c in ['age', 'age_low', 'age_high']:
    sites[c] = sites[c].astype('Int64')

A = sites[(sites['dir_comp_name'] == CA) & (sites['result_quality'] == 'g')
          & ~sites['site'].str.startswith('29')]
assert len(A) == 26 and A['dir_n_specimens'].sum() == 205


def cart(d, i):
    return np.array(pmag.dir2cart([d, i, 1.0]))


def unified(sub, ref):
    """Invert the directions (and VGPs) more than 90 deg from ref; return dir and pole means."""
    r = cart(*ref)
    m = np.array([np.dot(cart(x.dir_dec, x.dir_inc), r) < 0 for x in sub.itertuples()])
    d = ipmag.fisher_mean(di_block=ipmag.make_di_block(
        list(np.where(m, (sub['dir_dec'] + 180) % 360, sub['dir_dec'])),
        list(np.where(m, -sub['dir_inc'], sub['dir_inc']))))
    p = ipmag.fisher_mean(di_block=ipmag.make_di_block(
        list(np.where(m, (sub['vgp_lon'] + 180) % 360, sub['vgp_lon'])),
        list(np.where(m, -sub['vgp_lat'], sub['vgp_lat']))))
    return d, p, m


# published means (Table 2): N, D, I, alpha95, k
for name, codes, (n, D, I, a95, k) in [
        ('all', None, (26, 82.2, 82.7, 3.1, 83)),
        ('syenites', ('SM', 'SF', 'T'), (5, 68.2, 82.4, 7.8, 97)),
        ('lamprophyres', ('L',), (16, 95.1, 82.9, 4.0, 86)),
        ('fenites', ('FI', 'FM', 'FO'), (5, 61.3, 80.8, 9.9, 61))]:
    ids = {t[0] for t in TAB if codes is None or t[1] in codes}
    d, _, _ = unified(A[A['site'].isin(ids)], (82.2, 82.7))
    got = (int(d['n']), round(d['dec'], 1), round(d['inc'], 1), round(d['alpha95'], 1),
           round(d['k']))
    assert got == (n, D, I, a95, k), (name, got)

GEOG = dict(continent_ocean='North America', country='Canada', state_province='Ontario',
            terranes='Laurentia')
MC_LOC = 'SO-SUN:LP-DIR-AF:LP-DIR-T:DE-BFL:DE-BFP:DE-FM-LP:DE-FM:DE-VGP:ST-C'
DESC_A = (
    "Paleomagnetic pole from the characteristic A magnetization of the Callander alkaline "
    "carbonatite complex, a circular pluton ca. 3 km across beneath Callander Bay at the east "
    "end of Lake Nipissing, Ontario, emplaced into Grenville anorthositic and granitic "
    "gneisses along the Nipissing branch of the Ottawa-Bonnechere graben (Symons and "
    "Chiasson, 1991). The accepted sites are 16 radiating lamprophyre dikes, five "
    "syenitic-core sites (nepheline syenite and potassium trachyte dikes) and five sites in "
    "the fenite aureole; the three groups have indistinguishable means. The authors report "
    "that the complex has not been metamorphosed or tilted since emplacement (untilted "
    "aeromagnetic anomaly; flat-lying Middle Ordovician limestone on the nearby Manitou "
    "Islands complex), so directions are in situ. The A component is carried by magnetite "
    "(and hematite at the trachyte site 7 and in some fenite) with pseudo-single-domain to "
    "multidomain and single-domain grain sizes, and was isolated by alternating-field and "
    "thermal demagnetization with line fits and, at a few sites, remagnetization circles. "
    "The pole is the Fisher mean of the site virtual geomagnetic poles after unifying "
    "polarity, reported in the northern hemisphere at the antipode of the authors' 46.3 S, "
    "121.4 E. A contact test with a lamprophyre dike (site 19) and Grenville granite gneiss "
    "(site 29) is positive: the gneiss 0.1 m from the contact is remagnetized with the dike "
    "direction, whereas 0.6-4.1 m away it retains an older hematite-borne B component, which "
    "the authors correlate with the ca. 980 Ma Grenville segment of the apparent polar wander "
    "path; the gneiss blocks are included as host-rock rows. The authors' reversal test is "
    "positive but, in their words, 'not robust' with three reversed sites. The age of the "
    "complex is estimated by the authors from K-Ar biotite and Pb-Pb mineral isochron dates; "
    "a concordant U-Pb zircon date of 577 +/- 1 Ma on nepheline syenite (S.L. Kamo, "
    "unpublished data) is reported by Kamo et al. (1995).")

d, p, m = unified(A, (82.2, 82.7))
if p['inc'] < 0:
    p['dec'], p['inc'] = (p['dec'] + 180) % 360, -p['inc']
locs = pd.DataFrame([{
    'location': LOC, 'location_type': 'Outcrop',
    'result_name': 'Callander Alkaline Complex A component pole',
    'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(A['site']),
    'method_codes': MC_LOC, 'citations': f'{CIT}:{CIT_KAMO}',
    'geologic_classes': 'Intrusive:Metamorphic',
    'lithologies': 'Lamprophyre:Nepheline syenite:Trachyte:Fenite',
    'lat_s': LAT, 'lat_n': LAT, 'lon_w': LON, 'lon_e': LON,
    'age': AGE_A[0], 'age_low': AGE_A[1], 'age_high': AGE_A[2], 'age_unit': 'Ma',
    'dir_tilt_correction': 0,
    'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
    'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1),
    'dir_n_sites': int(d['n']), 'pole_comp_name': CA,
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
    'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
    'pole_n_sites': int(p['n']), 'pole_reversed_perc': round(100 * m.sum() / len(m), 1),
    'contact_test': 'C+', 'external_database_ids': 'GPMDB[6458]',
    'description': DESC_A, **GEOG}])

# ---- ages table: measured dates with stated uncertainties ------------------------
ages = pd.DataFrame([
    {'location': LOC, 'age': 577, 'age_low': 576, 'age_high': 578, 'age_unit': 'Ma',
     'method_codes': 'GM-UPB', 'citations': CIT_KAMO,
     'description': ('Concordant U-Pb zircon age of a nepheline syenite sample from the '
                     'Callander Complex, 577 +/- 1 Ma (S.L. Kamo, unpublished data, reported '
                     'by Kamo et al., 1995; uncertainty level not stated)')},
    {'location': LOC, 'age': 576, 'age_low': 573, 'age_high': 579, 'age_unit': 'Ma',
     'method_codes': 'GM-PBPB-I', 'citations': CIT,
     'description': ('Pb-Pb isochron age from sphene, brookite, anatase, apatite and rutile '
                     'separates of a lamprophyre dike, 576 +/- 3 Ma (Kamo et al., 1989, as '
                     'presented; the abstract gave 550 +/- 4 Ma), as reported by Symons and '
                     'Chiasson (1991); uncertainty level not stated')},
    {'location': LOC, 'age': 564, 'age_low': 550, 'age_high': 578, 'age_unit': 'Ma',
     'method_codes': 'GM-KAR', 'citations': CIT,
     'description': ('Two identical K-Ar biotite ages of lamprophyre dikes, 564 +/- 14 Ma '
                     '(M. Shafiqullah, in Ferguson and Currie, 1972; current decay constants), '
                     'as reported by Symons and Chiasson (1991); uncertainty level not stated. '
                     'Currie (1976) gave K-Ar biotite ages of 568 Ma (lamprophyre dike) and '
                     '575 Ma (nepheline syenite) without uncertainties')},
])
ages['timescale_eon'] = 'Proterozoic'
ages['timescale_era'] = 'Neoproterozoic'
ages['timescale_period'] = 'Ediacaran'


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
upload = os.path.join(OUT, f'Symons1991_Callander_{stamp}.txt')
with open(upload, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {upload}')
print(f'    {len(sites)} site rows ({len(A)} accepted A sites, {A["dir_n_specimens"].sum()} '
      f'specimens), {len(locs)} pole row, {len(ages)} ages')
r = locs.iloc[0]
print(f"  pole {r['pole_lat']}/{r['pole_lon']} A95 {r['pole_alpha95']} K {r['pole_k']} "
      f"N {r['pole_n_sites']} rev% {r['pole_reversed_perc']}; dir {r['dir_dec']}/{r['dir_inc']} "
      f"a95 {r['dir_alpha95']} k {r['dir_k']}")
print('  published (Table 2): 82.2/82.7 a95 3.1 k 83, pole 46.3S 121.4E dp 5.9 dm 6.1')

# pole of the Swanson-Hysell et al. (2026) compilation (lat, lon, A95, B, N specimens)
EXPECT = (46.5, 300.8, 5.8, 26, 205)
got = (r['pole_lat'], r['pole_lon'], r['pole_alpha95'], int(r['pole_n_sites']),
       int(A['dir_n_specimens'].sum()))
assert got == EXPECT, got
print('-I- pole matches the compilation')

sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True)
