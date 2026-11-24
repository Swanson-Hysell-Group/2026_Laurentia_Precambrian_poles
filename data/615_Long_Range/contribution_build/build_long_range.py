"""Build the MagIC contribution for the Long Range dykes of Murthy et al. (1992).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Murthy, G. S., Gower, C., Tubrett, M., & Patzold, R. (1992). Paleomagnetism of Eocambrian Long
Range dykes and Double Mer Formation from Labrador, Canada. Canadian Journal of Earth Sciences
29, 1224-1234. doi:10.1139/e92-098

Dates: Kamo, S. L., & Gower, C. F. (1994). Note: U-Pb baddeleyite dating clarifies age of
characteristic paleomagnetic remanence of Long Range dykes, southeastern Labrador. Atlantic
Geology 30, 259-262. doi:10.4138/2133. Kamo, S. L., Gower, C. F., & Krogh, T. E. (1989).
Birthdate for the Iapetus Ocean? A precise U-Pb zircon and baddeleyite age for the Long Range
dikes, southeast Labrador. Geology 17, 602-605 (its DOI contains a colon, which would split
MagIC's colon-delimited citations list, so it is named in descriptions only). K-Ar ages:
Grasty et al. (1969, doi:10.1139/e69-030) and Wanless et al. (1970, GSC Paper 69-2A), as cited
by Murthy et al. (1992).

Sources and checks
------------------
- Table 1 (p. 1229) was transcribed from the page image (the text layer of the scan is
  unreliable): 17 sampling sites in six dykes, with N (samples collected, excluding the
  host-rock samples of the contact test), dyke trend and width, D, I, k, alpha95, n1
  (specimens with stable end points) and n2 (specimens from great circles, vector
  subtraction and other procedures), and a mean for each dyke.
- Counts: N -> dir_n_total_samples, n1 -> dir_n_specimens_lines, n2 ->
  dir_n_specimens_planes, n1+n2 -> dir_n_specimens. The number of samples used is not
  given, so dir_n_samples is blank.
- Site k/alpha95 are stored as printed (the site means combine lines and great circles).
  The printed dyke means are reproduced from the site means (asserted below; dyke 6:
  124.4/67.5 calculated, 124.6/67.5 printed, no k or alpha95 printed for this two-site
  mean), and the three-dyke mean 124.8/55.5 (k 48, alpha95 18.0) from the dyke means.
- The paper gives no site coordinates. PR84-191/199/206 are at the U-Pb sample locality of
  Kamo & Gower (1994; 53 deg 46.2' N, 56 deg 37.1' W), which they state corresponds to these
  sites. The other Sandwich Bay sites are read from the Fig. 2 sketch map (corners labelled
  54 deg 00' N, 57 deg 35' W and 53 deg 22' N, 56 deg 14' W; reading Fig. 2 at the Kamo &
  Gower locality agrees within 1 km), and dyke 6 (PR86-008, PR86-015) from Fig. 1 by
  reference to the Fig. 2 frame drawn on it (ca. 5 km). Dyke-mean VGPs are calculated at the
  mean coordinates of each dyke's sites. The site coordinates are mapped by
  plot_site_coordinates.py (Long_Range_site_coordinates.png), which can also plot them on
  Fig. 2 of the paper as a check. With one reference locality (53.7 N, 303.3 E) for every
  dyke instead, the pole is -12.4/165.4, A95 20.8.
- The paper also reports paleomagnetic results from the Double Mer Formation red beds
  (Table 2), which are not part of this contribution.

Structure and interpretation
----------------------------
- One location. The sites table holds 17 sampling-site rows (result_type 'i') and six
  dyke-mean rows (result_type 'a', 'Dyke 1' to 'Dyke 6'); each dyke is one cooling unit and
  the dyke means are the averaging unit of the pole. Dykes 3 and 5 were sampled at one site
  each, so their dyke-mean rows repeat the site mean.
- Dykes 2, 4 and 6 share the southeast-and-down direction that the authors interpret as
  primary ('g'). Dykes 1, 3 and 5 carry the 'anomalous' directions that the authors
  interpret as disturbed and exclude ('b').
- The pole is the Fisher mean of the three primary dyke-mean VGPs, reported in the southern
  hemisphere as in the paper; it is the pole used in the Swanson-Hysell et al. (2026)
  compilation. Directions are in situ (near-vertical dykes).
- Methods (p. 1226): samples oriented by solar and magnetic compass (SO-SM); AF and thermal
  demagnetization (AF for most samples); PCA (MAD < 10 deg) for end-point specimens
  (DE-BFL) and remagnetization circles (Halls, 1976, 1978) or vector subtraction (DE-BFP);
  Fisher site means (DE-FM).
- Baked-contact test on dyke 2 (Fig. 5): host pre-Grenvillian gabbro at the contact
  (sample 108) carries the dyke direction, a hybrid-zone sample at 8.5 m (109) carries both,
  and a sample at 13 m (107) the Grenvillian direction. The authors conclude that the dyke
  remanence is primary -> ST-C, contact_test 'C+'. The host data are shown only in Fig. 5,
  so there are no host-rock rows.
- Ages: site and location rows carry the compilation age, 615 Ma (613-617 Ma), from the
  615 +/- 2 Ma U-Pb date of the swarm (Kamo et al., 1989, on dyke 1; an assigned age that
  satisfies the zircon upper intercept, 615.3 +4.3/-2.1 Ma, as a minimum and the mean
  baddeleyite 207Pb/206Pb age, 615.9 +/- 1.2 Ma, as a maximum, with 2-sigma errors); Kamo & Gower (1994)
  date the primary-remanence dyke 4 at 614 +6/-4 Ma. The ages table holds the two U-Pb
  dates. The K-Ar ages of dykes 2 (514 +/- 8 Ma) and 4 (553 +/- 22 Ma), which the authors
  interpret as post-emplacement cooling ages, are given in the dyke descriptions only.
"""
import os, sys
from datetime import date
import numpy as np
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

CIT = '10.1139/e92-098'
CIT_KG94 = '10.4138/2133'
# Kamo et al. (1989), Geology 17, 602-605: DOI 10.1130/0091-7613(1989)017<0602:BFTLOA>2.3.CO;2
# (verified with Crossref). It contains a colon, which would break MagIC's
# colon-delimited citations list, so it is named in descriptions only.
LOCATION = 'Long Range Dykes'
AGE, AGE_LOW, AGE_HIGH = 615, 613, 617
PRIMARY = {'2', '4', '6'}

# Table 1 (p. 1229): dyke, site, samples, N (collected), trend, width (m),
# D, I, k, a95, n1, n2
T1 = [
    ('2', 'PR84-056', '56-64', 9, '000', '50', 111.2, 50.9, 31, 8.0, 7, 12),
    ('2', 'PR84-098', '98-109', 8, '028', '30', 106.2, 49.9, 45, 7.9, 5, 10),
    ('2', 'PR85-631', '631-645', 5, '020', '50', 131.7, 39.7, 66, 8.5, 3, 8),
    ('4', 'PR84-154', '154-161', 6, '020', '>50', 114.2, 52.4, 60, 9.5, 3, 6),
    ('4', 'PR84-171', '171-177', 5, '000', '35', 158.0, 53.2, 32, 22.3, 3, 0),
    ('4', 'PR84-191', '191-198', 6, '045', '50', 129.0, 48.6, 44, 18.8, 3, 0),
    ('4', 'PR84-199', '199-205', 5, '017', '0.25', 120.5, 51.7, 16, 15.6, 4, 8),
    ('4', 'PR84-206', '206-213', 5, '020', '>20', 136.5, 64.3, 32, 10.4, 3, 10),
    ('4', 'PR85-565', '565-569', 5, '000', '>1', 137.1, 40.0, 104, 6.6, 4, 4),
    ('4', 'PR85-570', '570-576', 5, '025', '>50', 136.8, 43.3, 11, 20.1, 4, 6),
    ('6', 'PR86-008', '008-014', 5, '026', '3', 117.0, 66.0, 79, 5.8, 6, 6),
    ('6', 'PR86-015', '015-021', 5, '020', '1', 132.6, 68.6, 257, 3.1, 7, 5),
    ('1', 'PR84-025', '025-030', 5, '025', '>50', 201.4, -72.1, 103, 7.2, 0, 3),
    ('1', 'PR85-528', '528-534', 4, '025', '>50', 220.5, -85.6, 6, 28.0, 0, 4),
    ('1', 'PR84-051', '051-055', 8, '025', '0.3', 225.5, -67.9, 73, 12.0, 0, 7),
    ('3', 'PR84-162', '162-170', 6, '035', '>10', 249.8, -49.2, 43, 16.4, 2, 4),
    ('5', 'PR84-178', '178-184', 5, '023', '60', 279.1, 28.1, 31, 7.3, 8, 15),
]
# Table 1 dyke means: D, I, k, a95 (None where not printed)
DYKE_MEANS = {'2': (117.3, 47.4, 52, 17.3), '4': (133.1, 51.2, 48, 8.8),
              '6': (124.6, 67.5, None, None), '1': (215.3, -75.5, 67, 15.2),
              '3': (249.8, -49.2, 43, 16.4), '5': (279.1, 28.1, 31, 7.3)}
MEAN3 = (124.8, 55.5, 48, 18.0)            # 'Mean characteristic remanence (N = 3 dykes)'

# site coordinates (deg N, deg E) -- see docstring for how each was obtained
KG94 = (round(53 + 46.2 / 60, 3), round(360 - (56 + 37.1 / 60), 3))
COORD = {
    'PR84-056': (53.652, 302.851, 'fig2'), 'PR84-098': (53.900, 303.020, 'fig2'),
    'PR85-631': (53.447, 302.709, 'fig2'), 'PR84-154': (53.760, 303.353, 'fig2'),
    'PR84-171': (53.780, 303.358, 'fig2'),
    'PR84-191': KG94 + ('kg94',), 'PR84-199': KG94 + ('kg94',), 'PR84-206': KG94 + ('kg94',),
    'PR85-565': (53.680, 303.310, 'fig2'), 'PR85-570': (53.680, 303.310, 'fig2'),
    'PR86-008': (53.30, 304.11, 'fig1'), 'PR86-015': (53.30, 304.11, 'fig1'),
    'PR84-025': (53.610, 302.534, 'fig2'), 'PR85-528': (53.610, 302.534, 'fig2'),
    'PR84-051': (53.610, 302.534, 'fig2'), 'PR84-162': (53.728, 303.309, 'fig2'),
    'PR84-178': (53.769, 303.445, 'fig2'),
}
COORD_NOTE = {
    'fig2': 'Coordinates read from the Fig. 2 sketch map (approximate, ca. 1-2 km).',
    'fig1': 'Coordinates read from the Fig. 1 map (approximate, ca. 5 km).',
    'kg94': ('Coordinates of the U-Pb sample locality of Kamo & Gower (1994), which '
             'they state corresponds to sites PR84-191, PR84-199 and PR84-206.'),
}
DYKE_NOTE = {
    '1': ('Dyke 1, the longest dyke of the swarm in Labrador; coarse main phase '
          '(PR84-025, PR85-528) intruded by a fine-grained porphyritic phase of the '
          'same trend (PR84-051); cut by a Carboniferous dyke, near which the coarse '
          'phase carries a Carboniferous overprint (specimens distant from it give the '
          'tabulated direction). Dated by U-Pb at 615 +/- 2 Ma (Kamo et al., 1989).'),
    '2': 'Dyke 2, the dyke of the baked-contact test; K-Ar whole-rock age 514 +/- 8 Ma.',
    '3': 'Dyke 3.',
    '4': ('Dyke 4, the most thoroughly sampled dyke; U-Pb baddeleyite-zircon age '
          '614 +6/-4 Ma (Kamo & Gower, 1994) and K-Ar biotite age 553 +/- 22 Ma.'),
    '5': 'Dyke 5; samples showed great stability.',
    '6': ('Dyke 6, the most mineralogically modified dyke (foliated, olivine altered, '
          'quartz-epidote veins), with minor compositional differences from dykes 2 and '
          '4 attributed by the authors to low-grade metamorphism.'),
}
PRIM = ('Remanence interpreted by the authors as primary (common SE-and-down '
        'direction of dykes 2, 4 and 6).')
ANOM = ('Anomalous remanence, interpreted by the authors as a disturbed signature; '
        'excluded from the pole.')


def vgp(dec, inc, a95, lat, lon):
    plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95 if a95 else 0.0, lat, lon)
    r = lambda x: round(float(x), 1)
    return r(plat), r(plon), (r(dp) if a95 else None), (r(dm) if a95 else None)


GEO = dict(geologic_classes='Intrusive', geologic_types='Volcanic Dike',
           lithologies='Diabase')
rows = []
for (dyke, site, samples, n, trend, width, dec, inc, k, a95, n1, n2) in T1:
    lat, lon, src = COORD[site]
    plat, plon, dp, dm = vgp(dec, inc, a95, lat, lon)
    codes = ['SO-SM', 'LP-DIR-AF', 'LP-DIR-T'] + (['DE-BFL'] if n1 else []) + \
            (['DE-BFP'] if n2 else []) + ['DE-FM']
    good = dyke in PRIMARY
    desc = (f'Sampling site in dyke {dyke} (samples {samples}); dyke trend {trend} deg, '
            f'width {width} m. ' + (PRIM if good else ANOM) + ' ' + COORD_NOTE[src])
    rows.append({'site': site, 'location': LOCATION, 'result_type': 'i',
                 'result_quality': 'g' if good else 'b', 'method_codes': ':'.join(codes),
                 'citations': CIT, **GEO, 'lat': lat, 'lon': lon,
                 'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
                 'dir_tilt_correction': 0, 'dir_dec': dec, 'dir_inc': inc,
                 'dir_alpha95': a95, 'dir_k': k,
                 'dir_n_specimens': n1 + n2, 'dir_n_specimens_lines': n1,
                 'dir_n_specimens_planes': n2, 'dir_n_total_samples': n,
                 'vgp_lat': plat, 'vgp_lon': plon, 'vgp_dp': dp, 'vgp_dm': dm,
                 'description': desc})
sites_i = pd.DataFrame(rows)

# dyke-mean rows (the averaging unit of the pole)
ORDER = ['SO-SM', 'LP-DIR-AF', 'LP-DIR-T', 'DE-BFL', 'DE-BFP', 'DE-FM']
drows = []
for dyke in ['1', '2', '3', '4', '5', '6']:
    sub = sites_i[sites_i['site'].isin([t[1] for t in T1 if t[0] == dyke])]
    dec, inc, k, a95 = DYKE_MEANS[dyke]
    lat, lon = round(sub['lat'].mean(), 3), round(sub['lon'].mean(), 3)
    plat, plon, dp, dm = vgp(dec, inc, a95, lat, lon)
    # check the printed dyke mean against the Fisher mean of its site means
    f = pmag.fisher_mean(sub[['dir_dec', 'dir_inc']].values.tolist())
    if len(sub) > 1:
        assert abs(f['dec'] - dec) < 0.3 and abs(f['inc'] - inc) < 0.1, (dyke, f)
    good = dyke in PRIMARY
    nsite = len(sub)
    desc = (f'{DYKE_NOTE[dyke]} Dyke mean of {nsite} sampling site(s) '
            f'({", ".join(sub["site"])}), as printed in Table 1'
            + ('; no k or alpha95 printed for this two-site mean' if a95 is None else '')
            + '. VGP computed at the mean coordinates of its sites. '
            + (PRIM if good else ANOM))
    drows.append({'site': f'Dyke {dyke}', 'location': LOCATION, 'result_type': 'a',
                  'result_quality': 'g' if good else 'b',
                  'method_codes': ':'.join(c for c in ORDER if any(
                      c in m.split(':') for m in sub['method_codes'])), 'citations': CIT,
                  **GEO, 'lat': lat, 'lon': lon,
                  'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
                  'dir_tilt_correction': 0, 'dir_dec': dec, 'dir_inc': inc,
                  'dir_alpha95': a95, 'dir_k': k,
                  'dir_n_specimens': int(sub['dir_n_specimens'].sum()),
                  'dir_n_specimens_lines': int(sub['dir_n_specimens_lines'].sum()),
                  'dir_n_specimens_planes': int(sub['dir_n_specimens_planes'].sum()),
                  'dir_n_total_samples': int(sub['dir_n_total_samples'].sum()),
                  'vgp_lat': plat, 'vgp_lon': plon, 'vgp_dp': dp, 'vgp_dm': dm,
                  'description': desc})
sites = pd.concat([pd.DataFrame(drows), sites_i], ignore_index=True)
for c in ['dir_n_specimens', 'dir_n_specimens_lines', 'dir_n_specimens_planes',
          'dir_n_total_samples']:
    sites[c] = sites[c].astype('Int64')

pole_rows = sites[(sites['result_type'] == 'a') & (sites['result_quality'] == 'g')]
assert list(pole_rows['site']) == ['Dyke 2', 'Dyke 4', 'Dyke 6']
blk = pmag.flip(ipmag.make_di_block(pole_rows['vgp_lon'].tolist(), pole_rows['vgp_lat'].tolist()),
                combine=True)
p = ipmag.fisher_mean(di_block=blk)
if p['inc'] > 0:                     # report at the southern pole, as in the paper
    p = ipmag.fisher_mean(di_block=ipmag.do_flip(di_block=blk))
d = ipmag.fisher_mean(di_block=pole_rows[['dir_dec', 'dir_inc']].values.tolist())
assert abs(d['dec'] - MEAN3[0]) < 0.1 and abs(d['inc'] - MEAN3[1]) < 0.1

locs = pd.DataFrame([{
    'location': LOCATION, 'location_type': 'Region',
    'result_name': 'Long Range dykes ca. 615 Ma pole (dykes 2, 4, 6)',
    'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(pole_rows['site']),
    'method_codes': 'SO-SM:LP-DIR-AF:LP-DIR-T:DE-BFL:DE-BFP:DE-FM:DE-VGP:ST-C',
    'citations': f'{CIT}:{CIT_KG94}',
    'geologic_classes': 'Intrusive', 'lithologies': 'Diabase',
    'lat_s': sites['lat'].min(), 'lat_n': sites['lat'].max(),
    'lon_w': sites['lon'].min(), 'lon_e': sites['lon'].max(),
    'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
    'dir_tilt_correction': 0, 'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
    'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1),
    'dir_n_sites': int(d['n']),
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
    'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
    'pole_n_sites': int(p['n']), 'pole_reversed_perc': 0, 'contact_test': 'C+',
    'external_database_ids': 'GPMDB[6934]',
    'continent_ocean': 'North America', 'country': 'Canada',
    'state_province': 'Newfoundland and Labrador', 'terranes': 'Laurentia',
    'description': (
        "Paleomagnetic pole for the Long Range dykes of the Sandwich Bay area, southeast "
        "Labrador (Grenville Province), a swarm of near-vertical, NNE-trending, coarse, "
        "ophitic mafic dykes correlated with the Long Range dykes of northwest "
        "Newfoundland (Murthy et al., 1992). The pole is the Fisher mean of the virtual "
        "geomagnetic poles of the three dykes (2, 4 and 6) that share a southeast-and-"
        "down characteristic direction interpreted by the authors as primary; each dyke "
        "is treated as one cooling unit and its VGP is computed from the dyke-mean "
        "direction. It is reported at the southern-hemisphere pole as in the paper. "
        "Three other dykes (1, 3 and 5) carry mutually different 'anomalous' directions "
        "that the authors interpret as disturbed paleomagnetic signatures and exclude. "
        "The remanence is carried by magnetite (Curie point "
        "570 C; multidomain behaviour in the Lowrie-Fuller test) and was isolated mostly "
        "by alternating-field demagnetization, with thermal demagnetization of some "
        "samples, principal-component analysis and remagnetization circles. A "
        "baked-contact test on dyke 2 (one host-gabbro sample each at the contact, in a "
        "hybrid zone at 8.5 m and at 13 m) is interpreted by the authors as showing the "
        "dyke remanence to be primary. Site coordinates are not given in the paper and "
        "were read from its maps. U-Pb baddeleyite-zircon dating of dyke 4 (Kamo & Gower, 1994) ties "
        "the primary direction to the ca. 615 Ma emplacement age of the swarm "
        "determined on dyke 1 (Kamo et al., 1989); K-Ar ages on dykes 2 and 4 are "
        "interpreted as post-emplacement cooling ages.")}])

# ---- ages table: measured dates -------------------------------------------------
AGE_ROWS = [
    dict(site='Dyke 4', age=614, age_low=610, age_high=620, method_codes='GM-UPB-CC-T0',
         citations=CIT_KG94, period='Ediacaran', era='Neoproterozoic', eon='Proterozoic',
         description=('U-Pb upper-intercept age of three baddeleyite and two zircon '
                      'fractions from a coarse-grained olivine gabbro of dyke 4 at the '
                      'locality of sites PR84-191, PR84-199 and PR84-206, 614 +6/-4 Ma '
                      '(95% confidence) assuming Pb loss between 0 and 300 Ma; the '
                      'unconstrained regression gives 614 +22/-4 Ma (Kamo & Gower, 1994; '
                      'given as 615 +6/-4 Ma in their abstract and discussion).')),
    dict(site='Dyke 1', age=615, age_low=613, age_high=617, method_codes='GM-UPB',
         citations=f'{CIT}:{CIT_KG94}', period='Ediacaran', era='Neoproterozoic',
         eon='Proterozoic',
         description=('U-Pb zircon and baddeleyite age of dyke 1, the longest dyke of the '
                      'swarm in Labrador, 615 +/- 2 Ma (Kamo et al., 1989, Geology 17, '
                      '602-605). The assigned age satisfies both the zircon upper-intercept '
                      'age of 615.3 +4.3/-2.1 Ma, taken as a minimum, and the mean '
                      '207Pb/206Pb baddeleyite age of 615.9 +/- 1.2 Ma, taken as a maximum '
                      '(analytical errors 2 sigma); age_low and age_high are the assigned '
                      'range.')),
]
ages = pd.DataFrame([{
    'location': LOCATION, 'site': r['site'], 'age': r['age'], 'age_low': r['age_low'],
    'age_high': r['age_high'], 'age_unit': 'Ma', 'method_codes': r['method_codes'],
    'citations': r['citations'], 'timescale_eon': r['eon'], 'timescale_era': r['era'],
    'timescale_period': r['period'], 'description': r['description']} for r in AGE_ROWS])


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
combined = os.path.join(OUT, f'Murthy1992_Long_Range_dykes_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {combined}')

r = locs.iloc[0]
print(f"  pole {r.pole_lat}/{r.pole_lon} A95 {r.pole_alpha95} K {r.pole_k} N {r.pole_n_sites}; "
      f"dir {r.dir_dec}/{r.dir_inc} a95 {r.dir_alpha95} k {r.dir_k}")
print(pole_rows[['site', 'lat', 'lon', 'vgp_lat', 'vgp_lon']].to_string(index=False))

# pole of the Swanson-Hysell et al. (2026) compilation (lat, lon, A95, B, N); N is the
# number of samples collected in the three dykes (samples used are not reported)
EXPECT = (-12.3, 165.4, 20.0, 3, 69)
got = (r.pole_lat, r.pole_lon, r.pole_alpha95, int(r.pole_n_sites),
       int(pole_rows['dir_n_total_samples'].sum()))
assert got == EXPECT, got
print('-I- the pole matches the compilation')

from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(combined, tables=['locations', 'sites', 'ages'], server=True)
