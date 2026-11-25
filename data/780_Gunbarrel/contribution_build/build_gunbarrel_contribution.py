"""Build the MagIC contribution for the ca. 780 Ma Gunbarrel pole of Ding et al. (2025).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Ding, J., Evans, D. A. D., Kilian, T., Mitchell, R. N., Gong, Z., Chamberlain, K., Bleeker, W.,
& Zhang, S. (2025). Local rotations biasing supercontinental configurations: Revisiting a key
ca. 780 Ma paleomagnetic pole for Laurentia and Rodinia. Journal of Geophysical Research: Solid
Earth 130, e2025JB031762. doi:10.1029/2025JB031762

Earlier site results compiled by Ding et al. (2025): Harlan, S. S., Geissman, J. W., Snee, L. W.,
& Reynolds, R. L. (1997), USGS Professional Paper 1580, doi:10.3133/pp1580 (BT-8, MM1); Park,
J. K., Buchan, K. L., & Gandhi, S. S. (1995), Geological Survey of Canada Current Research 1995-C,
195-200 (Gunbarrel, Margaret and Calder sheets; no DOI).
Dates: Harlan, S. S., Heaman, L., LeCheminant, A. N., & Premo, W. R. (2003). Gunbarrel mafic
magmatic event: A key 780 Ma time marker for Rodinia plate reconstructions. Geology 31,
1053-1056. doi:10.1130/G19944.1. Uncertainties are 95% confidence (2 sigma); age_sigma is 1 sigma.

Sources and checks
------------------
- Ding et al. (2025) Table 1 (p. 9) is transcribed below from the born-digital text layer and
  checked against the page image: every row, including the AF sister-sample rows (italic, not
  used), the three Calder results averaged into the Calder mean, and the baked and unbaked
  host-rock sites of the two contact tests. The 14 bold rows are the units of the pole
  (L780: 3.2 N, 151.5 E, A95 8.0; reproduced from the tabulated VGPs below).
- Directions, VGPs and dp/dm are stored as printed. For the Mount Moran dike (MM1), Table 1
  pairs the tilt-corrected direction (285.5/1.5) with the VGP and dp/dm of the in-situ
  direction (15.2 N, 152.6 E, 5.3/10.4; Harlan et al., 1997, p. 7, give 11.6 N, 148.9 E for
  the tilt-corrected direction, whose tilt is inferred). The tabulated VGP is stored, so the
  pole reproduces L780; with the tilt-corrected VGP the pole would be 3.0 N, 151.2 E, A95 7.9.
  Margaret and the Calder mean print one value in the dp/dm column (14.8, 12.7) and no k;
  dp/dm are left blank for them.
- Bedding: the Beartooth Mountains sites are tilt-corrected with the Beartooth Butte strata
  (dip direction 198, dip 10; Harlan et al., 1997, as used by Ding et al., 2025), which
  reproduces every tabulated tilt-corrected direction there; that bedding is stored for them.
  MM1 is corrected with the inferred Teton Range tilt of Harlan et al. (1997; strike N10E,
  10 deg west, which they call problematic), stored as dip direction 280, dip 10. The
  Fishtrap Gabbro host (R10HL11A/B) is corrected with bedding that is not tabulated, so
  bed_dip is left blank; the Slave craton sheets and the Christmas Lake host dike have
  identical in-situ and tilt-corrected directions.
- Coordinates are those of Table 1 (degrees W converted to E). The Calder mean is printed at
  the R10HL9 coordinates.

Structure and interpretation
----------------------------
- One location, 'Gunbarrel LIP', holds every site. Its first row has no pole and lists every
  site; the second row is the pole of Ding et al. (2025), which is the pole used in the
  Swanson-Hysell et al. (2026) compilation: the Fisher mean of the 14 tabulated VGPs of 11
  Wyoming craton dikes (9 of Ding et al., 2025; BT-8 and MM1 of Harlan et al., 1997) and the
  Gunbarrel, Margaret and Calder sheets of the Slave craton.
- The Calder mean (result_type 'a') averages the Calder result of Park et al. (1995) and sites
  R10HL8 and R10HL9 of Ding et al. (2025) ('i' rows). The AF rows of four dikes, the Calder
  constituents and the host-rock sites of the contact tests are in no pole ('g' rows).
- contact_test C+ on the pole row: positive baked-contact tests at the Christmas Lake dike
  (J20B24; baked J20B25, unbaked J20B26 in an older NW-SE dike) and at the Calder sheet (R10HL8;
  baked R10HL11A, unbaked R10HL11B in the Fishtrap Gabbro). No reversal test (one polarity).
- Ages: every Gunbarrel site and the pole carry the Gunbarrel event age of the compilation,
  779.5 Ma (776.5-782.5 Ma). The host rocks of the contact tests carry rock-age brackets, not
  dates: the undated NW-SE dike 776.5-2600 Ma (older than the Christmas Lake dike; dike ages
  in the Archean basement uplifts of the Wyoming province range from ca. 2.6 to 0.75 Ga,
  Harlan et al., 1997) and the Fishtrap Gabbro 1875-1900 Ma (Hottah succession, ca. 1895 Ma,
  older than the 1875-1840 Ma Great Bear arc; Ding et al., 2025). The ages
  table holds the U-Pb baddeleyite dates of Harlan et al. (2003) on BT-8 (Christmas Lake dike),
  the Gunbarrel sheet and the Calder sheet. The 40Ar/39Ar dates of Harlan et al. (1997) are
  given in the site descriptions.
"""
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(OUT))

DING = '10.1029/2025JB031762'
HARLAN97 = '10.3133/pp1580'
HARLAN03 = '10.1130/G19944.1'
LOC = 'Gunbarrel LIP'
AGE = (779.5, 776.5, 782.5)
BEARTOOTH_BED = (198, 10)          # dip direction, dip (Beartooth Butte; Harlan et al., 1997)
MM1_BED = (280, 10)                # inferred Teton Range tilt, strike N10E, 10 W (Harlan et al., 1997)
HOST_AGE = {'J20B25': (None, 776.5, 2600), 'J20B26': (None, 776.5, 2600),
            'R10HL11A': (None, 1875, 1900), 'R10HL11B': (None, 1875, 1900)}

# Ding et al. (2025) Table 1: site, role, lat N, lon W, n, N, Dg, Ig, Ds, Is, k, a95,
# Plat, Plon, dp, dm
# role: 'pole' (bold, in L780), 'af' (italic AF sister samples), 'cal' (Calder constituent),
# 'calmean' (Calder mean, in L780), 'host' (contact-test host rock)
T1 = [
    ('BT50', 'pole', 45.028, 109.884, 8, 8, 265.2, -11.4, 267.4, -15.1, 38.1, 9.9, -7.3, 156.5, 5.2, 10.2),
    ('BT51', 'pole', 45.027, 109.883, 8, 8, 268.7, 5.1, 268.1, 1.7, 47, 8.2, -0.7, 162.1, 4.1, 8.2),
    ('BT51AF', 'af', 45.027, 109.883, 3, 4, 278.3, 9.5, 276.8, 7.7, 65.6, 15.3, 7.5, 158.0, 7.8, 15.4),
    ('BT54', 'pole', 45.102, 109.911, 7, 8, 267.2, -18.5, 270.7, -21.8, 19.5, 14, -7.5, 151.6, 7.8, 14.8),
    ('BT5658', 'pole', 45.102, 109.909, 6, 6, 275.5, 7.8, 274.4, 5.5, 44.1, 10.2, 5.1, 158.9, 5.1, 10.2),
    ('BT5658AF', 'af', 45.102, 109.909, 2, 2, 281.3, 14.7, 278.8, 13.3, 967.2, 8, 11.0, 158.6, 4.2, 8.2),
    ('BT64', 'pole', 45.095, 109.909, 4, 4, 264.2, 14.4, 262.2, 10.2, 61.1, 11.9, -1.8, 169.2, 6.1, 12.0),
    ('BT68', 'pole', 45.034, 109.662, 8, 8, 256.5, 5.5, 256.1, 0.2, 35.2, 9.5, -9.7, 170.3, 4.8, 9.5),
    ('BT69', 'pole', 45.004, 109.631, 6, 8, 289, -35.3, 296, -34.5, 137.3, 5.7, 3.6, 128.8, 3.8, 6.5),
    ('BT69AF', 'af', 45.004, 109.631, 3, 5, 281.1, -34.8, 288.1, -35.4, 22.1, 26.9, -1.7, 134.0, 18.0, 31.1),
    ('BT72', 'pole', 45.034, 109.661, 8, 8, 263.7, 5.9, 263.1, 1.7, 30.1, 10.3, -4.3, 165.8, 5.1, 10.3),
    ('J20B24', 'pole', 44.960, 109.388, 7, 8, 289.6, -29.6, 295.2, -28.8, 169.2, 4.7, 5.9, 131.9, 2.8, 5.1),
    ('J20B24AF', 'af', 44.960, 109.388, 5, 6, 304.4, -24.7, 308.5, -21.5, 22.4, 16.5, 17.2, 124.1, 9.2, 17.4),
    ('BT-8', 'pole', 44.971, 109.445, 9, None, 288.3, -21.8, 292.3, -21.4, 44, 7.9, 7.4, 137.0, 4.4, 8.3),
    ('MM1', 'pole', 43.818, 110.742, 7, None, 285.6, 11.5, 285.5, 1.5, 36, 10.2, 15.2, 152.6, 5.3, 10.4),
    ('Gun', 'pole', 65.667, 118.143, 7, 7, 281, 3, 281, 3, 95, 6, 5.9, 142.4, 3.0, 6.0),
    ('Mar', 'pole', 64.401, 117.132, 4, 4, 280.4, 16.1, 280.4, 16.1, None, 20, 11.9, 147.0, None, None),
    ('Cal', 'cal', 65.527, 117.305, 6, 6, 292, 27, 292, 27, 36, 10, 22.0, 138.4, 5.9, 10.9),
    ('R10HL8', 'cal', 65.011, 118.430, 8, 8, 284, 36.7, 284, 36.7, 80.4, 6.2, 24.3, 147.9, 4.2, 7.2),
    ('R10HL9', 'cal', 65.011, 118.463, 4, 8, 275.5, 27.1, 275.5, 27.1, 226.7, 6.1, 15.3, 152.8, 3.6, 6.6),
    ('Cal_mean', 'calmean', 65.011, 118.463, 3, 3, 283.5, 30.6, 283.5, 30.6, None, 15.3, 20.6, 146.4, None, None),
    ('J20B25', 'host', 44.960, 109.388, 8, 8, 295.5, -42.6, 295.5, -42.6, 161.4, 4.4, None, None, None, None),
    ('J20B26', 'host', 44.961, 109.388, 5, 8, 203.8, 43.4, 203.8, 43.4, 15.5, 20, None, None, None, None),
    ('R10HL11A', 'host', 65.017, 118.397, 4, 6, 274.3, 14.9, 276.3, 5, 20.9, 20.6, None, None, None, None),
    ('R10HL11B', 'host', 65.017, 118.397, 4, 10, 358.9, -30.8, 5.8, -42.7, 24.4, 19, None, None, None, None),
]
BEARTOOTH = {'BT50', 'BT51', 'BT51AF', 'BT54', 'BT5658', 'BT5658AF', 'BT64', 'BT68', 'BT69',
             'BT69AF', 'BT72', 'J20B24', 'J20B24AF', 'BT-8'}
DING_SITES = {'BT50', 'BT51', 'BT51AF', 'BT54', 'BT5658', 'BT5658AF', 'BT64', 'BT68', 'BT69',
              'BT69AF', 'BT72', 'J20B24', 'J20B24AF', 'R10HL8', 'R10HL9', 'Cal_mean', 'J20B25',
              'J20B26', 'R10HL11A', 'R10HL11B'}
PARK95 = {'Gun', 'Mar', 'Cal'}
HARLAN = {'BT-8', 'MM1'}

MC_DING = 'SO-SUN:SO-MAG:LP-DIR-T:DE-BFL:DE-FM'
MC_DING_AF = 'SO-SUN:SO-MAG:LP-DIR-AF:DE-BFL:DE-FM'
MC_LEGACY = 'DE-FM'

DESC = {
    'BT50': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'BT51': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'BT54': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'BT5658': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'BT64': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'BT68': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'BT69': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'BT72': 'W-E to WNW-ESE diabase dike, Beartooth Mountains, Wyoming craton (Ding et al., 2025).',
    'J20B24': ('Christmas Lake dike, ca. 10 m wide, Beartooth Mountains, Wyoming craton (Ding et al., '
               '2025); the dike of the baked-contact test with J20B25 (baked) and J20B26 (unbaked).'),
    'BT-8': ('Christmas Lake dike, Beartooth Mountains, Wyoming craton; paleomagnetic result of Harlan '
             'et al. (1997) as compiled by Ding et al. (2025). U-Pb baddeleyite date 779.5 +/- 3.0 Ma '
             '(Harlan et al., 2003); 40Ar/39Ar hornblende plateau date of Harlan et al. (1997).'),
    'MM1': ('Mount Moran dike, Teton Range, Wyoming; paleomagnetic result of Harlan et al. (1997) as '
            'compiled by Ding et al. (2025). 40Ar/39Ar date 775 +/- 10 Ma (Harlan et al., 1997; ca. '
            '772 Ma as recalculated by Harlan et al., 2003). Tilt-corrected for the inferred 10 deg '
            'west tilt of the Teton Range (Harlan et al., 1997), which they call problematic. The '
            'VGP is as tabulated by Ding et al. (2025) and is that of the in-situ direction; '
            'Harlan et al. (1997) give 11.6 N, 148.9 E for the tilt-corrected direction.'),
    'Gun': ('Gunbarrel gabbro sheet (Hottah sheets), Wopmay orogen, Slave craton; paleomagnetic result '
            'of Park et al. (1995) as compiled by Ding et al. (2025). U-Pb baddeleyite date 779.6 +/- '
            '1.4 Ma (Harlan et al., 2003).'),
    'Mar': ('Margaret gabbro sheet (Hottah sheets), Wopmay orogen, Slave craton; paleomagnetic result '
            'of Park et al. (1995) as compiled by Ding et al. (2025), which print no k and a single '
            'pole uncertainty of 14.8 deg.'),
    'Cal': ('Calder gabbro sheet (Hottah sheets), Hottah Lake, Slave craton; paleomagnetic result of '
            'Park et al. (1995) as compiled by Ding et al. (2025). U-Pb baddeleyite date 779.5 +/- 1.8 '
            'Ma (Harlan et al., 2003). Averaged into the Calder mean.'),
    'R10HL8': ('Calder gabbro-diabase sheet, Hottah Lake, Slave craton (Ding et al., 2025); the sheet of '
               'the baked-contact test with R10HL11A (baked) and R10HL11B (unbaked). Averaged into the '
               'Calder mean.'),
    'R10HL9': ('Calder gabbro-diabase sheet, Hottah Lake, Slave craton (Ding et al., 2025). Averaged into '
               'the Calder mean.'),
    'Cal_mean': ('Mean of the three Calder sheet results (Park et al., 1995; R10HL8 and R10HL9 of Ding '
                 'et al., 2025), as calculated by Ding et al. (2025), who print it at the R10HL9 '
                 'coordinates, with no k and a single pole uncertainty of 12.7 deg.'),
    'J20B25': ('Baked host rock of the Christmas Lake dike (J20B24): older NW-SE mafic dike sampled ca. '
               '1-3 m from the dike margin (Ding et al., 2025). The direction resembles that of the '
               'Christmas Lake dike. This baked host rock of the contact test is not included in a pole.'),
    'J20B26': ('Unbaked host rock of the Christmas Lake dike test: the older NW-SE mafic dike sampled ca. '
               '130 m from the intersection (Ding et al., 2025), with a direction that Ding et al. '
               '(2025) compare with ca. 1500 Ma Laurentian results; the dike is undated. This host rock '
               'of the contact test is not included in a pole.'),
    'R10HL11A': ('Baked host rock of the Calder sheet (R10HL8): Fishtrap Gabbro (ca. 1895 Ma; Davis et '
                 'al., 2015; Ootes et al., 2015, as cited by Ding et al., 2025) within ca. 4.5 m below '
                 'the sheet. Tilt-corrected with bedding that is not tabulated. This baked host rock of '
                 'the contact test is not included in a pole.'),
    'R10HL11B': ('Unbaked host rock of the Calder sheet test: Fishtrap Gabbro (ca. 1895 Ma; Davis et al., '
                 '2015; Ootes et al., 2015, as cited by Ding et al., 2025) ca. 26-52 m below the sheet. '
                 'Tilt-corrected with bedding that is not tabulated. This host rock of the contact test '
                 'is not included in a pole.'),
}
for s in ('BT51AF', 'BT5658AF', 'BT69AF', 'J20B24AF'):
    DESC[s] = (f'Alternating-field demagnetization of sister samples of dike {s[:-2]} (Ding et al., '
               f'2025); not used in the pole, which uses the thermal result.')


def r(x):
    return None if x is None else x


rows = []
for (s, role, lat, lonw, n, n0, dg, ig, ds, i_s, k, a95, plat, plon, dp, dm) in T1:
    lon = round(360 - lonw, 3)
    if s in DING_SITES:
        cit, mc = DING, (MC_DING_AF if role == 'af' else MC_DING)
    elif s in HARLAN:
        cit, mc = f'{HARLAN97}:{DING}', MC_LEGACY
    else:
        cit, mc = DING, MC_LEGACY      # Park et al. (1995) has no DOI
    if role == 'host':
        gclass, gtype, lith = 'Intrusive', ('Volcanic Dike' if s.startswith('J20') else 'Pluton'), \
            ('Diabase' if s.startswith('J20') else 'Gabbro')
        age = HOST_AGE[s]
    else:
        gclass = 'Intrusive'
        if s in PARK95 or s.startswith('R10HL') or s == 'Cal_mean':
            gtype, lith = 'Sill', 'Gabbro:Diabase'
        else:
            gtype, lith = 'Volcanic Dike', 'Diabase'
        age = AGE
    bdd, bd = BEARTOOTH_BED if s in BEARTOOTH else (MM1_BED if s == 'MM1' else (None, None))
    rtype = 'a' if role == 'calmean' else 'i'
    for tc, d, i in ((0, dg, ig), (100, ds, i_s)):
        rows.append({
            'site': s, 'location': LOC, 'result_type': rtype, 'result_quality': 'g',
            'method_codes': mc, 'citations': cit, 'geologic_classes': gclass,
            'geologic_types': gtype, 'lithologies': lith, 'bed_dip_direction': bdd,
            'bed_dip': bd, 'lat': lat, 'lon': lon, 'age': age[0], 'age_low': age[1],
            'age_high': age[2], 'age_unit': 'Ma', 'dir_tilt_correction': tc,
            'dir_dec': d, 'dir_inc': i, 'dir_alpha95': a95, 'dir_k': k, 'dir_n_samples': n,
            'dir_n_total_samples': n0,
            'vgp_lat': plat if tc == 100 else None, 'vgp_lon': plon if tc == 100 else None,
            'vgp_dp': dp if tc == 100 else None, 'vgp_dm': dm if tc == 100 else None,
            'description': DESC[s]})
sites = pd.DataFrame(rows)
for c in ['dir_n_samples', 'dir_n_total_samples', 'bed_dip', 'bed_dip_direction',
          'dir_tilt_correction']:
    sites[c] = pd.to_numeric(sites[c]).astype('Int64')

# the Beartooth bedding reproduces the tabulated tilt-corrected directions of those sites
for s in BEARTOOTH:
    g = sites[(sites.site == s) & (sites.dir_tilt_correction == 0)].iloc[0]
    t = sites[(sites.site == s) & (sites.dir_tilt_correction == 100)].iloc[0]
    d, i = pmag.dotilt(g.dir_dec, g.dir_inc, *BEARTOOTH_BED)
    assert pmag.angle([d, i], [t.dir_dec, t.dir_inc])[0] < 0.6, s

POLE_UNITS = [x[0] for x in T1 if x[1] in ('pole', 'calmean')]
assert len(POLE_UNITS) == 14
tc = sites[(sites.dir_tilt_correction == 100) & sites.site.isin(POLE_UNITS)]
p = pmag.fisher_mean([[lo, la, 1.] for la, lo in zip(tc.vgp_lat, tc.vgp_lon)])
assert (round(p['inc'], 1), round(p['dec'], 1), round(p['alpha95'], 1)) == (3.2, 151.5, 8.0), p
# with the VGP of MM1's tilt-corrected direction instead, the pole would be 3.0/151.2, A95 7.9
alt = [[lo, la, 1.] for la, lo in zip(tc.vgp_lat, tc.vgp_lon)]
alt[list(tc.site).index('MM1')] = [148.9, 11.6, 1.]
p_alt = pmag.fisher_mean(alt)
assert (round(p_alt['inc'], 1), round(p_alt['dec'], 1), round(p_alt['alpha95'], 1)) == (3.0, 151.2, 7.9)


def r1(x):
    return round(float(x), 1)


GEOG = dict(continent_ocean='North America', terranes='Laurentia')
ALL_CODES = ':'.join(dict.fromkeys(':'.join(sites.method_codes).split(':')))
all_sites = list(dict.fromkeys(sites.site))
locs = pd.DataFrame([
    {'location': LOC, 'location_type': 'Region', 'sites': ':'.join(all_sites),
     'method_codes': ALL_CODES, 'citations': f'{DING}:{HARLAN97}:{HARLAN03}',
     'geologic_classes': 'Intrusive', 'lithologies': 'Diabase:Gabbro',
     'lat_s': sites.lat.min(), 'lat_n': sites.lat.max(), 'lon_w': sites.lon.min(),
     'lon_e': sites.lon.max(), 'age': AGE[0], 'age_low': AGE[1], 'age_high': AGE[2],
     'age_unit': 'Ma', **GEOG,
     'description': (
         "Ca. 780 Ma Gunbarrel intrusions of the Wyoming craton (Beartooth Mountains, Montana and "
         "Wyoming; Mount Moran, Teton Range, Wyoming) and the Slave craton (Hottah sheets, "
         "Northwest Territories, Canada), from Table 1 of Ding et al. (2025), including earlier "
         "results of Harlan et al. (1997) and Park et al. (1995) compiled there. Sites in no pole: "
         "the alternating-field results of four dikes (the pole uses their thermal results), the "
         "three Calder sheet results that are averaged into the Calder mean, and the baked and "
         "unbaked host rocks of the two contact tests, whose ages are rock-age brackets rather than "
         "dates (the NW-SE host dike is undated; the Fishtrap Gabbro is part of the ca. 1.9 Ga "
         "Hottah succession).")},
    {'location': LOC, 'location_type': 'Region', 'result_name': 'Gunbarrel LIP ca. 780 Ma pole',
     'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(POLE_UNITS),
     'method_codes': ALL_CODES + ':DE-VGP:ST-C', 'citations': f'{DING}:{HARLAN97}:{HARLAN03}',
     'geologic_classes': 'Intrusive', 'lithologies': 'Diabase:Gabbro',
     'lat_s': tc.lat.min(), 'lat_n': tc.lat.max(), 'lon_w': tc.lon.min(), 'lon_e': tc.lon.max(),
     'age': AGE[0], 'age_low': AGE[1], 'age_high': AGE[2], 'age_unit': 'Ma',
     'dir_tilt_correction': 100,
     'pole_lat': r1(p['inc']), 'pole_lon': r1(p['dec']), 'pole_alpha95': r1(p['alpha95']),
     'pole_k': r1(p['k']), 'pole_n_sites': 14, 'pole_reversed_perc': 0, 'contact_test': 'C+',
     **GEOG,
     'description': (
         "Paleomagnetic pole of Ding et al. (2025, L780) for the ca. 780 Ma Gunbarrel magmatic "
         "event: the Fisher mean of the virtual geomagnetic poles of 14 intrusions, 11 dikes of "
         "the Wyoming craton (nine of Ding et al., 2025; BT-8 and MM1 of Harlan et al., 1997) and "
         "the Gunbarrel, Margaret and Calder sheets of the Slave craton (Park et al., 1995; the "
         "Calder mean adds two sites of Ding et al., 2025). The site VGPs are those tabulated by "
         "Ding et al. (2025), in tilt-corrected coordinates except for the Mount Moran dike, whose "
         "tabulated VGP is that of its in-situ direction (its tilt correction is inferred; with the "
         "VGP of the tilt-corrected direction the pole would be 3.0 N, 151.2 E, A95 7.9). The "
         "high-temperature component is "
         "carried by low-titanium magnetite and is of single, westerly polarity. Primary origin "
         "is supported by positive baked-contact tests at the Christmas Lake dike (Beartooth "
         "Mountains) and the Calder sheet (Hottah Lake). Ding et al. (2025) exclude results from "
         "the Tobacco Root Mountains and the Mackenzie Mountains, which they interpret as "
         "affected by local vertical-axis rotations. The age is that of the Gunbarrel event from "
         "U-Pb baddeleyite dates (Harlan et al., 2003).")},
])
for c in ['pole_n_sites', 'dir_tilt_correction']:
    locs[c] = pd.to_numeric(locs[c]).astype('Int64')
LOC_ORDER = ['location', 'location_type', 'result_name', 'result_type', 'result_quality', 'sites',
             'method_codes', 'citations', 'geologic_classes', 'lithologies', 'lat_s', 'lat_n',
             'lon_w', 'lon_e', 'age', 'age_low', 'age_high', 'age_unit', 'dir_tilt_correction',
             'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k', 'pole_n_sites',
             'pole_reversed_perc', 'contact_test', 'continent_ocean', 'terranes', 'description']
assert set(LOC_ORDER) == set(locs.columns)
locs = locs[LOC_ORDER]

AGE_ROWS = [  # site, date, 2 sigma, text
    ('BT-8', 779.5, 3.0, 'Christmas Lake dike, Beartooth Mountains'),
    ('Gun', 779.6, 1.4, 'Gunbarrel gabbro sheet, Slave craton (mean of two analyses)'),
    ('Cal', 779.5, 1.8, 'Calder gabbro sheet, Slave craton (mean of three fractions)'),
]
ages = pd.DataFrame([{
    'location': LOC, 'site': s, 'age': a, 'age_sigma': e / 2, 'age_unit': 'Ma',
    'method_codes': 'GM-UPB', 'citations': HARLAN03, 'timescale_eon': 'Proterozoic',
    'timescale_era': 'Neoproterozoic', 'timescale_period': 'Tonian',
    'description': (f'U-Pb ID-TIMS baddeleyite weighted-mean 207Pb/206Pb date, {txt}; {a} +/- {e} '
                    f'Ma (95% confidence; Harlan et al., 2003).')} for s, a, e, txt in AGE_ROWS])


def table_text(df, kind):
    df = df.dropna(axis=1, how='all')
    return f'tab delimited\t{kind}\n' + df.to_csv(sep='\t', index=False)


TABLES = (('locations', locs), ('sites', sites), ('ages', ages))
for kind, df in TABLES:
    with open(os.path.join(OUT, f'{kind}.txt'), 'w', encoding='utf-8') as f:
        f.write(table_text(df, kind))
stamp = date.today().strftime('%d.%b.%Y')
upload = os.path.join(OUT, f'Ding2025_Gunbarrel_LIP_{stamp}.txt')
with open(upload, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {upload}: {len(locs)} locations, {len(sites)} site rows '
      f'({sites.site.nunique()} sites), {len(ages)} ages')
print(f"  pole {r1(p['inc'])}/{r1(p['dec'])} A95 {r1(p['alpha95'])} K {r1(p['k'])} N 14")

sys.path.insert(0, os.path.join(REPO, 'scripts'))
from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True)
