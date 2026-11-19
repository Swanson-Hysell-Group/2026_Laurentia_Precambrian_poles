"""Build the MagIC contribution for the Abitibi dyke swarm (Ernst and Buchan, 1993).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Ernst, R. E., & Buchan, K. L. (1993). Paleomagnetism of the Abitibi dyke swarm, southern
Superior Province, and implications for the Logan Loop. Can. J. Earth Sci. 30,
1886-1897. doi:10.1139/e93-167. Table 1 (ages; uncertainties 2 sigma), Table 2 (site
results), Table 3 (dyke means and poles). The K-Ar biotite dates compiled in Table 1
are not entered; the dykes have U-Pb dates.

Also used:
  - Halls, H. C., Stott, G. M., & Davis, D. W. (2005). Paleomagnetism, geochronology and
    geochemistry of several Proterozoic mafic dike swarms in northwestern Ontario.
    Ontario Geological Survey Open File Report 6171 (no DOI). Table 4 (site SL13), the
    U-Pb date of the Hillsport dike (site DP4 = dyke A1 of Ernst and Buchan), which
    places dyke A1 in the ca. 2171 Ma Biscotasing swarm.
  - Krogh et al. (1987), Geological Association of Canada Special Paper 34, 147-152
    (no DOI): 1140.6 +/- 2 Ma (95%) U-Pb baddeleyite, Great Abitibi dyke (site 17).
  - Piispa et al. (2018), doi:10.1029/2018JB015992, Table 1: 1142 +/- 3 Ma U-Pb
    baddeleyite of dyke A3 (L. M. Heaman, unpublished, cited in Sage, 2000).

Sources and checks
------------------
- Table 2 was transcribed from the page image; it agrees with an independent
  transcription for every coordinate, direction, alpha95 and k. alpha95
  recalculated from k and n agrees with every site (asserted below), and each dyke mean
  of Table 3 is reproduced from the site D'/I' (asserted below).
- VGPs are calculated from the directions and stored to 0.01 deg in the polarity of the
  direction. The dyke-mean VGPs use the common reference locality 48.0 N, 81.0 W of the
  dyke means and reproduce the poles of Table 3 within 0.07 deg.

Structure and interpretation
----------------------------
- Site rows are the cooling units: the dyke means of Table 3 (at the reference locality),
  SL13, and the baked and unbaked host rocks of the baked-contact test. The sampling-site
  results of Table 2 are listed in the description of each dyke mean.
- Location 'Abitibi dyke swarm': an umbrella row and the eight-unit Abitibi pole used in
  the Swanson-Hysell et al. (2026) compilation: the Fisher mean of the VGPs of dykes A3,
  A4, A5, A7, A8, units 1 and 2 of the Great Abitibi dyke (A6), and site SL13 of Halls
  et al. (2005). The published normal-polarity pole and all-dyke pole of Ernst and Buchan
  (1993) are described, not entered: the first is a single-polarity subset and the
  second includes dyke A1.
- Location 'Hillsport dyke (Biscotasing swarm)': dyke A1, dated at 2170.7 +/- 1.1 Ma by
  Halls et al. (2005).
- Site ages are ages of the rocks: dated dykes carry their U-Pb dates; the host rocks of
  the baked-contact test carry the age of the Matachewan swarm. The Keweenawan overprint
  poles that Ernst and Buchan report from dyke A3 are described in its rows, not entered
  as poles.
"""
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

PIISPA = '10.1029/2018JB015992'
HALLS = 'Halls et al. (2005) Ontario Geological Survey Open File Report 6171'
KROGH = 'Krogh et al. (1987) Geological Association of Canada Special Paper 34, 147-152'
LOC = 'Abitibi dyke swarm'
LOC_A1 = 'Hillsport dyke (Biscotasing swarm)'
REF_LAT, REF_LON = 48.0, 279.0      # common reference locality of the dyke means, 48.0 N, 81.0 W

AGE_SWARM = dict(age=1141, age_low=1139, age_high=1143)
AGE_A6 = dict(age=1140.6, age_sigma=1.0)
AGE_A3 = dict(age=1142, age_low=1139, age_high=1145)
AGE_A1 = dict(age=2170.7, age_sigma=0.55)
AGE_MATACHEWAN = dict(age_low=2443, age_high=2489)
AGE_TEXT = {
    'swarm': ('Age of the rock: the Abitibi dyke swarm, U-Pb baddeleyite 1140.6 +/- 2 Ma (95%; '
              'Krogh et al., 1987) and 1142 +/- 3 Ma (dyke A3), given as age 1141 Ma, '
              'age_low/age_high 1139-1143 Ma.'),
    'A6': ('Age of the rock: U-Pb baddeleyite 1140.6 +/- 2 Ma (95% confidence; age_sigma is 1 '
           'sigma), Great Abitibi dyke at site 17 (Krogh et al., 1987, Geological Association of '
           'Canada Special Paper 34, 147-152).'),
    'A3': ('Age of the rock: U-Pb baddeleyite 1142 +/- 3 Ma, Wawa Abitibi dyke (L. M. Heaman, '
           'unpublished, cited in Sage, 2000; listed in Piispa et al., 2018, Table 1); the '
           'confidence level is not stated, so it is given as age_low/age_high 1139-1145 Ma.'),
    'A1': ('Age of the rock: U-Pb baddeleyite 2170.7 +/- 1.1 Ma (95% confidence; age_sigma is 1 '
           'sigma), Hillsport dyke, site DP4 of Halls et al. (2005).'),
    'Matachewan': ('Age of the rock: the Matachewan dyke swarm, U-Pb 2446 +/- 3 and 2473 +16/-9 '
                   'Ma (Heaman, 1997), given as age_low/age_high 2443-2489 Ma.'),
}

# ---- Ernst & Buchan (1993) Table 2 -------------------------------------------
# site no., dyke, N_A, n_A, N_T, n_T, lat, long (W), D, I, D', I', k, a95
# (None = not tabulated). 'B'/'U' rows are baked/unbaked Matachewan host rock.
T2 = [
    ('01', 'A1', 10, 4, 8, 0, 49.459, 85.625, 95.3, -49.0, 99.3, -45.2, 30, 17.1),
    ('02', 'A2', 15, 0, 1, 0, 49.022, 84.772, None, None, None, None, None, None),
    ('03', 'A3', 7, 0, 7, 6, 48.205, 84.845, 305.0, 58.5, 306.8, 56.4, 230, 4.4),
    ('04', 'A3', 7, 0, 7, 6, 48.233, 84.817, 288.7, 50.8, 291.2, 48.0, 98, 6.8),
    ('05', 'A3', 9, 0, 10, 0, 50.145, 82.211, None, None, None, None, None, None),
    ('06', 'A4', 7, 7, 1, 0, 48.387, 81.298, 124.7, -71.4, 125.2, -71.2, 129, 5.3),
    ('07', 'A4', 7, 7, 1, 0, 48.533, 80.937, 127.1, -72.0, 127.8, -71.8, 154, 4.9),
    ('08', 'A5', 5, 0, None, None, 48.079, 82.290, None, None, None, None, None, None),
    ('09', 'A5', 5, 5, 1, 0, 48.218, 81.832, 313.6, 49.1, 314.1, 48.5, 131, 6.7),
    ('10', 'A5', 6, 5, None, None, 48.415, 81.202, 316.7, 66.5, 317.1, 66.2, 42, 11.9),
    ('11', 'A5', 10, 8, None, None, 48.429, 81.175, 304.3, 54.6, 304.6, 54.3, 55, 7.6),
    ('12', 'A5', 7, 7, None, None, 48.441, 81.117, 306.0, 67.2, 306.5, 66.9, 42, 9.5),
    ('13', 'A6_1', 10, 8, None, None, 47.878, 82.620, 275.0, 65.2, 275.9, 64.4, 115, 5.2),
    ('14', 'A6_1', 9, 9, None, None, 47.978, 82.325, 292.4, 53.0, 293.2, 52.2, 59, 6.8),
    ('15', 'A6_1', 12, 12, None, None, 48.002, 82.250, 288.5, 57.5, 289.2, 56.7, 82, 4.8),
    ('16', 'A6_2', 1, 1, None, None, 48.072, 82.052, 311.3, 63.5, 311.7, 63.1, None, None),
    ('16B', 'A6_2', 4, 4, None, None, 48.072, 82.052, 310.0, 66.3, 310.3, 65.8, 125, 8.3),
    ('16U', 'A6_2', 3, 3, None, None, 48.072, 82.052, 218.6, -15.7, 219.3, -16.4, 47, 18.2),
    ('17', 'A6_2', 19, 19, 7, 0, 48.185, 81.583, 297.0, 70.7, 297.4, 70.4, 100, 3.4),
    ('18', 'A6_2', 5, 5, None, None, 48.263, 81.333, 299.2, 64.9, 299.6, 64.6, 518, 3.4),
    ('19', 'A6_2', 10, 10, None, None, 48.530, 80.640, 285.5, 69.5, 286.0, 69.5, 67, 5.9),
    ('20', 'A6_2', 18, 17, None, None, 48.550, 80.593, 316.4, 72.7, 317.1, 72.5, 146, 3.0),
    ('21B', 'A6_2', 4, 4, None, None, 48.553, 80.575, 295.3, 75.6, 296.2, 75.6, 3333, 1.6),
    ('21U', 'A6_2', 4, 4, None, None, 48.553, 80.575, 212.9, -9.2, 212.6, -8.0, 150, 7.5),
    ('22', 'A6_2', 5, 4, None, None, 48.562, 80.553, 328.8, 71.6, 329.3, 71.3, 169, 7.1),
    ('23', 'A6_2', 5, 5, None, None, 48.888, 79.505, 300.2, 75.3, 301.5, 75.5, 244, 4.9),
    ('24', 'A6_2', 17, 16, None, None, 48.918, 79.418, 292.8, 68.7, 293.2, 69.2, 39, 6.0),
    ('25', 'A6_1', 4, 4, None, None, 49.585, 77.102, 276.0, 51.1, 274.3, 53.7, 277, 5.5),
    ('26', 'A6_1', 6, 6, None, None, 49.880, 75.780, 283.1, 55.3, 281.1, 58.1, 547, 2.9),
    ('27', 'A7', 5, 5, None, None, 47.075, 83.745, 146.3, -78.1, 143.5, -78.0, 223, 5.1),
    ('28', 'A7', 7, 7, 1, 0, 47.075, 83.745, 157.9, -83.4, 151.1, -83.4, 345, 3.3),
    ('29', 'A8', 5, 4, None, None, 48.151, 79.148, 268.9, 64.6, 267.7, 65.5, 288, 5.4),
    ('30', 'A8', 5, 4, None, None, 48.214, 79.977, 277.5, 51.4, 276.9, 52.1, 95, 9.5),
    ('31', 'A8', 5, 4, None, None, 48.616, 78.221, 261.5, 66.1, 259.8, 67.5, 88, 9.8),
    ('32', 'A8', 5, 4, None, None, 48.676, 78.126, 258.5, 46.0, 256.5, 48.4, 24, 19.1),
    ('33', 'A8', 5, 3, 1, 0, 49.32, 77.18, 264.5, 63.5, 262.8, 65.6, 192, 8.9),
]
THICK = {'01': '~90', '02': '>15', '03': '~50', '04': '>30', '05': '~50', '06': '>35',
         '10': '>30', '12': '>80', '13': '89', '14': '149', '15': '138', '16': '109',
         '17': '210', '18': '231', '19': '186', '20': '187', '21B': '143', '22': '161',
         '23': '185', '24': '220', '25': '148', '26': '85', '27': '80', '28': '80',
         '29': '~165', '30': '~130', '31': '~140', '32': '~85'}
QUEBEC_EAST_OF = 360 - 79.52   # Ontario-Quebec boundary north of Lake Timiskaming

# ---- Ernst & Buchan (1993) Table 3: dyke means at 48.0 N, 81.0 W ----------------
# dyke, N measured, n used, D', I', k, a95, pole lat (N), pole long (W; reversed dykes
# at the normal antipode), dm, dp
T3 = [
    ('A1', 1, 1, 99.3, -45.2, None, None, 25.5, 158.6, None, None),
    ('A3', 3, 2, 298.3, 52.4, None, None, 42.1, 165.6, None, None),
    ('A4', 2, 2, 126.4, -71.5, None, None, 57.0, 136.2, None, None),
    ('A5', 5, 4, 310.5, 59.1, 72, 10.9, 54.1, 165.3, 16.3, 12.2),
    ('A6_1', 5, 5, 283.1, 57.2, 153, 6.2, 35.2, 151.2, 9.1, 6.6),
    ('A6_2', 9, 7, 302.7, 70.9, 189, 4.4, 54.8, 137.1, 7.7, 6.7),
    ('A7', 2, 2, 146.2, -80.7, None, None, 61.5, 102.3, None, None),
    ('A8', 5, 5, 264.9, 60.1, 68, 9.4, 26.3, 138.0, 14.2, 10.7),
]
NOT_IN_DYKE_MEAN = {'02', '05', '08', '16', '16B', '16U', '21B', '21U'}

# ---- Halls et al. (2005) Table 4: SL13 ---------------------------------------------
# N, SLAT, SLON(W), strike, widths (m), DEC, INC, PLAT, PLON(E), k, a95
SL13 = (14, 48.35, 84.67, '020', '5 and 8', 315.9, 62.7, 59.6, 194.1, 61, 5.1)

UNIT_NAME = {'A1': 'dyke A1 (Hillsport dyke)', 'A2': 'dyke A2', 'A3': 'dyke A3',
             'A4': 'dyke A4', 'A5': 'dyke A5', 'A6_1': 'unit 1 of the Great Abitibi dyke (A6)',
             'A6_2': 'unit 2 of the Great Abitibi dyke (A6)', 'A7': 'dyke A7', 'A8': 'dyke A8'}
A3_OVERPRINT = ('Thermal demagnetization of the samples of sites 03 and 04 resolves, below the '
                'characteristic component (510-575 C), a reversed (RK, 475-510 C) and a normal '
                '(NK, 250-475 C) component that the authors interpret as thermal overprints from '
                'the Keweenawan Midcontinent Rift; their poles are RK 49.9 N, 201.8 E (dp 6.5, dm '
                '8.3, 5 samples) and NK 35.7 N, 181.2 E (dp 3.2, dm 5.2, 12 samples).')


def fisher_a95(k, n):
    R = n - (n - 1) / k
    return np.degrees(np.arccos(1 - (n - R) / R * (20 ** (1 / (n - 1)) - 1)))


def angle(lat1, lon1, lat2, lon2):
    return pmag.angle([lon1, lat1], [lon2, lat2])[0]


def site_lith(site, dyke):
    if dyke == 'A1':
        return 'Monzogabbro'
    if dyke == 'A6_2' or site in ('11', '12'):
        return 'Olivine Gabbro:Monzogabbro'
    return 'Olivine Gabbro'


def dyke_age(dyke):
    key = 'A1' if dyke == 'A1' else 'A6' if dyke.startswith('A6') else 'A3' if dyke == 'A3' \
        else 'swarm'
    return {'A1': AGE_A1, 'A6': AGE_A6, 'A3': AGE_A3, 'swarm': AGE_SWARM}[key], AGE_TEXT[key]


for t in T2:   # transcription check: alpha95 from k and n
    site, dyke, NA, nA, NT, nT, lat, lonw, D, I, Dp, Ip, k, a95 = t
    n = (nA or 0) + (nT or 0)
    if k is not None and n > 1:
        assert abs(fisher_a95(k, n) - a95) <= 0.15, site

GEOGRAPHIC = 'Geographic coordinates (vertical dykes; no structural correction).'


def site_record(t):
    """One Table 2 sampling site as text for the description of its dyke."""
    site, dyke, NA, nA, NT, nT, lat, lonw, D, I, Dp, Ip, k, a95 = t
    n = (nA or 0) + (nT or 0)
    s = f'site {site} ({lat} N, {lonw} W' + (f'; dyke thickness {THICK[site]} m' if site in THICK
                                              else '')
    s += '; in Quebec' if 360 - lonw > QUEBEC_EAST_OF else ''
    s += ')'
    if D is None:
        s += ': no stable direction'
        if site == '05':
            s += (', on the probable extension of the dyke (overprint and characteristic '
                  'components could not be separated)')
    else:
        s += f': D {D}, I {I} (D\' {Dp}, I\' {Ip})'
        s += f', k {k}, alpha95 {a95}' if k is not None else ''
        s += f', n {n}'
        notes = []
        if site in ('03', '04'):
            notes.append('thermal demagnetization')
        if site == '16':
            notes.append('single sample from the chilled margin, not used in the dyke mean')
        if site in ('09', '17'):
            notes.append('weak unresolved secondary Keweenawan component')
        if site == '17':
            notes.append('locality of the U-Pb date of Krogh et al., 1987')
        if site == '33':
            notes.append('coordinates approximate in the paper')
        if notes:
            s += ' (' + '; '.join(notes) + ')'
    s += f'; samples measured {NA} by alternating field' + (f', {NT} thermally' if NT else '')
    return s


rows = []
for t in T2:   # host rocks of the baked-contact test
    site, dyke, NA, nA, NT, nT, lat, lonw, D, I, Dp, Ip, k, a95 = t
    if site[-1] not in 'BU':
        continue
    num = site.rstrip('BU')
    baked = site.endswith('B')
    meas = f'Samples measured: {NA} by alternating field.'
    if site == '16B':
        text = ('Baked host rock: north-trending Matachewan diabase dyke at the contact of the '
                'Great Abitibi dyke (A6) at site 16; carries the unit 2 Abitibi direction, as '
                'does the single sample from the chilled margin of the Abitibi dyke; positive '
                'baked-contact test (Ernst and Buchan, 1993).')
    elif site == '21B':
        text = ('Baked host rock: Matachewan diabase dyke 50 m from the contact of the Great '
                'Abitibi dyke (A6) at site 21, where the dyke does not outcrop and is located '
                'by its aeromagnetic signature; carries the unit 2 Abitibi direction. No dyke '
                'material was sampled at this locality, so it is not a formal baked-contact '
                'test, but it supports a primary unit 2 remanence (Ernst and Buchan, 1993).')
    else:
        dist = '240' if site == '16U' else '370'
        text = (f'Unbaked host rock: Matachewan diabase dyke {dist} m from the contact of the '
                f'Great Abitibi dyke (A6) at site {num}; carries the shallow Matachewan '
                f'direction (Ernst and Buchan, 1993).')
    text += f' {meas} Not part of any pole.'
    if baked:
        text += (' The remanence of the baked rock dates from the intrusion of the Abitibi '
                 'dyke (ca. 1141 Ma).')
    text += f' {GEOGRAPHIC} {AGE_TEXT["Matachewan"]}'
    rows.append(dict(site=f'{dyke}-{int(num)}{site[len(num):]}', site_alternatives=site,
                     location=LOC, result_type='i', result_quality='b',
                     method_codes='LP-DIR-AF:DE-BFL:DE-FM:ST-C', citations='This study',
                     geologic_classes='Intrusive',
                     geologic_types='Baked Contact' if baked else 'Volcanic Dike',
                     lithologies='Diabase', lat=lat, lon=round(360 - lonw, 3), age_unit='Ma',
                     **AGE_MATACHEWAN, dir_tilt_correction=0, dir_dec=D, dir_inc=I,
                     dir_alpha95=a95, dir_k=k, dir_n_samples=(nA or 0) + (nT or 0),
                     description=text))

# ---- dyke means (Table 3) -----------------------------------------------------------------
for (dyke, N, n, Dp, Ip, k, a95, plat, plonw, dm, dp) in T3:
    members = [t for t in T2 if t[1] == dyke and t[0] not in NOT_IN_DYKE_MEAN and t[8] is not None]
    assert len(members) == n, dyke
    if n > 1:
        fm = pmag.fisher_mean([[t[10], t[11]] for t in members])
        assert angle(fm['inc'], fm['dec'], Ip, Dp) < 0.3, dyke
        if k is not None:
            assert abs(fm['k'] - k) / k < 0.02 and abs(fm['alpha95'] - a95) < 0.15, dyke
    clon, clat, _, _ = pmag.dia_vgp(Dp, Ip, 0, REF_LAT, REF_LON)
    tlat, tlon = (-plat, (360 - plonw + 180) % 360) if Ip < 0 else (plat, (360 - plonw) % 360)
    assert angle(clat, clon, tlat, tlon) < 0.1, dyke
    age, age_text = dyke_age(dyke)
    sites_used = ', '.join(t[0] for t in members)
    text = (f'Mean of {n} site mean directions for {UNIT_NAME[dyke]} (sites {sites_used}; {N} '
            f'sites sampled), recalculated by the authors to the common reference locality '
            f'48.0 N, 81.0 W (Ernst and Buchan, 1993, Table 3); the coordinates are those of the '
            f'reference locality. vgp_lat/vgp_lon are calculated from the mean direction at that '
            f'locality; Table 3 lists the pole as {plat} N, {plonw} W'
            + (' (for the normal-polarity direction)' if Ip < 0 else '')
            + (f', dp {dp}, dm {dm}' if dp else '') + '.')
    if dyke == 'A6_1':
        text += (' Unit 1 is the earlier, geochemically distinct pulse of the Great Abitibi dyke; '
                 'its direction is from sites 13-15, 25 and 26, and unit 1 samples at sites 16-24 '
                 'carry the unit 2 direction, which the authors interpret as an overprint by the '
                 'intrusion of unit 2. The two units are treated as distinct cooling units.')
    if dyke == 'A6_2':
        text += (' Unit 2 is the later interior pulse of the Great Abitibi dyke; its primary '
                 'remanence is supported by a positive baked-contact test at site 16.')
    if dyke == 'A3':
        text += ' ' + A3_OVERPRINT
    if dyke == 'A1':
        text += (' Dyke A1 is the 2170.7 +/- 1.1 Ma Hillsport dyke of the Biscotasing swarm '
                 '(Halls et al., 2005). The direction is from four reversely magnetized '
                 'samples; of the others, two may carry normal Abitibi directions and four are '
                 'scattered, and the authors advise caution in its use.')
    dyke_sites = [t for t in T2 if t[1] == dyke and t[0][-1] not in 'BU']
    text += (' Sampling sites (Ernst and Buchan, 1993, Table 2; D/I in situ, D\'/I\' at the '
             'reference locality, n samples used): ' + '; '.join(site_record(t) for t in dyke_sites)
             + '.')
    text += f' {GEOGRAPHIC} {age_text}'
    rows.append(dict(
        site=dyke, location=LOC_A1 if dyke == 'A1' else LOC, result_type='a', result_quality='g',
        method_codes=('LP-DIR-T' if dyke == 'A3' else 'LP-DIR-AF') + ':DE-BFL:DE-FM',
        citations='This study' + (f':{HALLS}' if dyke == 'A1' else ''),
        geologic_classes='Intrusive', geologic_types='Volcanic Dike',
        lithologies=site_lith('', dyke), lat=REF_LAT, lon=REF_LON, age_unit='Ma', **age,
        dir_tilt_correction=0, dir_comp_name='characteristic', dir_dec=Dp, dir_inc=Ip,
        dir_alpha95=a95, dir_k=k, dir_n_samples=sum((t[3] or 0) + (t[5] or 0) for t in members),
        dir_polarity='n' if Ip > 0 else 'r', vgp_lat=round(clat, 2), vgp_lon=round(clon, 2),
        vgp_dp=dp, vgp_dm=dm, description=text))

# ---- SL13 (Halls et al., 2005, Table 4) ---------------------------------------------------
N, slat, slonw, strike, widths, D, I, plat, plon, k, a95 = SL13
sl_lon = round(360 - slonw, 2)
clon, clat, cdp, cdm = pmag.dia_vgp(D, I, a95, slat, sl_lon)
assert angle(clat, clon, plat, plon) < 0.1 and abs(fisher_a95(k, N) - a95) < 0.15
rows.append(dict(
    site='SL13', location=LOC, result_type='i', result_quality='g',
    method_codes='SO-SM:LP-DIR-AF:DE-BFL:DE-FM:ST-C', citations=HALLS,
    geologic_classes='Intrusive', geologic_types='Volcanic Dike', lithologies='Mafic Dike',
    lat=slat, lon=sl_lon, age_unit='Ma', **AGE_SWARM, dir_tilt_correction=0,
    dir_comp_name='characteristic', dir_dec=D, dir_inc=I, dir_alpha95=a95, dir_k=k,
    dir_n_samples=N, dir_polarity='n', vgp_lat=round(clat, 2), vgp_lon=round(clon, 2),
    vgp_dp=round(cdp, 1), vgp_dm=round(cdm, 1),
    description=(
        f'Site SL13 of Halls et al. (2005, Table 4): a younger set of mafic dykes (widths '
        f'{widths} m, strike {strike}) cutting a Matachewan dyke, listed by the authors under '
        f'Abitibi; a positive baked-contact test shows that the remanence of the younger dykes is '
        f'primary. The authors consider the swarm affinity conjectural: the geochemistry is '
        f'anomalous and the dykes may be apophyses of a large Abitibi dyke that passes the site '
        f'immediately to the south. Alternating-field demagnetization with principal component '
        f'analysis; sun and magnetic compass. vgp_lat/vgp_lon are calculated from the site '
        f'direction and coordinates; Table 4 lists the VGP as {plat} N, {plon} E. {GEOGRAPHIC} '
        f'{AGE_TEXT["swarm"]}')))

sites = pd.DataFrame(rows)
for c in ('dir_n_samples', 'dir_tilt_correction'):
    sites[c] = sites[c].astype('Int64')
assert sites['site'].is_unique

# ---- locations ----------------------------------------------------------------------------
POLE_UNITS = ['A3', 'A4', 'A5', 'A6_1', 'A6_2', 'A7', 'A8', 'SL13']
units = sites.set_index('site').loc[POLE_UNITS].reset_index()
blk = [[(lo + 180) % 360, -la] if la < 0 else [lo, la]
       for la, lo in zip(units['vgp_lat'], units['vgp_lon'])]
p = pmag.fisher_mean(blk)
normal = units['dir_inc'] > 0
d = pmag.fisher_mean([[dd, ii] if ii > 0 else [(dd + 180) % 360, -ii]
                      for dd, ii in zip(units['dir_dec'], units['dir_inc'])])
coords = pd.DataFrame([dict(dyke=t[1], lat=t[6], lon=round(360 - t[7], 3),
                            used=t[0] not in NOT_IN_DYKE_MEAN and t[8] is not None) for t in T2]
                      + [dict(dyke='SL13', lat=slat, lon=sl_lon, used=True)])
abitibi = coords[coords['dyke'] != 'A1']                      # every Abitibi sampling site
contrib = coords[coords['dyke'].isin(POLE_UNITS) & coords['used']]   # sites in the pole
GEO = dict(location_type='Region', continent_ocean='North America', country='Canada',
           terranes='Laurentia:Superior Province', age_unit='Ma')


def r1(x):
    return round(float(x), 1)


locations = pd.DataFrame([
    dict(location=LOC, **GEO, sites=':'.join(sites.loc[sites['location'] == LOC, 'site']),
         lat_s=abitibi['lat'].min(), lat_n=abitibi['lat'].max(), lon_w=abitibi['lon'].min(),
         lon_e=abitibi['lon'].max(), age_low=1139, age_high=2489,
         geologic_classes='Intrusive',
         lithologies=':'.join(sorted({x for s in sites.loc[sites['location'] == LOC, 'lithologies']
                                      for x in s.split(':')})),
         method_codes='LP-DIR-AF:LP-DIR-T', citations=f'This study:{HALLS}', contact_test='C+',
         description=(
             'Northeast-trending Abitibi dyke swarm of olivine gabbro to olivine-bearing '
             'monzogabbro dykes, southern Superior Province, Ontario and Quebec (Ernst and Buchan, '
             '1993). The site rows are the cooling units: the dyke means of Table 3 (dykes A3, A4, '
             'A5, A7, A8 and units 1 and 2 of the Great Abitibi dyke, A6), each listing its '
             'sampling sites in its description, and site SL13 of Halls et al. (2005); with the '
             'baked and unbaked Matachewan host rocks of the baked-contact test at sites 16 and '
             '21. Dyke A2 (site 02, 49.022 N, 84.772 W) gave no stable direction and has no dyke '
             'mean. Dyke A1 of Ernst and Buchan (1993), a Biscotasing '
             'dyke, is in the location Hillsport dyke (Biscotasing swarm). Alternating-field '
             'demagnetization (about 15 steps to 100 mT) of two specimens per sample and thermal '
             'demagnetization at selected sites, with least-squares line fits; the remanence is '
             'carried by fine-grained low-Ti titanomagnetite. The age range spans the Abitibi '
             'dykes and the Matachewan host rocks. Sites in Ontario and Quebec.')),
    dict(location=LOC, **GEO, result_name='Abitibi dyke swarm pole',
         result_type='a', result_quality='g', sites=':'.join(POLE_UNITS),
         lat_s=contrib['lat'].min(), lat_n=contrib['lat'].max(), lon_w=contrib['lon'].min(),
         lon_e=contrib['lon'].max(), **AGE_SWARM, geologic_classes='Intrusive',
         lithologies='Mafic Dike:Monzogabbro:Olivine Gabbro',
         method_codes='LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP',
         citations=f'This study:{HALLS}:{PIISPA}', contact_test='C+',
         dir_tilt_correction=0, dir_dec=r1(d['dec']), dir_inc=r1(d['inc']),
         dir_alpha95=r1(d['alpha95']), dir_k=r1(d['k']), dir_n_sites=int(d['n']),
         pole_comp_name='characteristic', pole_lat=r1(p['inc']), pole_lon=r1(p['dec']),
         pole_alpha95=r1(p['alpha95']), pole_k=r1(p['k']), pole_n_sites=int(p['n']),
         pole_reversed_perc=round(100 * (~normal).sum() / len(units), 1),
         description=(
             'Fisher mean of eight cooling-unit VGPs after polarity unification: the dyke means '
             'of dykes A3, A4, A5, A7 and A8 and of units 1 and 2 of the Great Abitibi dyke (A6) '
             '(Ernst and Buchan, 1993, Table 3, at the common reference locality 48.0 N, 81.0 W), '
             'and site SL13 of Halls et al. (2005), whose swarm affinity those authors consider '
             'conjectural. Dyke A1 is excluded as a 2170.7 Ma Biscotasing dyke (Halls et al., '
             '2005). The same eight units form the modified Abitibi mean of Piispa et al. (2018). '
             'The normal-polarity unit 2 remanence of dyke A6 is primary on the basis of a '
             'positive baked-contact test against a Matachewan dyke at site 16. Two dykes (A4, '
             'A7) are reversed; their directions are steeper than the normal mean, and Ernst and '
             'Buchan (1993) regarded the reversed data as too few to average secular variation '
             'or to test reversal symmetry. Ernst and Buchan (1993) also give a five-dyke '
             'normal-polarity pole (42.8 N, 151.5 W, dp 12.5, dm 16.3) and an all-dyke pole that '
             'includes A1 (44.4 N, 148.6 W, dp 12.1, dm 15.4; GPMDB 7193, MagIC legacy '
             'contribution 17341), both poles of mean directions. '
             'The mean direction is given in the majority (normal) polarity and combines the dyke '
             'means at the reference locality with SL13 at its own site. Age: U-Pb baddeleyite '
             '1140.6 +/- 2 Ma (95%; Krogh et al., 1987) of the Great Abitibi dyke and 1142 +/- 3 '
             'Ma of dyke A3. Sites in Ontario and Quebec.')),
    dict(location=LOC_A1, **GEO, state_province='Ontario', sites='A1', lat_s=49.459,
         lat_n=49.459, lon_w=round(360 - 85.625, 3), lon_e=round(360 - 85.625, 3), **AGE_A1,
         geologic_classes='Intrusive', lithologies='Monzogabbro', method_codes='LP-DIR-AF',
         citations=f'This study:{HALLS}',
         description=(
             'Dyke A1 of Ernst and Buchan (1993), sampled at site 01 as part of the Abitibi swarm. '
             'Halls et al. (2005) dated the same dyke (Hillsport dyke, their site DP4) by U-Pb '
             'baddeleyite at 2170.7 +/- 1.1 Ma (95%) and assigned it to the Biscotasing swarm, so '
             'its reversed direction is not part of the Abitibi pole. The Ernst and Buchan (1993) '
             'direction rests on four of ten samples, and the authors advise caution in its use.')),
])
for c in ('dir_n_sites', 'pole_n_sites', 'dir_tilt_correction'):
    locations[c] = locations[c].astype('Int64')

# ---- ages table: the U-Pb dates (Krogh et al. 1987 and Halls et al. 2005, 95%; dyke A3) -------
PER = dict(age_unit='Ma', timescale_eon='Proterozoic')
age_rows = [
    dict(location=LOC, site='A6_2', age=1140.6, age_sigma=1.0, method_codes='GM-UPB',
         citations=KROGH, timescale_era='Mesoproterozoic', timescale_period='Stenian',
         description='U-Pb baddeleyite 1140.6 +/- 2 Ma (95%; age_sigma is 1 sigma), average '
                     '207Pb/206Pb age of the baddeleyite fractions, Great Abitibi dyke (A6) on '
                     'Highway 144 south of Timmins, site 17 (Krogh et al., 1987, Geological '
                     'Association of Canada Special Paper 34, 147-152).'),
    dict(location=LOC, site='A3', age=1142, age_low=1139, age_high=1145, method_codes='GM-UPB',
         citations=PIISPA, timescale_era='Mesoproterozoic', timescale_period='Stenian',
         description='U-Pb baddeleyite 1142 +/- 3 Ma, Wawa Abitibi dyke (dyke A3; L. M. Heaman, '
                     'unpublished, cited in Sage, 2000, Ontario Geological Survey Open File Report '
                     '6106; listed in Piispa et al., 2018, Table 1); confidence level not stated.'),
    dict(location=LOC_A1, site='A1', age=2170.7, age_sigma=0.55, method_codes='GM-UPB',
         citations=HALLS, timescale_era='Paleoproterozoic', timescale_period='Rhyacian',
         description='U-Pb baddeleyite 2170.7 +/- 1.1 Ma (95%; age_sigma is 1 sigma), Hillsport '
                     'dyke, site DP4 of Halls et al. (2005) = dyke A1, site 01 of Ernst and '
                     'Buchan (1993).'),
]
ages = pd.DataFrame(age_rows)
for k_, v_ in PER.items():
    ages[k_] = v_


# ---- write ----------------------------------------------------------------------------------
def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False)


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locations, 'locations', os.path.join(OUT, 'locations.txt'))
write_magic(ages, 'ages', os.path.join(OUT, 'ages.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Ernst1993_Abitibi_dykes_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    for i, (kind, df) in enumerate([('locations', locations), ('sites', sites), ('ages', ages)]):
        if i:
            f.write('>>>>>>>>>>\n')
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False)

print(f'{len(sites)} site rows; pole {p["inc"]:.2f}N {p["dec"]:.2f}E A95 {p["alpha95"]:.2f} '
      f'K {p["k"]:.1f} N {p["n"]}; dir {d["dec"]:.1f}/{d["inc"]:.1f} a95 {d["alpha95"]:.1f} '
      f'k {d["k"]:.1f}; samples in pole {int(units["dir_n_samples"].sum())}')
print(f'wrote {combined}')

from validate_magic_contribution import validate_upload_file  # noqa: E402
if not validate_upload_file(combined, tables=['locations', 'sites', 'ages']):
    raise SystemExit('validation failed')
