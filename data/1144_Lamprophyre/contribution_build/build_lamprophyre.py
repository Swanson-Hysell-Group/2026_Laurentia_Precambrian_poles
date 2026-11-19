"""Build the MagIC contribution for the ca. 1144 Ma lamprophyre dykes of Piispa et al. (2018).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Piispa, E. J., Smirnov, A. V., Pesonen, L. J., & Mitchell, R. H. (2018). Paleomagnetism
and geochemistry of ~1144-Ma lamprophyre dikes, northwestern Ontario: Implications for
the North American polar wander and plate velocities. J. Geophys. Res. Solid Earth 123,
6195-6214. doi:10.1029/2018JB015992

Dates: Queen, M., Heaman, L. M., Hanes, J. A., Archibald, D. A., & Farrar, E. (1996).
40Ar/39Ar phlogopite and U-Pb perovskite dating of lamprophyre dykes from the eastern
Lake Superior region: evidence for a 1.14 Ga magmatic precursor to Midcontinent Rift
volcanism. Can. J. Earth Sci. 33, 958-965. doi:10.1139/e96-072. U-Pb errors are 2 sigma
(p. 960) and the 40Ar/39Ar plateau errors 2 sigma (Fig. 2); age_sigma is 1 sigma.

Sources and checks
------------------
- Table 3 (site means, host-rock directions of the baked-contact tests, the Abitibi dyke
  means and all group means and poles) was transcribed twice, independently, from the
  rendered page; the two transcriptions agree on every coordinate, count, direction,
  alpha95 and k.
- Lamprophyre site VGPs are calculated from each site's direction and coordinates and
  stored to 0.01 deg, in the polarity of the site direction; the VGP listed in Table 3
  (to 0.1 deg, for the normal-polarity direction) is given in each site description.
- The VGPs listed for the Abitibi dyke means were calculated by the authors at
  coordinates other than the approximate ones in Table 3, so they are stored as listed
  (in the polarity of the direction).
- With these VGPs both poles are reproduced at their printed precision (asserted
  below); the mean directions are recalculated from the site directions.
- The reversal tests of the paper are reproduced from the site means: lamprophyre R vs N
  class C, combined lamprophyre and Abitibi R vs N class B.

Structure and interpretation
----------------------------
- One location: an umbrella row and the two poles of the paper that combine both
  polarities, Lc (lamprophyre dykes) and LAc (lamprophyre and Abitibi dykes). Pole Lc is
  the pole used in the Swanson-Hysell et al. (2026) compilation. The single-polarity and
  Abitibi-only means of Table 3 are not entered as poles; the reversal tests use the
  site data.
- Site and pole ages are ages of the rocks. Dated lamprophyre dykes carry their U-Pb
  perovskite dates (D9, D15); the others the age of the swarm, 1144 Ma (1137-1151 Ma),
  from the 1144 +/- 7 Ma 40Ar/39Ar phlogopite and 1144 +6/-4 Ma composite U-Pb perovskite
  dates. Abitibi dyke means carry their U-Pb dates where dated (A3, A6) and otherwise
  1141 Ma (1139-1143 Ma). The host rocks of the baked-contact tests carry their own ages:
  the Matachewan dyke swarm (2443-2489 Ma) and the Archean greenstone belt.
"""
import contextlib
import io
import os
import sys
from datetime import date

import matplotlib
matplotlib.use('Agg')  # ipmag.common_mean_watson calls plt.show(), which blocks on a GUI backend
import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

LOCATION = 'Lamprophyre and Abitibi dykes, northeast of Lake Superior'
QUEEN = '10.1139/e96-072'
EB = '10.1139/e93-167'
KROGH = 'Krogh et al. (1987) Geological Association of Canada Special Paper 34, 147-152'
HALLS = ('Halls, H. C., Stott, G. M., & Davis, D. W. (2005). Paleomagnetism, geochronology '
         'and geochemistry of several Proterozoic mafic dike swarms in northwestern Ontario. '
         'Ontario Geological Survey Open File Report 6171')

AGE_LAMP = dict(age=1144, age_low=1137, age_high=1151)
AGE_ABI = dict(age=1141, age_low=1139, age_high=1143)
AGE_COMB = dict(age=1143, age_low=1137, age_high=1151)

# Table 3, lamprophyre dykes: site as printed, lat (N), lon (W), method, quality class,
# n specimens, N samples, D, I, a95, k, listed VGP lat, lon (for the normal-polarity
# direction), width (m), trend
T3 = """
ML2   48.8112 86.7228 AF+T   A 12  6 103.8 -61.4  5.5 125.0 38.6 206.9 0.3 NNW
D15   48.8027 86.7095 AF+T   A  9  5 126.2 -59.4 10.2  46.2 51.5 191.2 0.8 NW
D16_H 48.8034 86.7066 AF+T   A 11  8 118.8 -57.2  7.6  47.2 45.4 192.9 0.2 NW
MM3   48.8081 86.7076 AF+T   A  6  4 112.3 -85.9  9.9  65.0 51.3 261.2 1.0 E
D14   48.8024 86.7115 AF+T_c B  9  6 124.7 -67.4  9.8  48.0 54.9 207.1 0.4 NW
MM1   48.8029 86.7080 T      B  3  3 106.3 -69.8 12.9  61.6 45.7 218.8 0.2 NW
D10_H 48.8022 86.7145 AF+T_c B  9  7 184.6 -71.4 11.4  45.9 82.2 292.6 1.2 NW
D11   48.8022 86.7145 AF+T_c B 10  6 147.7 -55.6 11.7  25.5 63.3 167.1 0.1 NW
D12   48.8021 86.7141 AF+T_c B 10  6 123.1 -64.9 14.0  19.6 52.7 202.5 0.3 NW
MM2   48.8029 86.7080 AF+T   B  3  2 110.8 -59.7 13.3  58.4 41.8 200.9 0.2 NW
D1    47.7043 83.3191 AF+T_c B 18  8 182.7 -78.1  7.6  46.6 70.5 279.8 0.3 NE
D3_H  47.7041 83.3187 AF+T_c B  5  5 130.0 -68.1 14.2  76.0 57.9 211.9 0.2 NE
D4    48.0584 82.7135 AF+T_c B 18 15 117.9 -80.4 10.7  20.2 53.6 248.8 1.0 NE
D7_H  48.0373 84.5434 AF+T_c B 12  8 102.9 -79.0 14.0  13.7 48.3 243.4 0.7 ENE
D2_H  47.7042 83.3182 AF+T   A  5  5 259.4  81.1  8.6  64.7 42.0 253.4 0.1 NE
D5    48.0590 82.7131 AF+T   A  9  5 331.6  77.2 14.8  11.7 67.0 247.1 0.6 NE
D9    48.0373 84.5434 AF+T   A 10  6 271.2  61.6 14.3  19.0 31.0 216.5 0.2 NE
D13   48.8023 86.7120 AF+T   A 14 11 301.1  81.7 14.6   9.8 54.8 248.7 0.2 NW
ML1   48.8112 86.7228 AF+T   A  8  4   9.5  70.5 11.2  51.1 81.7 314.6 0.3 NNW
"""
# Host rocks of the baked-contact tests (Table 3; no coordinates, VGP or quality class):
# site, row as printed, dyke, n specimens, N samples, D, I, a95, k
HOST = [
    ('MM2 Baked', 'MM2 Baked', 'MM2', 2, 1, 170.7, -82.4, None, None),
    ('MM2 Unbaked', 'MM2 Unbaked', 'MM2', 2, 1, 74.7, 86.0, None, None),
    ('D1 Baked', 'D1 Baked', 'D1', 4, 2, 60.9, -75.6, 7.9, 503.0),
    ('D1 Unbaked', 'DI Unbaked', 'D1', 10, 5, 25.0, 33.0, 3.2, 465.0),
]
# Abitibi dyke means as tabulated in Table 3: site, dyke, approx lat (N), lon (W), method,
# B sites, N samples, D, I, a95, k, listed VGP lat, lon (for the normal-polarity direction)
ABITIBI = [
    ('A4a', 'dyke A4', 48.45, 81.10, 'AF', 2, 14, 125.9, -71.7, None, None, 57.0, 223.7),
    ('A7abi', 'dyke A7', 47.08, 83.75, 'AF', 2, 12, 150.4, -80.8, None, None, 61.5, 257.8),
    ('A3-abi', 'dyke A3', 48.20, 84.80, 'T', 2, 12, 296.1, 54.9, None, None, 42.5, 194.4),
    ('A5-abi', 'dyke A5', 48.40, 81.20, 'AF', 4, 25, 310.1, 59.5, 10.7, 74.0, 54.5, 196.0),
    ('A6-1-a', 'unit 1 of the Great Abitibi dyke (A6)', 47.90, 82.40, 'AF', 5, 39, 283.3,
     56.6, 6.5, 138.8, 35.2, 209.2),
    ('A6-2-a', 'unit 2 of the Great Abitibi dyke (A6)', 48.50, 79.50, 'AF', 7, 76, 302.0,
     71.0, 4.4, 193.3, 54.9, 223.3),
    ('A8-abi', 'dyke A8', 48.50, 78.20, 'AF', 5, 19, 266.3, 58.5, 9.5, 66.4, 26.8, 222.8),
    ('SL13', 'site SL13', 48.35, 84.67, 'AF', 1, 14, 315.9, 62.7, 5.1, 61.0, 59.6, 194.1),
]
QUEBEC = {'A6-2-a', 'A8-abi'}   # approximate coordinates east of the Ontario-Quebec boundary

lamp = pd.DataFrame([line.split() for line in T3.strip().splitlines()],
                    columns=['site', 'lat', 'lon_w', 'method', 'q', 'n_spec', 'n_samp', 'dec',
                             'inc', 'a95', 'k', 'tvgp_lat', 'tvgp_lon', 'width', 'trend'])
lamp['site'] = lamp['site'].str.replace('_', ' ')
lamp['method'] = lamp['method'].str.replace('_', ' ')
for c in lamp.columns.drop(['site', 'method', 'q', 'trend']):
    lamp[c] = pd.to_numeric(lamp[c])
assert len(lamp) == 19 and (lamp['inc'] < 0).sum() == 14


def antipode(lat, lon):
    return -lat, (lon + 180) % 360


def method_codes(method):
    codes = ['SO-SM', 'LP-LT']
    codes += ['LP-DIR-AF'] if 'AF' in method else []
    codes += ['LP-DIR-T'] if 'T' in method.replace('AF', '') else []
    codes += ['DE-BFL', 'DE-BFP', 'DE-FM-LP'] if method.endswith(' c') else ['DE-BFL', 'DE-FM']
    return ':'.join(codes)


CLASS_TEXT = ('quality class {q} of the authors (A: four or more samples, line fits only; B: '
              'partly or wholly from remagnetization-circle fits and/or two or three samples)')
GEOGRAPHIC = ('Geographic coordinates (vertical dykes in the Canadian Shield; no structural '
              'correction).')
SWARM_AGE = ('Age of the rock: the 1144 Ma lamprophyre swarm (Queen et al., 1996): 40Ar/39Ar '
             'phlogopite 1144 +/- 7 Ma and composite U-Pb perovskite 1144 +6/-4 Ma (2 sigma), '
             'given as age 1144 Ma, age_low/age_high 1137-1151 Ma.')

sites = []
for r in lamp.itertuples():
    rev = r.inc < 0
    lon = round(360 - r.lon_w, 4)
    plon, plat, dp, dm = pmag.dia_vgp(r.dec, r.inc, r.a95, r.lat, lon)
    tlat, tlon = antipode(r.tvgp_lat, r.tvgp_lon) if rev else (r.tvgp_lat, r.tvgp_lon)
    assert pmag.angle([plon, plat], [tlon, tlat])[0] < 0.15, r.site
    area = 'East Lake Superior area' if r.lat < 48.5 else 'Marathon area'
    text = [f'Ultramafic lamprophyre dyke, {area}, Ontario; width {r.width} m, trend {r.trend} '
            f'(Table 3).']
    if r.site.endswith(' H'):
        text.append(f"Printed as '{r.site}' in Table 3; in the first author's 2015 Michigan Tech "
                    f"dissertation the 'H' marks the high-coercivity component of dykes that "
                    f"also carry a low-coercivity component.")
    text.append(f'Site mean of the characteristic remanence, {CLASS_TEXT.format(q=r.q)}; '
                f'liquid-nitrogen treatment then demagnetization {r.method.replace(" c", "")}'
                + ('; includes remagnetization circles combined with line fits (McFadden and '
                   'McElhinny, 1988)' if r.method.endswith(' c') else '')
                + f'; n/N = {r.n_spec}/{r.n_samp} specimens/samples.')
    text.append(f'vgp_lat/vgp_lon are calculated from the site direction and coordinates; '
                f'Table 3 lists the VGP as {r.tvgp_lat}N, {r.tvgp_lon}E (for the normal-polarity '
                f'direction).')
    if r.site == 'D1':
        text.append('The 0.3 m wide dyke cuts a ~10 m wide Matachewan dyke; baked-contact test '
                    'positive (host-rock rows D1 Baked and D1 Unbaked).')
    if r.site == 'MM2':
        text.append('Cuts Archean Schreiber-White River greenstone-belt rocks; the authors consider '
                    'this baked-contact test only tentatively positive (one baked and one unbaked '
                    'sample; host-rock rows MM2 Baked and MM2 Unbaked).')
    text.append(GEOGRAPHIC)
    if r.site == 'D9':
        age = dict(age=1143, age_sigma=6)
        text.append('Age of the rock: U-Pb perovskite 1143 +/- 12 Ma (2 sigma; age_sigma is 1 '
                    'sigma), sample LH87-10 from this dyke (Queen et al., 1996).')
    elif r.site == 'D15':
        age = dict(age=1145, age_low=1135, age_high=1160)
        text.append('Age of the rock: U-Pb perovskite 1145 +15/-10 Ma (2 sigma), sample LH85-4B '
                    'from this dyke (Queen et al., 1996), given as age_low/age_high 1135-1160 Ma.')
    else:
        age = AGE_LAMP
        text.append(SWARM_AGE)
    sites.append(dict(
        site=r.site, location=LOCATION, result_type='i', result_quality='g',
        method_codes=method_codes(r.method), citations='This study',
        geologic_classes='Intrusive', geologic_types='Volcanic Dike', lithologies='Lamprophyre',
        lat=r.lat, lon=lon, age_unit='Ma', **age, dir_tilt_correction=0,
        dir_comp_name='characteristic', dir_dec=r.dec, dir_inc=r.inc, dir_alpha95=r.a95,
        dir_k=r.k, dir_n_samples=r.n_samp, dir_n_specimens=r.n_spec,
        dir_polarity='r' if rev else 'n', vgp_lat=round(plat, 2), vgp_lon=round(plon, 2),
        vgp_dp=round(dp, 1), vgp_dm=round(dm, 1), description=' '.join(text)))

coords = {s['site']: (s['lat'], s['lon']) for s in sites}
for site, printed, dyke, n_spec, n_samp, D, I, a95, k in HOST:
    baked = 'Unbaked' not in site
    lat, lon = coords[dyke]
    if dyke == 'D1':
        text = [f"{'Baked' if baked else 'Unbaked'} host rock of a baked-contact test: a ~10 m "
                f"wide Matachewan dyke cut by lamprophyre dyke D1, East Lake Superior area; "
                + ('samples within 5-7 cm of the contact, steep reversed like D1.' if baked else
                   'samples ~20-160 cm from the contact, close to the normal-polarity Matachewan '
                   'dyke mean of Evans and Halls (2010).')]
        cls, typ, lith = 'Intrusive', 'Baked Contact' if baked else 'Volcanic Dike', 'Diabase'
        age = dict(age_low=2443, age_high=2489)
        age_text = ('Age of the rock: the Matachewan dyke swarm, U-Pb 2446 +/- 3 and 2473 +16/-9 '
                    'Ma (Heaman, 1997), given as age_low/age_high 2443-2489 Ma.')
    else:
        text = [f"{'Baked' if baked else 'Unbaked'} host rock of a baked-contact test: Archean "
                f"Schreiber-White River greenstone-belt rock cut by lamprophyre dyke MM2, Marathon "
                f"area; " + ('one sample 5 cm from the contact, broadly similar to MM2.' if baked
                             else 'one sample, different from MM2.')
                + ' The authors consider this test only tentatively positive.']
        cls, typ, lith = 'Metamorphic', 'Baked Contact' if baked else 'Not Specified', 'Greenstone'
        age = dict(age_low=2500, age_high=4031)
        age_text = ('Age of the rock: Archean, given as age_low/age_high 2500-4031 Ma (base of the '
                    'Archean in the ICS International Chronostratigraphic Chart v2026/06).')
    if printed != site:
        text.append(f"Printed as '{printed}' in Table 3.")
    text.append(f'n/N = {n_spec}/{n_samp} specimens/samples; site coordinates are those of dyke '
                f'{dyke}. Not part of any pole.')
    if baked:
        text.append('The remanence of the baked rock dates from the intrusion of the lamprophyre '
                    'dyke (ca. 1144 Ma).')
    text += [GEOGRAPHIC, age_text]
    sites.append(dict(
        site=site, location=LOCATION, result_type='i', result_quality='b',
        method_codes=method_codes('AF + T') + ':ST-C', citations='This study',
        geologic_classes=cls, geologic_types=typ, lithologies=lith, lat=lat, lon=lon,
        age_unit='Ma', **age, dir_tilt_correction=0, dir_dec=D, dir_inc=I, dir_alpha95=a95, dir_k=k,
        dir_n_samples=n_samp, dir_n_specimens=n_spec, description=' '.join(text)))

for site, dyke, lat, lon_w, meth, b, n_samp, D, I, a95, k, tlat, tlon in ABITIBI:
    rev = I < 0
    vlat, vlon = antipode(tlat, tlon) if rev else (tlat, tlon)
    if site == 'SL13':
        text = [f'Abitibi dyke {dyke} of Halls et al. (2005; {HALLS}), a single-site mean, as '
                f'tabulated by Piispa et al. (2018); Halls et al. consider its swarm affinity '
                f'conjectural.']
        cit, lith = 'This study', 'Diabase'
    else:
        text = [f'Abitibi dyke swarm, {dyke}: mean of {b} site means ({n_samp} samples) of Ernst and Buchan '
                f'(1993), as tabulated by Piispa et al. (2018); olivine gabbro dykes (Ernst and '
                f'Buchan, 1993).']
        cit, lith = f'This study:{EB}', 'Olivine Gabbro'
    if a95 is None:
        text.append('Table 3 gives no alpha95 or k for this two-site mean.')
    text.append('Coordinates approximate, derived by Piispa et al. from the site locations of the '
                'original study (Table 3 footnote b); the VGP is as listed in Table 3 (stored in '
                'the polarity of the direction), which the authors calculated at the original '
                'site coordinates.')
    if site in QUEBEC:
        text.append('The approximate coordinates lie in Quebec.')
    text.append(GEOGRAPHIC)
    if site == 'A3-abi':
        age = dict(age=1142, age_low=1139, age_high=1145)
        text.append('Age of the rock: U-Pb baddeleyite 1142 +/- 3 Ma, Wawa Abitibi dyke (L. M. '
                    'Heaman, unpublished, cited in Sage, 2000, Ontario Geological Survey Open '
                    'File Report 6106), as listed by Piispa et al. (2018, Table 1); the confidence '
                    'level is not stated, so it is given as age_low/age_high 1139-1145 Ma.')
    elif site.startswith('A6'):
        age = dict(age=1140.6, age_sigma=1.0)
        text.append('Age of the rock: U-Pb baddeleyite 1140.6 +/- 2 Ma (95% confidence; age_sigma '
                    'is 1 sigma), Great Abitibi dyke (Krogh et al., 1987, Geological Association '
                    'of Canada Special Paper 34, 147-152).')
    else:
        age = AGE_ABI
        text.append('Age of the rock: the Abitibi dyke swarm, U-Pb baddeleyite 1140.6 +/- 2 Ma '
                    '(Krogh et al., 1987) and 1142 +/- 3 Ma, given as age 1141 Ma, '
                    'age_low/age_high 1139-1143 Ma.')
    sites.append(dict(
        site=site, location=LOCATION, result_type='i' if site == 'SL13' else 'a',
        result_quality='g', method_codes=('LP-DIR-AF' if meth == 'AF' else 'LP-DIR-T') + ':DE-BFL:DE-FM',
        citations=cit, geologic_classes='Intrusive', geologic_types='Volcanic Dike',
        lithologies=lith, lat=lat, lon=round(360 - lon_w, 2), age_unit='Ma', **age,
        dir_tilt_correction=0, dir_comp_name='characteristic', dir_dec=D, dir_inc=I,
        dir_alpha95=a95, dir_k=k, dir_n_samples=n_samp, dir_polarity='r' if rev else 'n',
        vgp_lat=vlat, vgp_lon=round(vlon, 1), description=' '.join(text)))

sites = pd.DataFrame(sites)
for c in ('dir_n_samples', 'dir_n_specimens', 'dir_tilt_correction'):
    sites[c] = sites[c].astype('Int64')
assert sites['site'].is_unique and len(sites) == 31


# ---- statistics --------------------------------------------------------------------------
def pole_of(sub):
    """Fisher mean of polarity-unified site VGPs, reported in the northern hemisphere."""
    blk = [[(lo + 180) % 360, -la] if la < 0 else [lo, la]
           for la, lo in zip(sub['vgp_lat'], sub['vgp_lon'])]
    p = pmag.fisher_mean(blk)
    if p['inc'] < 0:
        p['dec'], p['inc'] = (p['dec'] + 180) % 360, -p['inc']
    return p


def direction_of(sub, polarity):
    """Fisher mean of site directions in the given polarity ('n' or 'r')."""
    up = polarity == 'r'
    blk = [[(d + 180) % 360, -i] if (i < 0) != up else [d, i]
           for d, i in zip(sub['dir_dec'], sub['dir_inc'])]
    return pmag.fisher_mean(blk)


def reversal(sub):
    r = sub[sub['dir_inc'] < 0]
    n = sub[sub['dir_inc'] > 0]
    a = [[d, i] for d, i in zip(r['dir_dec'], r['dir_inc'])]
    b = [[(d + 180) % 360, -i] for d, i in zip(n['dir_dec'], n['dir_inc'])]
    with contextlib.redirect_stdout(io.StringIO()):
        out = ipmag.common_mean_watson(a, b, NumSims=10000, plot='no', random_seed=23)
    return out[1], out[2]


lamp_r = [s for s, i in zip(lamp['site'], lamp['inc']) if i < 0]
lamp_n = [s for s, i in zip(lamp['site'], lamp['inc']) if i > 0]
abi_r = [a[0] for a in ABITIBI if a[8] < 0]
abi_n = [a[0] for a in ABITIBI if a[8] > 0]

POLES = [
    dict(name='Pole Lc (lamprophyre dykes)', sites=lamp_r + lamp_n, polarity='n',
         pub=(306.4, 72.1, 5.5, 37.8, 58.0, 223.3, 9.2, 14.2), age=AGE_LAMP, gpmdb=9886,
         tests=dict(contact_test='C+', reversal_test='Rc'),
         text=('The paper\'s pole for the 19 lamprophyre dykes (Table 3, combined mean R + N, 14 '
               'reversed and 5 normal; pole Lc of Table 4, PALEOMAGIA quality factor 6). Reversal '
               'test of the R and N lamprophyre means positive, classification C (gamma_o 8.4, '
               'gamma_c 14.3; McFadden and McElhinny, 1990); reproduced from the site means. '
               'Baked-contact test positive at dyke D1 (Matachewan dyke host); the authors '
               'consider the test at dyke MM2 only tentatively positive. VGP angular dispersion of '
               'the 14 reversed dykes Sb = 18.7 +/- 5.8 (section 4.3). The printed sample total '
               '(121) is one more than the sum of the site N values (120).')),
    dict(name='Pole LAc (lamprophyre and Abitibi dykes)', sites=lamp_r + abi_r + lamp_n
         + abi_n, polarity='n', pub=(303.1, 70.2, 4.5, 38.4, 55.8, 220.0, 7.3, 15.5),
         age=AGE_COMB, gpmdb=9887, pole_vv_q=7,
         tests=dict(contact_test='C+', reversal_test='Rb'),
         text=('The paper\'s grand-mean pole of the 19 lamprophyre dykes and the 8 Abitibi dyke means '
               '(Table 3, combined mean R + N, 27 units; pole LAc of Table 4, quality factor 6), '
               'which the authors interpret as one magmatic pulse or closely spaced events; seven '
               'of seven points on the Van der Voo (1990) scale (section 4.3). Reversal test of the '
               'combined R and N means positive, classification B (gamma_o 3.6, gamma_c 9.9); '
               'reproduced from the site means. Baked-contact tests positive for the lamprophyre '
               'dykes (D1) and for the normal-polarity Abitibi dykes (Ernst and Buchan, 1993). VGP '
               'angular dispersion Sb = 20.0 +/- 4.4. Age: the authors give ~1143 +/- 2 Ma; '
               'age_low/age_high span the dated lamprophyre and Abitibi dykes. The printed sample '
               'total (332) is one more than the sum of the site N values (331). Sites in Ontario '
               'and Quebec.')),
]
TESTS = {'Lc': reversal(sites[sites['site'].isin(POLES[0]['sites'])]),
         'LAc': reversal(sites[sites['site'].isin(POLES[1]['sites'])])}
assert TESTS['Lc'][0] < TESTS['Lc'][1] and TESTS['LAc'][0] < TESTS['LAc'][1], TESTS


def r1(x):
    return round(float(x), 1)


GEO = dict(location=LOCATION, location_type='Region', continent_ocean='North America',
           country='Canada', terranes='Laurentia:Superior Province', age_unit='Ma')
SAMPLING = ('19 of 22 sampled lamprophyre dykes gave accepted site means (3 rejected for '
            'scattered directions and not tabulated: D6, D8 and D17 by elimination of the site '
            'names); 158 oriented samples (drill cores; blocks at D17, MM1 and MM2) oriented by '
            'sun and magnetic compass; liquid-nitrogen treatment, then alternating-field and/or '
            'thermal demagnetization of sister specimens; remanence carried by nearly '
            'stoichiometric magnetite.')
CHECK = []
locations = [dict(
    **GEO, sites=':'.join(sites['site']), lat_s=sites['lat'].min(), lat_n=sites['lat'].max(),
    lon_w=sites['lon'].min(), lon_e=sites['lon'].max(), age_low=1135, age_high=4031,
    geologic_classes='Intrusive:Metamorphic', lithologies=':'.join(sorted(set(sites['lithologies']))),
    method_codes='SO-SM:LP-LT:LP-DIR-AF:LP-DIR-T', citations=f'This study:{QUEEN}:{EB}',
    contact_test='C+',
    description=(
        'Study area of Piispa et al. (2018), Superior Province northeast of Lake Superior, Ontario: '
        '(1) site means of ca. 1144 Ma ultramafic lamprophyre (aillikite) dykes, Marathon area (ML1-2, '
        'MM1-3, D10-D16) and East Lake Superior area (D1-D9); ' + SAMPLING + ' (2) the host rocks of '
        'two baked-contact tests (D1: Paleoproterozoic Matachewan dyke, positive; MM2: Archean '
        'greenstone belt, tentatively positive), whose ages set the upper bound of this row\'s age '
        'range; (3) 8 Abitibi dyke means of Ernst and Buchan (1993) and Halls et al. (2005) as '
        'tabulated by Piispa et al. (2018) with approximate coordinates. Ages of the dykes: U-Pb '
        'perovskite and 40Ar/39Ar phlogopite (Queen et al., 1996) and U-Pb baddeleyite (Krogh et '
        'al., 1987). Sites in Ontario and Quebec.'))]
for p in POLES:
    sub = sites[sites['site'].isin(p['sites'])]
    assert len(sub) == len(p['sites'])
    d = direction_of(sub, 'r' if p['polarity'] == 'r' else 'n')
    q = pole_of(sub)
    pdec, pinc, pa95, pk, plat, plon, A95, K = p['pub']
    assert pmag.angle([d['dec'], d['inc']], [pdec, pinc])[0] < 0.15, (p['name'], d)
    assert (r1(q['inc']), r1(q['dec']), r1(q['alpha95'])) == (plat, plon, A95), (p['name'], q)
    CHECK.append((p['name'], q['inc'], q['dec'], q['alpha95'], q['k'], plat, plon, A95))
    # mean direction: in the majority polarity of the sites
    maj = 'r' if (sub['dir_inc'] < 0).sum() > len(sub) / 2 else 'n'
    dm = direction_of(sub, maj)
    text = p['text']
    if maj == 'r' and p['polarity'] == 'n':
        text += (f' The mean direction is given in the majority (reversed) polarity; Table 3 gives '
                 f'it in normal polarity as {pdec}/{pinc}.')
    if p.get('tests', {}).get('reversal_test'):
        key = 'Lc' if p['name'].startswith('Pole Lc') else 'LAc'
        text += (f' Recalculated: angle between the R and N means {TESTS[key][0]:.1f} deg, critical '
                 f'angle {TESTS[key][1]:.1f} deg (Watson V simulation, 10000 draws).')
    text += (' Pole: Fisher mean of the polarity-unified site VGPs, reported in the northern '
             'hemisphere as in the paper.')
    row = dict(**GEO, result_name=p['name'], result_type='a', result_quality='g',
               sites=':'.join(p['sites']), lat_s=sub['lat'].min(), lat_n=sub['lat'].max(),
               lon_w=sub['lon'].min(), lon_e=sub['lon'].max(), **p['age'],
               geologic_classes='Intrusive', lithologies=':'.join(sorted(set(sub['lithologies']))),
               method_codes='SO-SM:LP-LT:LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP',
               citations='This study' + ('' if set(p['sites']) <= set(lamp['site'])
                                         else f':{EB}'),
               dir_tilt_correction=0, dir_dec=r1(dm['dec']), dir_inc=r1(dm['inc']),
               dir_alpha95=r1(dm['alpha95']), dir_k=r1(dm['k']), dir_n_sites=int(dm['n']),
               pole_comp_name='characteristic', pole_lat=r1(q['inc']), pole_lon=r1(q['dec']),
               pole_alpha95=r1(q['alpha95']), pole_k=r1(q['k']), pole_n_sites=int(q['n']),
               pole_reversed_perc=round(100 * (sub['dir_inc'] < 0).sum() / len(sub), 1),
               description=text, **p.get('tests', {}))
    if not set(p['sites']) & QUEBEC:
        row['state_province'] = 'Ontario'
    if p.get('gpmdb'):
        row['external_database_ids'] = f'GPMDB[{p["gpmdb"]}]'
    if p.get('pole_vv_q'):
        row['pole_vv_q'] = p['pole_vv_q']
    locations.append(row)
locations = pd.DataFrame(locations)
for c in ('dir_n_sites', 'pole_n_sites', 'pole_vv_q', 'dir_tilt_correction'):
    locations[c] = locations[c].astype('Int64')

# ---- ages table: the dates (Queen et al. 1996, 2 sigma; Krogh et al. 1987, 95%) ----------
PERIOD = dict(location=LOCATION, age_unit='Ma', timescale_eon='Proterozoic',
              timescale_era='Mesoproterozoic', timescale_period='Stenian')
ages = pd.DataFrame([
    dict(**PERIOD, site='D9', age=1143, age_sigma=6, method_codes='GM-UPB', citations=QUEEN,
         description='U-Pb perovskite 1143 +/- 12 Ma (2 sigma; age_sigma is 1 sigma), average '
                     '207Pb/206Pb age, sample LH87-10 from lamprophyre dyke D9 (Wawa area).'),
    dict(**PERIOD, site='D15', age=1145, age_low=1135, age_high=1160, method_codes='GM-UPB',
         citations=QUEEN,
         description='U-Pb perovskite upper-intercept age 1145 +15/-10 Ma (2 sigma), sample LH85-4B '
                     'from lamprophyre dyke D15 (Marathon).'),
    dict(**PERIOD, age=1141, age_low=1132, age_high=1150, method_codes='GM-UPB', citations=QUEEN,
         description='U-Pb perovskite 1141 +/- 9 Ma (2 sigma), sample LH91-8 from the Slate Islands '
                     'lamprophyre dyke (dyke D17, sampled for paleomagnetism but not used).'),
    dict(**PERIOD, age=1144, age_low=1140, age_high=1150, method_codes='GM-UPB', citations=QUEEN,
         description='Composite U-Pb perovskite upper-intercept age 1144 +6/-4 Ma (2 sigma) of the '
                     'Marathon and Wawa lamprophyre dykes.'),
    dict(**PERIOD, age=1144, age_sigma=3.5, method_codes='GM-ARAR', citations=QUEEN,
         description='40Ar/39Ar phlogopite 1144 +/- 7 Ma (2 sigma; age_sigma is 1 sigma), average of '
                     'the plateau segments of sample KL1, Kapuskasing Uplift lamprophyres (not '
                     'paleomagnetic sites).'),
    dict(**PERIOD, site='A6-1-a', age=1140.6, age_sigma=1.0, method_codes='GM-UPB',
         citations=KROGH,
         description='U-Pb baddeleyite 1140.6 +/- 2 Ma (95% confidence; age_sigma is 1 sigma), Great '
                     'Abitibi dyke (Krogh et al., 1987, Geological Association of Canada Special '
                     'Paper 34, 147-152).'),
    dict(**PERIOD, site='A6-2-a', age=1140.6, age_sigma=1.0, method_codes='GM-UPB',
         citations=KROGH,
         description='U-Pb baddeleyite 1140.6 +/- 2 Ma (95% confidence; age_sigma is 1 sigma), Great '
                     'Abitibi dyke (Krogh et al., 1987, Geological Association of Canada Special '
                     'Paper 34, 147-152).'),
    dict(**PERIOD, site='A3-abi', age=1142, age_low=1139, age_high=1145, method_codes='GM-UPB',
         citations='This study',
         description='U-Pb baddeleyite 1142 +/- 3 Ma, Wawa Abitibi dyke A3 (L. M. Heaman, unpublished, '
                     'cited in Sage, 2000, Ontario Geological Survey Open File Report 6106; listed in '
                     'Piispa et al., 2018, Table 1); confidence level not stated.'),
])


# ---- write ----------------------------------------------------------------------------
def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False)


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locations, 'locations', os.path.join(OUT, 'locations.txt'))
write_magic(ages, 'ages', os.path.join(OUT, 'ages.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Piispa2018_Lamprophyre_dykes_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    for i, (kind, df) in enumerate([('locations', locations), ('sites', sites), ('ages', ages)]):
        if i:
            f.write('>>>>>>>>>>\n')
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False)

print('reversal tests (angle, critical):', {k: (round(a, 2), round(c, 2)) for k, (a, c) in
                                            TESTS.items()})
for name, lat, lon, a95, kk, plat, plon, pa95 in CHECK:
    print(f'  {name:48s} {lat:6.3f}N {lon:7.3f}E A95 {a95:5.2f} K {kk:5.1f} | published '
          f'{plat}/{plon} A95 {pa95}')
print(f'wrote {combined}')

from validate_magic_contribution import validate_upload_file  # noqa: E402
if not validate_upload_file(combined, tables=['locations', 'sites', 'ages']):
    raise SystemExit('validation failed')
