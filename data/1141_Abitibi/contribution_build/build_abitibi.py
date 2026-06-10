"""Build the Abitibi dyke swarm MagIC contribution + notebook (1141_Abitibi).

Source: a student MagIC contribution (id 20639) for Ernst & Buchan (1993),
"Paleomagnetism of the Abitibi dyke swarm, southern Superior Province"
(doi:10.1139/e93-150), audited against the paper. The submission carried several
errors (W- vs E-longitude in the VGP fields, dp/dm swaps, an unreal `sites` list,
A1 age contradiction); the pole and per-dyke VGPs are recomputed cleanly here
with pmag.dia_vgp.

Pole: the dataset is updated to the eight-dyke group of Halls et al. (2005), who
reassigned the original A1 dyke to the ca. 2167 Ma Biscotasing swarm and added a
new normal-polarity Abitibi dyke (SL13). The eight-dyke mean (A3-A8 + SL13, A1
excluded) is 50.4 N / 213.4 E, A95 12.6, N=8 -- reproducing the Halls et al.
(2005) / Piispa et al. (2018) modified pole (50.5 N / 213.8 E, A95 12.5; D=297.4,
I=65.5), offset ~8.5 deg from the original Ernst & Buchan (1993) pole. The swarm
is of mixed polarity (A4, A7 reversed) and has a positive baked-contact test
(baked host margins carry the dyke direction; unbaked host differs). Age 1140.6
+/- 2 Ma (U-Pb baddeleyite; Krogh, 1987). SL13 direction/VGP from Piispa et al.
(2018) Table 3 (Halls et al., 2005).
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Ernst1993_Abitibi_magic_20639_source.txt')

data = open(SRC, encoding='latin-1').read()
sites = None
for b in data.split('>>>>>>>>>>'):
    if b.strip().startswith('tab') and 'sites' in b.strip().splitlines()[0]:
        sites = pd.read_csv(io.StringIO('\n'.join(b.strip().splitlines()[1:])), sep='\t')

# recompute signed VGPs for every row from its direction (fixes the W/E-long + dp/dm errors)
def vgp_row(r):
    if pd.isna(r['dir_dec']) or pd.isna(r['dir_alpha95']):
        a = 5.0
    else:
        a = r['dir_alpha95']
    if pd.isna(r['dir_dec']):
        return pd.Series({'vgp_lon': float('nan'), 'vgp_lat': float('nan'),
                          'vgp_dp': float('nan'), 'vgp_dm': float('nan')})
    plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], a, r['lat'], r['lon'])
    return pd.Series({'vgp_lon': round(plon, 1), 'vgp_lat': round(plat, 1),
                      'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1)})

sites = pd.concat([sites.drop(columns=[c for c in ['vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm'] if c in sites.columns]),
                   sites.apply(vgp_row, axis=1)], axis=1)
sites['location'] = 'southern Superior Province'

# Add the new normal-polarity Abitibi dyke reported by Halls et al. (2005), site
# SL13 (direction + VGP from Piispa et al., 2018, Table 3, footnote c). This brings
# the Ernst & Buchan (1993) dataset to the eight-dyke group of Halls et al. (2005),
# who also reassigned the original A1 dyke to the ca. 2167 Ma Biscotasing swarm.
_sl_lat, _sl_lon = 48.35, 275.33
_sp, _slat, _sdp, _sdm = pmag.dia_vgp(315.9, 62.7, 5.1, _sl_lat, _sl_lon)
sl13 = {c: '' for c in sites.columns}
sl13.update({'site': 'SL13', 'location': 'southern Superior Province', 'result_type': 'a',
             'result_quality': 'g', 'method_codes': 'LP-DIR-AF:DE-BFL:DE-FM:DE-VGP',
             'citations': 'Halls et al. 2005 (OGS OFR 6171)', 'geologic_classes': 'Igneous',
             'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
             'lat': _sl_lat, 'lon': _sl_lon, 'age': 1140.6, 'age_sigma': 1.0, 'age_unit': 'Ma',
             'dir_tilt_correction': 0, 'dir_dec': 315.9, 'dir_inc': 62.7, 'dir_alpha95': 5.1,
             'dir_k': 61.0, 'dir_n_samples': 1, 'dir_polarity': 'n',
             'description': 'Dyke SL13: new normal-polarity Abitibi dyke (Halls et al. 2005; Piispa et al. 2018 Table 3). With A1 reassigned to the ca. 2167 Ma Biscotasing swarm, this completes the eight-dyke modified group mean.',
             'vgp_lon': round(_sp, 1), 'vgp_lat': round(_slat, 1),
             'vgp_dp': round(_sdp, 1), 'vgp_dm': round(_sdm, 1)})
sites = pd.concat([sites, pd.DataFrame([sl13])], ignore_index=True)

# 8-dyke pole (dyke means, exclude A1; SL13 added from Halls et al. 2005)
dykes = sites[(sites['result_type'] == 'a') & (sites['site'] != 'A1')]
blk = pmag.flip(ipmag.make_di_block(dykes['vgp_lon'].tolist(), dykes['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)

locs = pd.DataFrame([{
    'location': 'southern Superior Province', 'location_type': 'Region',
    'result_name': 'Abitibi dyke swarm ca. 1141 Ma pole', 'result_type': 'a',
    'sites': ':'.join(dykes['site'].tolist()),
    'method_codes': 'LP-DIR-T:LP-DIR-AF:DE-BFL:DE-FM:DE-VGP:ST-C:GM-UPB-CC',
    'citations': '10.1139/e93-150', 'geologic_classes': 'Intrusive', 'lithologies': 'Diabase',
    'lat_s': sites['lat'].min(), 'lat_n': sites['lat'].max(), 'lon_w': sites['lon'].min(), 'lon_e': sites['lon'].max(),
    'age': 1141, 'age_low': 1139, 'age_high': 1143, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1), 'pole_alpha95': round(p['alpha95'], 1),
    'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'description': 'Abitibi dyke swarm pole, Fisher mean of 8 dyke-mean VGPs (A3-A8 + the new Halls et al. 2005 dyke SL13; the original A1 dyke reassigned to the ca. 2167 Ma Biscotasing swarm). Reproduces the Halls et al. (2005) / Piispa et al. (2018) modified pole (50.5N/213.8E, A95 12.5). Mixed polarity (A4, A7 reversed); positive baked-contact test. U-Pb baddeleyite 1140.6+/-2 Ma (Krogh, 1987).'}])

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])}')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/1141_Abitibi.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Abitibi dyke swarm ca. 1141 Ma paleomagnetic pole

## Geologic context

The Abitibi dykes are northeast-trending diabase dykes of the southern Superior
Province (Ontario/Quebec) emplaced ca. 1141 Ma that predates the initiation of
the Midcontinent Rift ca. 1109 Ma (Ernst & Buchan, 1993). A U-Pb baddeleyite age
of 1140.6 ± 2 Ma (Krogh, 1987) dates the swarm. The dykes carry a steep,
dominantly down-to-the-northwest magnetization of mixed polarity. It being a
primary thermal remanent magnetization is supported by a positive baked-contact
test.

## Pole

This notebook recreates the Abitibi pole at the dyke level from the per-dyke mean
directions of Ernst & Buchan (1993), recomputing signed VGPs. The dataset is
updated to the eight-dyke group of Halls et al. (2005): the original A1 dyke is
reassigned to the ca. 2167 Ma Biscotasing swarm and excluded, and a new
normal-polarity dyke (SL13) is added.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data and the eight-dyke pole

The contribution holds individual site results and per-dyke means
(`result_type == 'a'`). The pole uses the eight dyke means — A3-A8 plus the new
Halls et al. (2005) dyke SL13 — excluding A1 (reassigned to the Biscotasing
swarm). Dykes A4 and A7 carry the reversed polarity (a reversal test).""")

co("""sites = pd.read_csv('../data/1141_Abitibi/sites.txt', sep='\\t', skiprows=1)
dykes = sites[(sites['result_type'] == 'a') & (sites['site'] != 'A1')].reset_index(drop=True)
study_lat, study_lon = 48.0, 279.0
print(f'{len(dykes)} dyke means for the pole:', dykes['site'].tolist())
print('reversed-polarity dykes:', dykes[dykes['dir_inc'] < 0]['site'].tolist())
dykes[['site', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95', 'dir_n_samples', 'vgp_lat', 'vgp_lon']]""")

co("""ipmag.plot_net()
ipmag.plot_di(dykes['dir_dec'].tolist(), dykes['dir_inc'].tolist(), color='royalblue', marker='o')
plt.title('Abitibi dyke mean directions (in-situ; mixed polarity, 8 dykes)')
plt.show()

vgp_block, pole_mean = pt.compute_mean_pole(dykes, unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
print('\\nHalls et al. (2005) / Piispa et al. (2018) modified pole: 50.5 N / 213.8 E, A95 12.5, N=8')
print('original Ernst & Buchan (1993) pole (GPMDB 7193): 48.8 N / 215.5 E, A95 14.1, N=7')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('Abitibi dyke swarm pole')
plt.show()""")

md("""## The Halls et al. (2005) revision

Ernst & Buchan (1993) based their group-mean pole on five normally and two
reversely magnetized dykes, with two independent normal directions obtained from
distinct lithological units of the Great Abitibi Dike (site A6). The primary
nature of the normal-polarity directions was established by a baked-contact test
against a Matachewan dyke, and the presence of opposite-polarity directions in
neighbouring dykes argued against pervasive regional remagnetization. The
reversed-polarity directions are steeper than their normal counterparts, so the
primary origin of the reversed remanence remains somewhat tentative (Ernst &
Buchan, 1993; Piispa et al., 2018).

Subsequently, one of the normal-polarity dykes (A1) was shown to belong to the
ca. 2167 Ma Biscotasing swarm, and a new normal-polarity Abitibi dyke direction
was reported (Halls et al., 2005). With A1 removed and the new dyke (SL13, here
taken from Piispa et al., 2018, Table 3) added, the modified group-mean direction
is D = 297.4°, I = 65.5° (α95 = 8.3°, N = 8), corresponding to a pole at
50.5°N, 213.8°E (A95 = 12.5°) — offset ~8.5° from the original Ernst & Buchan
(1993) pole (Halls et al., 2005; Piispa et al., 2018). This notebook adopts that
revised eight-dyke dataset, and the recreated pole above reproduces it.""")

co("""dir_block, dir_mean = pt.compute_mean_direction(dykes, unify_polarity=True)
ipmag.print_direction_mean(dir_mean)""")

md("""## Reversal test (McFadden & McElhinny, 1990)

The swarm carries six normal-polarity dykes (A3, A5, A6_1, A6_2, A8, SL13) and
two reversed-polarity dykes (A4, A7). The McFadden & McElhinny (1990) reversal
test inverts the reversed group to normal polarity and applies a Watson V
common-mean test, returning a classification (A/B/C) or a negative/indeterminate
result. The two reversed dykes are steeper than the normal group, so the means
remain distinguishable and the test is negative/indeterminate at this resolution;
the presence of both polarities nonetheless documents that the swarm records
reversals (Meert R6).""")

co('''print('normal-polarity dykes  :', dykes[dykes['dir_inc'] > 0]['site'].tolist())
print('reversed-polarity dykes:', dykes[dykes['dir_inc'] < 0]['site'].tolist())

# McFadden & McElhinny (1990) reversal test: Watson V common-mean test + classification
ipmag.reversal_test_MM1990(dec=dykes['dir_dec'].tolist(), inc=dykes['dir_inc'].tolist(),
                           plot_stereo=True)
plt.show()''')

md("""## Field test: baked contact

The swarm includes baked host-rock margins (sites A6_2-16B, A6_2-21B) that carry
the dyke direction, while the unbaked host away from the contact (A6_2-16U,
A6_2-21U) carries a distinct shallow direction — a positive baked-contact test.""")

co("""bc = sites[sites['site'].isin(['A6_2-16B', 'A6_2-16U', 'A6_2-21B', 'A6_2-21U'])]
print(bc[['site', 'dir_dec', 'dir_inc', 'dir_alpha95', 'description']].to_string(index=False) if len(bc) else 'baked-contact rows not present')""")

md("""## Paleosecular variation""")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

md("## The Abitibi pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Abitibi Dykes', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=1141)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | U-Pb baddeleyite 1140.6 ± 2 Ma (Krogh, 1987). |
| 2 | Techniques and statistical analysis | **0** | AF + thermal demagnetization with PCA, but only B = 8 dykes and A95 = 12.6° is large; Piispa et al. (2018) note the pole is derived from a limited number of independently cooled units and does not represent the time-averaged field, so paleosecular variation is marginally sampled. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence characterized in the source study (Ernst & Buchan, 1993). |
| 4 | Field tests constrain age of magnetization | **C** | Positive baked-contact test (baked host margins carry the dyke direction; unbaked host differs). |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous southern Superior Province; vertical dykes need no tilt correction. |
| 6 | Presence of reversals | **1** | Mixed polarity (dykes A4 and A7 reversed); the reversal test is negative/indeterminate as the reversed directions are steeper than the normal group. |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles. |
| | Total | **6/7** | Grade A |""")

md("""## Nordic workshop summary

The Abitibi pole is recreated at the dyke level from the eight dyke means
(A3-A8 plus the new Halls et al. 2005 dyke SL13; A1 reassigned to the Biscotasing
swarm) of Ernst & Buchan (1993), with signed VGPs recomputed (the source carried
W/E-longitude and dp/dm errors). The recreation (50.4°N/213.4°E, A95 12.6°, N=8)
reproduces the Halls et al. (2005) / Piispa et al. (2018) modified pole
(50.5°N/213.8°E, A95 12.5°), which is offset ~8.5° from the original Ernst &
Buchan (1993) pole (GPMDB 7193). Positive baked-contact test (R4 = C) and mixed
polarity (R6 = 1); the limited number of cooling units keeps R2 = 0.""")

co("""abitibi_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Abitibi Dykes',
    sites=dykes,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Primary characteristic remanent magnetization (magnetite); mixed polarity; eight-dyke modified mean of Halls et al. (2005)',
    tests='C+ (positive baked-contact test) and R (mixed polarity; dykes A4, A7 reversed, reversal test negative/indeterminate)',
    gpmdb_number='7193',
    percent_reversed=25,
    demag_code=4,
    R1=1, R2=0, R3=1, R4='C', R5=1, R6=1, R7=1, Grade='A',
    nominal_age=1141, lomagage=1139, himagage=1143,
    REF_method='U-Pb baddeleyite age 1140.6 +/- 2 Ma on the Abitibi dyke swarm (Krogh, 1987; compiled by Ernst & Buchan, 1993).',
    POLE_AUTHORS='Ernst, R. E., & Buchan, K. L.',
    YEAR=1993,
    JOURNAL='Canadian Journal of Earth Sciences',
    VOLUME='30',
    VPAGES='1886-1897',
    TITLE='Paleomagnetism of the Abitibi dyke swarm, southern Superior Province, and implications for the Logan Loop',
    COMMENT='Abitibi dyke swarm pole recreated at the dyke level from the eight dyke means (A3-A8 + SL13) of the Halls et al. (2005) modified dataset: the original A1 dyke is reassigned to the ca. 2167 Ma Biscotasing swarm and excluded, and the new normal-polarity dyke SL13 (direction/VGP from Piispa et al. 2018 Table 3, after Halls et al. 2005) is added. Signed per-dyke VGPs recomputed with pmag.dia_vgp (source contribution had W/E-longitude and dp/dm errors). Recreated 50.4N/213.4E, A95 12.6, N=8 -- reproduces the Halls et al. (2005) / Piispa et al. (2018) modified pole 50.5N/213.8E (A95 12.5; D=297.4/I=65.5), offset ~8.5 deg from the original Ernst & Buchan (1993) pole (GPMDB 7193, 48.8N/215.5E). Positive baked-contact test (R4=C); mixed polarity A4/A7 (R6=1), reversal test negative/indeterminate (reversed dykes steeper than normal). U-Pb baddeleyite 1140.6+/-2 Ma (R1=1). Limited number of cooling units, Piispa et al. (2018) note the pole does not represent the time-averaged field -> R2=0. R=6, Grade A.'
)
pt.save_nordic_summary(abitibi_summary, '1141_Abitibi')
abitibi_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
