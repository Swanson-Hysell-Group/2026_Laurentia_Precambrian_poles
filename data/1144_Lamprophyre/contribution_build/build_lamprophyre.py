"""Build the NW Ontario lamprophyre dyke MagIC contribution + notebook (1144_Lamprophyre).

Source: a student MagIC contribution (id 20635) for Piispa, Smirnov, Pesonen &
Mitchell (2018), "Paleomagnetism and Geochemistry of ~1144-Ma Lamprophyre Dikes,
Northwestern Ontario" (JGR Solid Earth), audited against the paper.

Audited fixes (instructor review): the location pole_alpha95 was the directional
value (5.5); replaced by recomputing the VGP-Fisher mean (the published Pole Lc =
58.0 N / 223.3 E, A95 9.2, K 17.2, N=19); all site `location` values set to the
declared "Southern Superior Province"; Queen et al. (1996) age reference DOI
10.1139/e96-072 added.

The repo treats the lamprophyre dykes (this pole) and the Abitibi dykes
(1141_Abitibi) separately; the prior compilation combined them (GPMDB 9887:
55.8/220, N=27). Lamprophyre-only pole: 19 dykes, mixed polarity (a reversal
test), positive baked-contact test (dykes D1, MM2). Age 1144 +/- 7 Ma
(40Ar/39Ar phlogopite and U-Pb perovskite; Queen et al., 1996).
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Piispa2018_Lamprophyre_magic_20635_source.txt')

data = open(SRC, encoding='latin-1').read()
sites = None
for b in data.split('>>>>>>>>>>'):
    if b.strip().startswith('tab') and 'sites' in b.strip().splitlines()[0]:
        sites = pd.read_csv(io.StringIO('\n'.join(b.strip().splitlines()[1:])), sep='\t')

sites['location'] = 'Southern Superior Province'  # unify orphan "Northwestern Ontario"

ABITIBI = ['A4a', 'A7abi', 'A3-abi', 'A5-abi', 'A6-1-a', 'A6-2-a', 'A8-abi', 'SL13']
BAKED = ['MM2 Baked', 'MM2 Unbaked', 'D1 Baked', 'D1 Unbaked']
lamp = sites[(~sites['site'].isin(ABITIBI + BAKED)) & (sites['vgp_lat'].notna())]
blk = pmag.flip(ipmag.make_di_block(lamp['vgp_lon'].tolist(), lamp['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)

locs = pd.DataFrame([{
    'location': 'Southern Superior Province', 'location_type': 'Region',
    'result_name': 'NW Ontario lamprophyre dykes ca. 1144 Ma pole', 'result_type': 'a',
    'sites': ':'.join(lamp['site'].tolist()),
    'method_codes': 'LP-DIR-T:LP-DIR-AF:DE-BFL:DE-FM:DE-VGP:ST-C:GM-AR-AR:GM-UPB-CC',
    'citations': '10.1029/2018JB015992:10.1139/e96-072', 'geologic_classes': 'Intrusive',
    'lithologies': 'Lamprophyre',
    'lat_s': sites['lat'].min(), 'lat_n': sites['lat'].max(), 'lon_w': sites['lon'].min(), 'lon_e': sites['lon'].max(),
    'age': 1144, 'age_low': 1139, 'age_high': 1150, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1), 'pole_alpha95': round(p['alpha95'], 1),
    'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'description': 'NW Ontario lamprophyre dyke pole (Piispa et al., 2018 Pole Lc): VGP-Fisher mean of 19 dykes; mixed polarity (reversal test); positive baked-contact test (D1, MM2). Age 1144+/-7 Ma (40Ar/39Ar phlogopite + U-Pb perovskite; Queen et al., 1996).'}])

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])}')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/1144_Lamprophyre.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# NW Ontario lamprophyre dykes ca. 1144 Ma paleomagnetic pole

## Geologic context

The lamprophyre (and related alkaline) dykes of northwestern Ontario intrude the
southern Superior Province and are interpreted as early magmatism of the ca. 1.1
Ga Midcontinent Rift initiation, coeval with the Abitibi dyke swarm (Piispa et
al., 2018). A combined 40Ar/39Ar phlogopite and U-Pb perovskite age of
1144 ± 7 Ma (Queen et al., 1996) dates the dykes. The dykes carry a steep,
mixed-polarity magnetization shown to be primary by a positive baked-contact
test.

## Pole

This notebook recreates the lamprophyre dyke pole at the dyke level from the 19
dykes of Piispa et al. (2018), excluding the Abitibi dykes (treated separately in
`1141_Abitibi`) and the baked/unbaked contact-test sites.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data and the 19-dyke pole

Loaded from `../data/1144_Lamprophyre/`. The pole uses the 19 lamprophyre dykes;
most are reversed-polarity (steep down-to-southeast), several normal — a reversal
test.""")

co("""sites = pd.read_csv('../data/1144_Lamprophyre/sites.txt', sep='\\t', skiprows=1)
abitibi = ['A4a', 'A7abi', 'A3-abi', 'A5-abi', 'A6-1-a', 'A6-2-a', 'A8-abi', 'SL13']
baked = ['MM2 Baked', 'MM2 Unbaked', 'D1 Baked', 'D1 Unbaked']
lamp = sites[(~sites['site'].isin(abitibi + baked)) & (sites['vgp_lat'].notna())].reset_index(drop=True)
study_lat, study_lon = 48.4, 273.6
print(f'{len(lamp)} lamprophyre dykes; reversed:', lamp[lamp['dir_inc'] < 0]['site'].tolist())
lamp[['site', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95', 'dir_n_samples', 'vgp_lat', 'vgp_lon']]""")

co("""ipmag.plot_net()
ipmag.plot_di(lamp['dir_dec'].tolist(), lamp['dir_inc'].tolist(), color='royalblue', marker='o')
plt.title('Lamprophyre dyke directions (in-situ; mixed polarity)')
plt.show()

vgp_block, pole_mean = pt.compute_mean_pole(lamp, unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
print('\\nPiispa et al. (2018) Pole Lc: 58.0 N / 223.3 E, A95 9.2, N=19')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('NW Ontario lamprophyre dyke pole')
plt.show()""")

co("""dir_block, dir_mean = pt.compute_mean_direction(lamp, unify_polarity=True)
ipmag.print_direction_mean(dir_mean)""")

md("""## Reversal test (McFadden & McElhinny, 1990)

The 19 lamprophyre dykes are of mixed polarity (14 reversed, 5 normal). The
McFadden & McElhinny (1990) reversal test inverts the reversed group to normal
polarity and applies a Watson V common-mean test, returning a classification
(A/B/C) or a negative result.""")

co('''print('normal-polarity dykes  :', lamp[lamp['dir_inc'] > 0]['site'].tolist())
print('reversed-polarity dykes:', lamp[lamp['dir_inc'] < 0]['site'].tolist())

# McFadden & McElhinny (1990) reversal test: Watson V common-mean test + classification
ipmag.reversal_test_MM1990(dec=lamp['dir_dec'].tolist(), inc=lamp['dir_inc'].tolist(),
                           plot_stereo=True)
plt.show()''')

md("""## Field test: baked contact

Dykes D1 and MM2 expose baked host-rock margins that carry the dyke direction,
while the unbaked host away from the contact differs — a positive baked-contact
test.""")

co("""bc = sites[sites['site'].isin(['D1 Baked', 'D1 Unbaked', 'MM2 Baked', 'MM2 Unbaked'])]
print(bc[['site', 'dir_dec', 'dir_inc', 'dir_alpha95']].to_string(index=False))""")

md("## Paleosecular variation")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

md("## The lamprophyre pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('NW Ontario Lamprophyre Dykes', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=1144)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | 40Ar/39Ar phlogopite and U-Pb perovskite 1144 ± 7 Ma (Queen et al., 1996). |
| 2 | Techniques and statistical analysis | **1** | AF + thermal demagnetization, PCA; N = 19 dykes, A95 = 9.2° passes the Deenen et al. (2011) envelope. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence characterized by rock-magnetic experiments (Piispa et al., 2018). |
| 4 | Field tests constrain age of magnetization | **C** | Positive baked-contact test (dykes D1, MM2: baked margins carry the dyke direction, unbaked host differs). |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous southern Superior Province; vertical dykes need no tilt correction. |
| 6 | Presence of reversals | **1** | Mixed polarity (14 reversed, 5 normal); positive McFadden & McElhinny (1990) reversal test (class C). |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles. |
| | Total | **7/7** | Grade A |""")

md("""## Nordic workshop summary

The lamprophyre dyke pole is recreated from the 19 dykes of Piispa et al. (2018),
reproducing the published Pole Lc (58.0°N/223.3°E, A95 9.2°). The prior
compilation combined these dykes with the Abitibi swarm (GPMDB 9887); they are
reported separately here. Positive baked-contact test (R4 = C), mixed polarity
with a positive reversal test (R6 = 1), and the 1144 Ma age give a Grade A,
R = 7 pole.""")

co("""lamp_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='NW Ontario Lamprophyre Dykes',
    sites=lamp,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Primary characteristic remanent magnetization (magnetite); mixed polarity',
    tests='C+ (positive baked-contact test, dykes D1 and MM2) and R+ (positive reversal test, class C; 14 reversed, 5 normal)',
    gpmdb_number='9887',
    magic_id='20635',
    percent_reversed=74,
    demag_code=4,
    R1=1, R2=1, R3=1, R4='C', R5=1, R6=1, R7=1, Grade='A',
    nominal_age=1144, lomagage=1139, himagage=1150,
    REF_method='40Ar/39Ar phlogopite and U-Pb perovskite age 1144 +/- 7 Ma on the NW Ontario lamprophyre dykes (Queen et al., 1996); correlated with the ca. 1141 Ma Abitibi swarm (Piispa et al., 2018).',
    POLE_AUTHORS='Piispa, E. L., Smirnov, A. V., Pesonen, L. J., & Mitchell, R. H.',
    YEAR=2018,
    JOURNAL='Journal of Geophysical Research: Solid Earth',
    VOLUME='123',
    VPAGES='',
    TITLE='Paleomagnetism and Geochemistry of ~1144-Ma Lamprophyre Dikes, Northwestern Ontario: Implications for the Midcontinent Rift',
    COMMENT='NW Ontario lamprophyre dyke pole recreated from the 19 dykes of Piispa et al. (2018) (audited MagIC contribution 20635: pole_alpha95 corrected from the directional 5.5 to the VGP-Fisher A95 9.2; site locations unified to Southern Superior Province; Queen et al. 1996 DOI e96-072 added). Recreated 57.9N/223.3E, A95 9.2, N=19 -- reproduces the published Pole Lc (58.0/223.3). Reported separately from the Abitibi swarm (1141_Abitibi); the prior compilation combined them (GPMDB 9887). Positive baked-contact test (R4=C), mixed polarity (R6=1), 1144+/-7 Ma (R1=1). R=7, Grade A.'
)
pt.save_nordic_summary(lamp_summary, '1144_Lamprophyre')
lamp_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
