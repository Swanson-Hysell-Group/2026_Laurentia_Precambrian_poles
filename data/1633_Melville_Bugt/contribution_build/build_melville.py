"""Build the Melville Bugt dyke swarm MagIC contribution + notebook (1633_Melville).

Source: a student MagIC contribution for Halls, Hamilton & Denyszyn (2011),
"The Melville Bugt Dyke Swarm of Greenland" (Dyke Swarms book ch. 27,
doi:10.1007/978-3-642-12496-9_27), audited against the paper.

Audited fixes (instructor review): pole reported with an internally consistent
k/A95 pair (a true Fisher mean of the 9 dyke VGPs gives k=31.3, A95=9.3, not the
paper-A95 8.7 with a recomputed k); QT dir_n_samples 8 -> 9; the chilled-margin
site "OF Margin" is the margin of the OF dyke and is not counted as an
independent dyke in the 9-dyke pole. The mixed-polarity pole can be reported at
+5 (here) or its -5 antipode.
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Halls2011_MelvilleBugt_magic_source.txt')

data = open(SRC, encoding='latin-1').read()
sites = locs = None
for b in data.split('>>>>>>>>>>'):
    b = b.strip()
    if not b:
        continue
    kind = b.splitlines()[0].split('\t')[1].strip()
    df = pd.read_csv(io.StringIO('\n'.join(b.splitlines()[1:])), sep='\t')
    if kind == 'sites':
        sites = df
    elif kind == 'locations':
        locs = df

# QT dir_n_samples 8 -> 9
sites.loc[sites['site'] == 'QT', 'dir_n_samples'] = 9

# 9-dyke pole (primary, VGP present, exclude the OF chilled margin)
nine = sites[(sites['dir_comp_name'] == 'primary') & (sites['vgp_lat'].notna()) &
             (sites['site'] != 'OF Margin')]
blk = pmag.flip(ipmag.make_di_block(nine['vgp_lon'].tolist(), nine['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)
locs['pole_lat'] = round(p['inc'], 1)
locs['pole_lon'] = round(p['dec'], 1)
locs['pole_alpha95'] = round(p['alpha95'], 1)
locs['pole_k'] = round(p['k'], 1)
locs['pole_n_sites'] = int(p['n'])
locs['age'] = 1633
locs['description'] = (str(locs['description'].iloc[0]) +
    ' Pole is the Fisher mean of the 9 dyke VGPs (k=%.1f, A95=%.1f); the OF chilled '
    'margin is the margin of the OF dyke and is not counted separately. Mixed '
    'polarity; the pole may be reported at the -%.1f antipode.' % (p['k'], p['alpha95'], abs(p['inc'])))

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])}')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/1633_Melville.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Melville Bugt dyke swarm ca. 1633 Ma paleomagnetic pole

## Geologic context

The Melville Bugt dyke swarm trends north-northwest for ~1000 km along the west
coast of Greenland (Halls, Hamilton & Denyszyn, 2011). The dykes are alkaline
trachybasalts, 100-200 m wide, with a remarkably constant geochemistry along the
swarm. Halls et al. proposed the swarm as a possible Laurentia-Greenland
connection to the 1.5-1.6 Ga Fennoscandian rapakivi-granite province. U-Pb
baddeleyite dating gives ca. 1633 Ma (1635 ± 2.7, 1632 ± 1.1, 1629.4 ± 0.8 Ma).

## Pole

This notebook recreates the Melville Bugt pole at the site level from the 9 dykes
that carry a primary magnetite remanence (Halls et al., 2011). The swarm is of
mixed polarity (a positive reversal test), and chilled-margin sampling with
geochemical fingerprinting was used to avoid double-counting duplicate dykes.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data

Site means are loaded from `../data/1633_Melville_Bugt/`. The pole uses the 9
dykes with a stable primary direction; secondary-component and present-field
(PEF) sites, and the OF chilled margin (the margin of the OF dyke), are excluded.
Directions are in-situ (dykes).""")

co("""sites_geo, _ = pt.load_magic_sites('../data/1633_Melville_Bugt/sites.txt')
primary = sites_geo[(sites_geo['dir_comp_name'] == 'primary') &
                    (sites_geo['vgp_lat'].notna()) &
                    (sites_geo['site'] != 'OF Margin')].reset_index(drop=True)
study_lat = round(primary['lat'].mean(), 1)
study_lon = round(primary['lon'].mean(), 1)
print(f'{len(primary)} dykes for the pole; study locality ~{study_lat} N, {study_lon} E')
primary[['site', 'lat', 'lon', 'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k',
         'dir_n_samples', 'vgp_lat', 'vgp_lon']]""")

md("""## Mean pole and the reversal test

The swarm carries both normal and reversed dykes ~180° apart (a positive reversal
test); polarity is unified before the Fisher mean.""")

co("""ipmag.plot_net()
ipmag.plot_di(primary['dir_dec'].tolist(), primary['dir_inc'].tolist(),
              color='blue', marker='o')
plt.title('Melville Bugt dyke directions (in-situ; mixed polarity)')
plt.show()

vgp_block, pole_mean = pt.compute_mean_pole(primary, unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
print('\\nprior compilation (GPMDB 9495): 5 N / 273.8 E, A95 8.7, N=9 (reported at the +5 pole;')
print('Halls et al. 2011 report the -5 S antipode)')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('Melville Bugt dyke swarm pole')
plt.show()""")

co("""dir_block, dir_mean = pt.compute_mean_direction(primary, unify_polarity=True)
ipmag.print_direction_mean(dir_mean)""")

md("""## Paleosecular variation and VGP-shape diagnostics""")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

co("""fishqq_result = pt.fishqq_vgps(primary, unify_polarity=True)
fishqq_result""")

co('''try:
    svei_result = pt.svei_test_vgps(primary, study_lon, study_lat, model='TK03_GAD', plot=True)
except TypeError:
    svei_result = pt.svei_test_vgps(primary, study_lon, study_lat, model='TK03_GAD', plot=False)
    print('(SVEI elongation plot skipped: E below the TK03.GAD model minimum)')
print(f"paleolatitude = {svei_result['lat']:.1f} deg; elongation E = {svei_result['E']:.2f} "
      f"({'consistent' if svei_result['E_result'] else 'inconsistent'} with TK03.GAD)")''')

md("## The Melville Bugt pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Melville Bugt diabase dykes', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=1633)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | U-Pb baddeleyite 1633 Ma (1635 ± 2.7, 1632 ± 1.1, 1629.4 ± 0.8 Ma; Halls et al., 2011). |
| 2 | Techniques and statistical analysis | **1** | AF + thermal demagnetization with vector analysis; N = 9 dyke means, K = 31.3, A95 = 9.3° passes the Deenen et al. (2011) envelope. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence (thermomagnetic analysis; Denyszyn, 2008). |
| 4 | Field tests constrain age of magnetization | **0** | No baked-contact or fold test (chilled-margin sampling was used for geochemical dyke fingerprinting, not as a country-rock baked-contact test). |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous NW Greenland (Laurentia-Greenland) basement; vertical dykes need no tilt correction. |
| 6 | Presence of reversals | **1** | The swarm is of mixed polarity (normal and reversed dykes ~180° apart). |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles. |
| | Total | **6/7** | Grade B |""")

md("""## Nordic workshop summary

The Melville Bugt pole is recreated at the site level from the 9 primary-remanence
dykes of Halls et al. (2011). The recreation reproduces the prior compilation
pole (5°N/273.8°E, A95 8.7, GPMDB 9495) with an internally consistent k/A95 pair
(k 31.3, A95 9.3). Mixed polarity scores R6 = 1; the swarm is well dated by U-Pb
(R1 = 1). The pole is reported at the +5 position; Halls et al. give the −5
antipode.""")

co("""melville_summary = pt.make_nordic_summary(
    terrane='Laurentia-Greenland',
    rockname='Melville Bugt diabase dykes',
    sites=primary,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Primary characteristic remanent magnetization carried by magnetite; mixed polarity',
    tests='R+ (positive reversal test: normal and reversed dykes ~180 deg apart)',
    gpmdb_number='9495',
    percent_reversed=44,
    demag_code=4,
    R1=1, R2=1, R3=1, R4='', R5=1, R6=1, R7=1, Grade='B',
    nominal_age=1633, lomagage=1628, himagage=1638,
    REF_method='U-Pb baddeleyite ages on the Melville Bugt dyke swarm: 1635 +/- 2.7, 1632 +/- 1.1, and 1629.4 +/- 0.8 Ma (Halls, Hamilton & Denyszyn, 2011).',
    POLE_AUTHORS='Halls, H. C., Hamilton, M. A., & Denyszyn, S. W.',
    YEAR=2011,
    JOURNAL='Dyke Swarms: Keys for Geodynamic Interpretation (Springer)',
    VOLUME='',
    VPAGES='509-535',
    TITLE='The Melville Bugt Dyke Swarm of Greenland: A Connection to the 1.5-1.6 Ga Fennoscandian Rapakivi Granite Province?',
    COMMENT='Melville Bugt pole recreated at the site level from the 9 primary-remanence dykes of Halls et al. (2011) (audited contribution; QT dir_n_samples 8->9; OF chilled margin not counted as a separate dyke). Recreated 4.0N/274.5E with an internally consistent k=31.3/A95=9.3 (the source mixed the paper A95 8.7 with a recomputed k); reproduces the prior compilation pole (5N/273.8E, GPMDB 9495). Mixed polarity -> R6=1; U-Pb ca. 1633 Ma -> R1=1. Reported at +5; Halls et al. give the -5 antipode. No field test (R4=0). Grade B.'
)
pt.save_nordic_summary(melville_summary, '1633_Melville_Bugt')
melville_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
