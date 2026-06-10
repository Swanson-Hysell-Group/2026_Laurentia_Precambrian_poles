"""Build the Uinta Mountain Group MagIC contribution + notebook (754_Uinta).

Source: a student MagIC contribution (id 20680) for Weil, Geissman & Ashby
(2006), Precambrian Research 147, 234-259 (doi:10.1016/j.precamres.2006.01.017),
audited against the paper.

Audited fixes (instructor review): site result_type a -> i; the eastern
overprint localities (Irish Canyon, Cross Mountain, Lone Mountain, Juniper
Mountain), which carry only a north-directed recent VRM (poles at 80-88 N), are
relabeled as "UMG present-field overprint" rather than the headline UMG pole;
primary-locality result_name updated to the ca. 754 Ma UMG pole.

Weil et al. (2006) UMG ChRM pole: 0.8 N, 161.3 E, a95 4.6, N=9 sampling
localities (79 sites). Hematite-cemented sandstone/quartzite; dual-polarity ChRM
(west = normal, east = reverse), shallow. Age: Tonian (~717-766 Ma; U-Pb ash
742 +/- 6 Ma, Karlstrom et al. 2000; entire UMG <766 Ma, pre-717 Ma, Dehler).
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Weil2006_UMG_magic_20680_source.txt')

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

OVERPRINT = ['Irish Canyon Section', 'Cross Mountain', 'Lone Mountain', 'Juniper Mountain']
# relabel overprint vs primary
locs.loc[locs['location'].isin(OVERPRINT), 'result_name'] = 'UMG present-field overprint'
prim_mask = (~locs['location'].isin(OVERPRINT)) & (~locs['location'].str.contains('Uinta Mountain Group ca', na=False))
locs.loc[prim_mask | locs['location'].str.contains('Uinta Mountain Group ca', na=False),
         'result_name'] = 'Uinta Mountain Group ca. 754 Ma pole'
# site result_type a -> i
if 'result_type' in sites.columns:
    sites['result_type'] = 'i'

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))

# the 9 primary locality poles (each one VGP for the grand mean)
prim = locs[prim_mask].copy()
blk = ipmag.make_di_block(prim['pole_lon'].tolist(), prim['pole_lat'].tolist())
p = pmag.fisher_mean(blk)
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; UMG pole '
      f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])} ({len(prim)} localities)')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/754_Uinta.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Uinta Mountain Group ca. 754 Ma paleomagnetic pole

## Geologic context

The Uinta Mountain Group (UMG) is one of the thickest and best-preserved
Neoproterozoic sedimentary successions in North America — 4-7 km of quartzose to
arkosic sandstone and interbedded mudstone, much of it hematite-cemented, exposed
in the east-west Laramide Uinta uplift of northeast Utah, northwest Colorado, and
southwest Wyoming (Weil, Geissman & Ashby, 2006). The UMG is interpreted as a
Neoproterozoic aulacogen/rift fill lacking burial metamorphism, deposited in the
Tonian (~717-766 Ma; U-Pb detrital-zircon ash 742 ± 6 Ma, Karlstrom et al., 2000;
correlation with the Chuar Group places it pre-717 Ma).

## Pole

This notebook recreates the UMG characteristic-remanence pole at the site level
from the 9 sampling localities (79 sites) of Weil et al. (2006). The primary ChRM
is a west- or east-directed (dual-polarity), shallow, hematite-carried
magnetization; the four eastern localities that yielded only a north-directed
present-field viscous overprint are shown separately.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Locality poles and the UMG mean

The locations table holds a VGP (locality pole) for each sampling locality. The
UMG pole is the Fisher mean of the 9 primary-ChRM locality VGPs; the four
present-field overprint localities (Irish Canyon, Cross Mountain, Lone Mountain,
Juniper Mountain) are excluded.""")

co("""locs = pd.read_csv('../data/754_Uinta_Mountain_Group/locations.txt', sep='\\t', skiprows=1)
overprint = ['Irish Canyon Section', 'Cross Mountain', 'Lone Mountain', 'Juniper Mountain']
primary = locs[(~locs['location'].isin(overprint)) &
               (~locs['location'].str.contains('Uinta Mountain Group ca', na=False))].copy()
study_lat, study_lon = 40.8, 250.7

# treat each locality VGP as one site for the grand-mean pole
loc_sites = pd.DataFrame({'site': primary['location'], 'vgp_lat': primary['pole_lat'],
                          'vgp_lon': primary['pole_lon'], 'dir_tilt_correction': 100,
                          'dir_n_samples': primary['pole_n_sites'].fillna(1).astype(int)})
print(f'{len(loc_sites)} primary localities')
vgp_block, pole_mean = pt.compute_mean_pole(loc_sites, unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
print('\\nWeil et al. (2006) UMG pole: 0.8 N / 161.3 E, A95 4.6, N=9')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('Uinta Mountain Group pole (mean of 9 locality VGPs)')
plt.show()""")

md("""## Present-field overprint localities (context)

The four eastern localities carry a north-directed, steep, present-field viscous
remanence — their VGPs plot near the present pole, not with the UMG ChRM.""")

co("""op = locs[locs['location'].isin(overprint)]
print(op[['location', 'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_n_sites']].to_string(index=False))""")

md("""## Paleosecular variation

The dual-polarity ChRM (west-normal / east-reverse) and the spread of locality
VGPs indicate the magnetization averaged secular variation; the Deenen et al.
(2011) test is applied to the 9-locality mean.""")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

md("## The UMG pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Uinta Mountain Group', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=754)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **0** | Only a broad Tonian bracket (~717-766 Ma; U-Pb detrital-zircon ash 742 ± 6 Ma is a maximum depositional age; Karlstrom et al., 2000; Dehler). |
| 2 | Techniques and statistical analysis | **1** | Thermal + chemical demagnetization, PCA, IRM mineralogy; N = 9 localities (79 sites); the dual-polarity ChRM averages secular variation. |
| 3 | Magnetic mineralogy characterized | **1** | Remanence carried by hematite, characterized by high unblocking temperatures (660-680 °C) and three-axis IRM (Weil et al., 2006). |
| 4 | Field tests constrain age of magnetization | **0** | Fold tests were local/inconclusive; no decisive regional field test. |
| 5 | Structural control / tectonic coherence | **0** | Detrital hematite sedimentary remanence with potential inclination shallowing (uncorrected); within a Laramide uplift. |
| 6 | Presence of reversals | **0** | Dual polarity is present (west-normal / east-reverse ChRM) but was not formalized as a passing reversal test in the source scoring. |
| 7 | No resemblance to younger poles | **0** | The shallow, near-equatorial pole overlaps other mid-Neoproterozoic Laurentia poles. |
| | Total | **2/7** | Grade B |""")

md("""## Nordic workshop summary

The UMG pole is recreated as the Fisher mean of the 9 primary-ChRM locality VGPs
of Weil et al. (2006), reproducing the published pole (0.8°N/161.3°E, A95 4.6°,
prior compilation GPMDB 9290). The four eastern present-field overprint localities
are excluded (relabeled in the contribution). The R-scores follow the prior
compilation (R = 2, Grade B).

**Flag for review:** the paper documents a dual-polarity ChRM that could support
R6 = 1, and recent geochronology (Dehler; Karlstrom 742 ± 6 Ma) suggests refining
the age toward ~720-760 Ma; both are noted but the prior compilation's
conservative scoring is retained here.""")

co("""umg_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Uinta Mountain Group',
    sites=loc_sites,
    dir_mean={'dec': 270.2, 'inc': 2.0, 'k': 61.5, 'alpha95': 6.6},
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    lithology='sedimentary',
    component_comment='Primary hematite-carried characteristic remanent magnetization (west-normal / east-reverse), shallow; detrital/early-diagenetic',
    tests='',
    gpmdb_number='9290',
    percent_reversed=50,
    demag_code=4,
    R1=0, R2=1, R3=1, R4='', R5=0, R6=0, R7=0, Grade='B',
    nominal_age=754, lomagage=717, himagage=771,
    REF_method='Tonian depositional age: U-Pb detrital-zircon ash 742 +/- 6 Ma (Karlstrom et al., 2000) is a maximum depositional age midway through the section; correlation with the Chuar Group and Dehler et al. place the entire UMG <766 Ma and pre-717 Ma.',
    POLE_AUTHORS='Weil, A. B., Geissman, J. W., & Ashby, J. M.',
    YEAR=2006,
    JOURNAL='Precambrian Research',
    VOLUME='147',
    VPAGES='234-259',
    TITLE='A new paleomagnetic pole for the Neoproterozoic Uinta Mountain supergroup, Central Rocky Mountain States, USA',
    COMMENT='Uinta Mountain Group ChRM pole recreated as the Fisher mean of the 9 primary-ChRM locality VGPs of Weil et al. (2006) (audited MagIC contribution 20680: site result_type a->i; the 4 eastern present-field overprint localities relabeled and excluded). Reproduces the published pole (0.8N/161.3E, A95 4.6, GPMDB 9290). Sedimentary (detrital hematite) -> lithology=sedimentary blanket-f f-block. R=2 (Grade B) per the prior compilation: R1=0 (broad Tonian age), R5=0 (possible inclination shallowing), R7=0 (overlaps younger poles). Flag: dual-polarity ChRM may support R6=1; age likely refinable to ~720-760 Ma (Dehler/Karlstrom).'
)
pt.save_nordic_summary(umg_summary, '754_Uinta_Mountain_Group')
umg_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
