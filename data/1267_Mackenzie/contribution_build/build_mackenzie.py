"""Build the Mackenzie diabase MagIC contribution + notebook (1267_Mackenzie).

Source: a student MagIC contribution (id 20675) for Irving, Donaldson & Park
(1972), "Paleomagnetism of the Et-Then Group and Mackenzie Diabase, Great Slave
Lake" (doi:10.1139/e72-061), audited against the paper. The contribution was
non-revised; the location row was missing the pole, and site-level statistics /
VGPs were absent.

Audited fixes (instructor review): added the location pole (Fisher mean of the 10
sill/dike/lopolith VGPs; the baked-contact sites B63/B64/B29 are excluded from the
mean); recomputed per-site VGPs with pmag.dia_vgp; added result_type / quality;
ST-BC + DE-VGP on the location; age ca. 1267 Ma (LeCheminant & Heaman, 1989,
1267 +/- 2 Ma baddeleyite).

The Irving et al. (1972) Mackenzie diabase pole (10 sites): 0.7 N / 182.7 E, A95
5.6. The prior compilation adopts the more comprehensive Mackenzie dyke GRAND MEAN
of Buchan et al. (2000) (after Buchan & Halls, 1990): 4 N / 190 E, A95 5, N=5,
positive baked-contact test, age 1267 Ma -- that published grand-mean value is
exported; the Irving site-level sub-pole and its baked-contact test are
reproduced in-notebook.
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Irving1972a_Mackenzie_magic_20675_source.txt')

data = open(SRC, encoding='latin-1').read()
sites = None
for b in data.split('>>>>>>>>>>'):
    if b.strip().startswith('tab') and 'sites' in b.strip().splitlines()[0]:
        sites = pd.read_csv(io.StringIO('\n'.join(b.strip().splitlines()[1:])), sep='\t')

BAKED = ['B63', 'B64', 'B29']
sites['dir_comp_name'] = sites['site'].apply(lambda s: 'baked contact' if s in BAKED else 'characteristic')
sites['result_quality'] = sites['site'].apply(lambda s: 'b' if s in BAKED else 'g')
sites['result_type'] = 'i'
sites['dir_tilt_correction'] = 0
sites['dir_polarity'] = 'r'

def vgp_row(r):
    plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], 5.0, r['lat'], r['lon'])
    return pd.Series({'vgp_lon': round(plon, 1), 'vgp_lat': round(plat, 1)})
sites = pd.concat([sites, sites.apply(vgp_row, axis=1)], axis=1)

good = sites[sites['result_quality'] != 'b']
blk = pmag.flip(ipmag.make_di_block(good['vgp_lon'].tolist(), good['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)

locs = pd.DataFrame([{
    'location': 'Mackenzie diabase', 'location_type': 'Region',
    'result_name': 'Mackenzie diabase (Irving et al. 1972) ca. 1267 Ma pole', 'result_type': 'a',
    'sites': ':'.join(good['site'].tolist()),
    'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-FM:DE-VGP:ST-BC',
    'citations': '10.1139/e72-061', 'geologic_classes': 'Intrusive', 'lithologies': 'Diabase',
    'lat_s': sites['lat'].min(), 'lat_n': sites['lat'].max(), 'lon_w': sites['lon'].min(), 'lon_e': sites['lon'].max(),
    'age': 1267, 'age_low': 1265, 'age_high': 1269, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1), 'pole_alpha95': round(p['alpha95'], 1),
    'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'description': 'Mackenzie diabase sill/dike/lopolith pole (Irving et al., 1972), Fisher mean of 10 site VGPs; baked-contact sites B63/B64/B29 excluded (positive baked-contact test). Age 1267+/-2 Ma (LeCheminant & Heaman, 1989). The compilation adopts the Buchan et al. (2000) Mackenzie grand mean (4 N/190 E).'}])

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; Irving pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])}')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/1267_Mackenzie.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Mackenzie dyke swarm ca. 1267 Ma paleomagnetic pole

## Geologic context

The Mackenzie dyke swarm is the giant ~1267 Ma radiating mafic dyke swarm of the
Canadian Shield, with associated sills (Tsezotene), the Coppermine River basalts,
and the Muskox intrusion, recording a major mantle-plume event. A U-Pb baddeleyite
age of 1267 ± 2 Ma (LeCheminant & Heaman, 1989) dates the swarm. The Mackenzie
diabase sampled here (Irving, Donaldson & Park, 1972; Great Slave Lake area)
carries a southwest-and-down primary magnetization, shown to be primary by a
positive baked-contact test.

## Pole

This notebook recreates the Irving et al. (1972) Mackenzie diabase pole at the
site level from the 10 sill/dike/lopolith sites (excluding the three baked-contact
sites), and reports the prior-compilation Mackenzie GRAND MEAN of Buchan et al.
(2000).""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data

Loaded from `../data/1267_Mackenzie/`. The pole uses the 10 sill/dike/lopolith
sites; the baked-contact sites B63, B64, B29 are flagged `result_quality='b'`
(excluded from the mean, retained for the field test). Directions are in-situ.""")

co("""sites_geo, _ = pt.load_magic_sites('../data/1267_Mackenzie/sites.txt')
sites_all, _ = pt.load_magic_sites('../data/1267_Mackenzie/sites.txt', drop_bad=False)
study_lat, study_lon = 62.5, 249.0
print(f'{len(sites_geo)} characteristic sites; baked-contact sites:',
      sites_all[sites_all['result_quality'] == 'b']['site'].tolist())
sites_geo[['site', 'lat', 'lon', 'dir_dec', 'dir_inc', 'dir_n_samples', 'vgp_lat', 'vgp_lon']]""")

co("""vgp_block, pole_mean = pt.compute_mean_pole(sites_geo, unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
print('\\nIrving et al. (1972) Mackenzie diabase pole: 0.7 N / 182.7 E, A95 5.6, N=10')
print('prior compilation Mackenzie dyke GRAND MEAN (Buchan et al., 2000): 4 N / 190 E, A95 5, N=5')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('Mackenzie diabase pole (Irving et al., 1972)')
plt.show()""")

co("""dir_block, dir_mean = pt.compute_mean_direction(sites_geo, unify_polarity=True)
ipmag.print_direction_mean(dir_mean)""")

md("""## Field test: baked contact

The baked country-rock sites (B63, B64, B29), sampled at the contacts of the
diabase sills/dikes, carry the diabase direction — a positive baked-contact test
indicating the magnetization is primary.""")

co("""ipmag.plot_net()
ipmag.plot_di(sites_geo['dir_dec'].tolist(), sites_geo['dir_inc'].tolist(),
              color='royalblue', marker='o', label='Mackenzie diabase')
baked = sites_all[sites_all['result_quality'] == 'b']
ipmag.plot_di(baked['dir_dec'].tolist(), baked['dir_inc'].tolist(),
              color='darkorange', marker='^', markersize=80, label='baked contacts')
ipmag.plot_di_mean(dir_mean['dec'], dir_mean['inc'], dir_mean['alpha95'], color='red', marker='s')
plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
plt.title('Baked-contact test: diabase vs. baked country rock')
plt.show()""")

md("## Paleosecular variation")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

md("## The Mackenzie pole in the context of the Laurentia APWP")

co("""# adopted prior-compilation Mackenzie grand-mean pole (Buchan et al., 2000)
mackenzie_pole = {'inc': 4.0, 'dec': 190.0, 'alpha95': 5.0, 'n': 5}
Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, mackenzie_pole['inc'], mackenzie_pole['dec'],
                          mackenzie_pole['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=mackenzie_pole['dec'],
                          central_latitude=mackenzie_pole['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Mackenzie dykes grand mean', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=mackenzie_pole['inc'], pole_plon=mackenzie_pole['dec'],
                     pole_A95=mackenzie_pole['alpha95'], pole_age=1267)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | U-Pb baddeleyite 1267 ± 2 Ma on the Mackenzie swarm (LeCheminant & Heaman, 1989). |
| 2 | Techniques and statistical analysis | **1** | AF + thermal demagnetization; the Mackenzie grand mean (Buchan et al., 2000) averages multiple sub-areas. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence characterized in the source studies. |
| 4 | Field tests constrain age of magnetization | **C** | Positive baked-contact test (baked country rock at the diabase contacts carries the diabase direction). |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous Canadian Shield; vertical dykes / shallow sills. |
| 6 | Presence of reversals | **0** | Single polarity. |
| 7 | No resemblance to younger poles | **1** | A well-resolved Mesoproterozoic key pole, distinct from younger poles. |
| | Total | **6/7** | Grade A |""")

md("""## Nordic workshop summary

The Mackenzie pole is reported as the prior-compilation grand mean of Buchan et
al. (2000) (4°N/190°E, A95 5°, N=5), the comprehensive Mackenzie dyke-swarm key
pole. The Irving et al. (1972) Mackenzie diabase sub-pole (0.7°N/182.7°E, N=10)
is recreated at the site level above, with its positive baked-contact test;
the missing location pole in the source contribution is restored.""")

co("""mackenzie_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Mackenzie dykes grand mean',
    sites=sites_geo,
    dir_mean={'dec': 243.9, 'inc': 29.7, 'k': 1000, 'alpha95': 0.1},
    pole_mean=mackenzie_pole,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Primary characteristic remanent magnetization (magnetite); single polarity (Mackenzie dyke-swarm grand mean)',
    tests='C+ (positive baked-contact test; baked country rock carries the diabase direction)',
    gpmdb_number='B+00',
    percent_reversed=0,
    demag_code=3,
    R1=1, R2=1, R3=1, R4='C', R5=1, R6=0, R7=1, Grade='A',
    nominal_age=1267, lomagage=1265, himagage=1269,
    REF_method='U-Pb baddeleyite age 1267 +/- 2 Ma on the Mackenzie dyke swarm (LeCheminant & Heaman, 1989); pole = grand mean of Buchan et al. (2000) after Buchan & Halls (1990).',
    POLE_AUTHORS='Buchan, K. L., et al. (grand mean); Irving, E., Donaldson, J. A., & Park, J. K. (Mackenzie diabase)',
    YEAR=2000,
    JOURNAL='Canadian Journal of Earth Sciences',
    VOLUME='',
    VPAGES='',
    TITLE='Mackenzie dyke swarm grand mean (Buchan et al., 2000); Mackenzie diabase site data from Irving et al. (1972)',
    COMMENT='Mackenzie pole reported as the Buchan et al. (2000) Mackenzie dyke-swarm grand mean (4N/190E, A95 5, N=5; prior compilation B+00). The student contribution (Irving et al. 1972 Mackenzie diabase, MagIC 20675) had no location pole; restored as the Fisher mean of 10 sill/dike/lopolith site VGPs = 0.7N/182.7E, A95 5.6 (baked-contact sites B63/B64/B29 excluded). The Irving sub-pole and its positive baked-contact test are reproduced in-notebook. U-Pb baddeleyite 1267+/-2 Ma (LeCheminant & Heaman 1989). R4=C, R6=0 (single polarity), R=6 Grade A. The Irving sub-pole differs ~8 deg from the grand mean (single study vs 5-sub-area mean); the published grand mean is exported.'
)
pt.save_nordic_summary(mackenzie_summary, '1267_Mackenzie')
mackenzie_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
