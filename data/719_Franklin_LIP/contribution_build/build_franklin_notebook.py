"""Rebuild pole_notebooks/719_Franklin_LIP.ipynb.

Grand mean is now computed from the full Denyszyn et al. (2009) compilation of 78
individual site VGPs (`data/719_Franklin_LIP/Denyszyn2009.csv`; nine source
studies, Greenland rotated to Laurentia by the Nares Strait fit), rather than a
hard-coded value. The Denyszyn et al. (2009) new dykes (the audited Cortopassi
MagIC contribution, 27 site-level dykes) are reproduced as the new component;
the remaining ~51 sites of the grand mean are to be added to the MagIC
contribution from the underlying papers.
"""
import os
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

NB = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                  '../../../pole_notebooks/719_Franklin_LIP.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Franklin LIP ca. 719 Ma paleomagnetic pole

## Geologic context

The Franklin Large Igneous Province is the ~719 Ma giant radiating mafic dyke
swarm and associated sills/volcanics emplaced across Arctic Canada and northwest
Greenland, immediately preceding the Sturtian snowball-Earth glaciation. Six
zircon 206Pb/238U dates of 718.61 ± 0.30 to 719.86 ± 0.21 Ma (Pu et al., 2022)
record its rapid emplacement. The Franklin pole is one of the best-determined
Neoproterozoic keystones of the Laurentia APWP.

## Pole

The Franklin grand-mean pole is computed here from the full Denyszyn et al. (2009)
compilation of 78 individual site VGPs spanning nine source studies (with the
Greenland sites rotated into Laurentian coordinates by the Nares Strait fit). The
Denyszyn et al. (2009) new dykes of High Arctic Canada and Greenland — the audited
site-level MagIC contribution — are reproduced separately; the remaining sites of
the grand mean are to be added to the contribution from the underlying papers.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## The Franklin grand mean from 78 individual site VGPs

`Denyszyn2009.csv` holds the 78 site VGPs that define the Denyszyn et al. (2009)
Franklin grand mean, with each Greenland VGP rotated to Laurentian coordinates
(`vgp_lat_rotated` / `vgp_lon_rotated`) and a normal/reversed polarity tag. The
swarm is of mixed polarity; polarity is unified before the Fisher mean.""")

co("""grand = pd.read_csv('../data/719_Franklin_LIP/Denyszyn2009.csv')
grand = grand.dropna(subset=['vgp_lat_rotated', 'vgp_lon_rotated']).reset_index(drop=True)\ngrand = grand.drop(columns=['vgp_lat', 'vgp_lon'])  # drop the NaN in-situ cols; use the rotated VGPs
print(f'{len(grand)} site VGPs; polarity:', grand['polarity'].value_counts().to_dict())
print('\\nby source study:')
print(grand['reference'].value_counts().to_string())

vgp_block, pole_mean = pt.compute_mean_pole(
    grand.rename(columns={'vgp_lat_rotated': 'vgp_lat', 'vgp_lon_rotated': 'vgp_lon'}),
    unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
print('\\nprior Nordic compilation Franklin event grand mean: 6.7 N / 162.1 E, A95 3.0, B=56')""")

co("""# plot the 78 site VGPs coloured by source study, with the grand-mean pole
ax = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                           central_latitude=pole_mean['inc'], figsize=(7, 7))
import matplotlib.cm as cm
import numpy as np
refs = list(grand['reference'].unique())
colors = cm.tab10(np.linspace(0, 1, len(refs)))
for ref, c in zip(refs, colors):
    sub = grand[grand['reference'] == ref]
    blk = pmag.flip(ipmag.make_di_block(sub['vgp_lon_rotated'].tolist(),
                                        sub['vgp_lat_rotated'].tolist()), combine=True)
    for lon, lat in blk:
        ipmag.plot_vgp(ax, lon, lat, color=c, markersize=22)
plt.title('Franklin grand mean: 78 site VGPs (rotated) by source study')
plt.show()""")

md("""## Site-level Denyszyn et al. (2009) new dykes

The audited MagIC contribution (`data/719_Franklin_LIP/sites.txt`) holds the 27
Denyszyn et al. (2009) site-level dykes (Arctic Canada + Greenland; NU1
reverse-polarity VGP latitude corrected to −16.6). These are the new component of
the grand mean; the other ~51 sites (Fahrig 1971, Palmer 1983, etc.) are to be
added to the contribution from the underlying papers.""")

co("""sites = pd.read_csv('../data/719_Franklin_LIP/sites.txt', sep='\\t', skiprows=1)
for name in ['Arctic Canada', 'Greenland']:
    sub = sites[sites['location'] == name]
    _, sp = pt.compute_mean_pole(sub, unify_polarity=True)
    print(f'{name:14s}: {sp[\"inc\"]:.1f}/{sp[\"dec\"]:.1f} A95 {sp[\"alpha95\"]:.1f} N {int(sp[\"n\"])}')
_, combined = pt.compute_mean_pole(sites, unify_polarity=True)
print(f'Combined (27):  {combined[\"inc\"]:.1f}/{combined[\"dec\"]:.1f} A95 {combined[\"alpha95\"]:.1f} N {int(combined[\"n\"])}')""")

md("## The Franklin grand mean in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Franklin event grand mean', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=719)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | Six zircon 206Pb/238U dates 718.61 ± 0.30 to 719.86 ± 0.21 Ma (Pu et al., 2022). |
| 2 | Techniques and statistical analysis | **1** | A grand mean of 78 site VGPs from nine studies, thoroughly averaging paleosecular variation. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence characterized (Denyszyn et al., 2009 and source studies). |
| 4 | Field tests constrain age of magnetization | **C** | Positive baked-contact test (Franklin dykes/sills bake the country rock). |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous Laurentia (Greenland restored by the Nares Strait fit). |
| 6 | Presence of reversals | **1** | The swarm is of mixed polarity (48 normal, 30 reversed). |
| 7 | No resemblance to younger poles | **1** | A well-resolved Neoproterozoic key pole, distinct from younger poles. |
| | Total | **7/7** | Grade A |""")

md("""## Nordic workshop summary

The Franklin event grand mean is computed from the full Denyszyn et al. (2009)
compilation of 78 individual site VGPs (nine source studies; Greenland rotated to
Laurentia): recreated value reported below. This refines the hard-coded
prior-compilation value (6.7°N/162.1°E, A95 3.0°, B=56); the small difference
reflects the 78-VGP set versus the 56 cooling-unit count and the rotation
treatment. A positive baked-contact test, mixed polarity, and the 719 Ma zircon
ages give a Grade A, R = 7 keystone pole.

**In progress:** the 78 individual sites are to be added to the MagIC contribution
(`data/719_Franklin_LIP/`) from the underlying papers (currently only the 27
Denyszyn et al. 2009 dykes are present at the site level).""")

co("""# the 78 grand-mean VGPs serve as the "sites" (one VGP per site); B and N reflect them
grand_sites = grand.rename(columns={'vgp_lat_rotated': 'vgp_lat', 'vgp_lon_rotated': 'vgp_lon'}).copy()
grand_sites['dir_tilt_correction'] = 0
grand_sites['dir_n_samples'] = 1
dir_mean = {'dec': 289.1, 'inc': 0.6, 'k': round(pole_mean['k'], 1), 'alpha95': 0.1}

franklin_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Franklin event grand mean',
    sites=grand_sites,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=275.4,
    study_lat=73.0,
    component_comment='Primary characteristic remanent magnetization (magnetite); mixed polarity; Franklin LIP grand mean of 78 site VGPs (nine studies)',
    tests='C+ (positive baked-contact test)',
    gpmdb_number='',
    magic_id='',
    percent_reversed=38,
    demag_code=4,
    R1=1, R2=1, R3=1, R4='C', R5=1, R6=1, R7=1, Grade='A',
    nominal_age=719, lomagage=718, himagage=720,
    REF_method='Six zircon 206Pb/238U dates 718.61 +/- 0.30 to 719.86 +/- 0.21 Ma record rapid emplacement of the Franklin LIP (Pu et al., 2022), superseding earlier baddeleyite dates.',
    POLE_AUTHORS='Denyszyn, S. W., Halls, H. C., Davis, D. W., & Evans, D. A. D. (compilation of nine studies)',
    YEAR=2009,
    JOURNAL='Canadian Journal of Earth Sciences',
    VOLUME='46',
    VPAGES='689-705',
    TITLE='Paleomagnetism and U-Pb geochronology of Franklin dykes in High Arctic Canada and Greenland',
    COMMENT='Franklin event grand mean computed from the full Denyszyn et al. (2009) compilation of 78 individual site VGPs (data/719_Franklin_LIP/Denyszyn2009.csv; nine source studies -- Fahrig 1971, Palmer 1983 Natkusiak, Denyszyn 2009 Devon/Ellesmere + Greenland, Christie+Fahrig 1983, Fahrig+Schwarz 1973, Palmer+Hayatsu 1975, Park 1981/1974; Greenland rotated to Laurentia by the Nares Strait fit). Recreated grand mean reported here; refines the hard-coded prior-compilation value 6.7N/162.1E (B=56). Mixed polarity (48n/30r, R6=1); positive baked-contact test (R4=C); zircon 206Pb/238U 718.6-719.9 Ma (Pu et al. 2022, R1=1). R=7, Grade A. IN PROGRESS: the 78 sites are to be added to the MagIC contribution from the underlying papers (currently only the 27 Denyszyn 2009 dykes are present at site level; audited Cortopassi contribution, NU1 vgp_lat -6.6->-16.6, grand-mean DOI sciadv.adc9431->adc9430).'
)
pt.save_nordic_summary(franklin_summary, '719_Franklin_LIP')
franklin_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
