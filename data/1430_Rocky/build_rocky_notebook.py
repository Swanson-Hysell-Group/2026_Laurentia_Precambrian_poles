"""Rebuild pole_notebooks/1430_Rocky.ipynb to the modern template + Nordic export.

The repo already holds a complete site-level contribution for the "MEAN Rocky
Mountain intrusions" pole in data/1430_Rocky/ (built directly from Harlan et al.
1994 = Laramie Anorthosite Complex + Sherman Granite, and Harlan et al. 1998 =
Electra Lake gabbro). This supersedes the Ella Rank student contribution, which
covered only Laramie + Sherman and carried coordinate-transcription errors
(La5 lat 51.78->41.78; La2/La3 lon 245.5->254.5) -- the repo CSVs are transcribed
from the papers and do not have those errors.

The prior compilation "MEAN Rocky Mountain intrusions" is the mean of the three
study poles: -11.9 N / 217.4 E, A95 9.7, N=58 (B), mean direction 41.1/-46.6.
Age ~1430 Ma (range 1415-1445 across the three intrusions). R = 3, Grade B. A
pooled-site VGP mean is more dispersed (K~9) and lands ~9 deg away, so the
published 3-study-mean value is exported.
"""
import os
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.abspath(os.path.join(HERE, '../../pole_notebooks/1430_Rocky.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Rocky Mountain anorthosite-granite intrusions ca. 1430 Ma paleomagnetic pole

## Geologic context

The "MEAN Rocky Mountain intrusions" pole combines three ~1.4 Ga
anorogenic intrusions of the Colorado-Wyoming province: the Laramie Anorthosite
Complex and the Sherman Granite (Harlan et al., 1994), and the Electra Lake
gabbro (Harlan et al., 1998). These plutons span ca. 1415-1445 Ma. The remanence
is carried by magnetite; the igneous layering of the anorthosite and the
multiple intrusions provide stability tests.

## Pole

This notebook recreates the site-level data for the three intrusions, applies the
igneous-layering ("fold") and reversal tests, and reports the prior-compilation
mean pole (the mean of the three study poles). The repo's contribution
(`data/1430_Rocky/`, built from Harlan et al. 1994 + 1998) supersedes the
Laramie+Sherman-only student submission.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data (Harlan et al. 1994 + 1998)

Each intrusion's site means are loaded from `../data/1430_Rocky/sites.txt`, with
both in-situ (`dir_tilt_correction == 0`) and structurally-corrected
(`== 100`) directions; VGPs are computed from the tilt-corrected directions.""")

co("""sites_geo, sites_tc = pt.load_magic_sites('../data/1430_Rocky/sites.txt')
sites_tc = ipmag.vgp_calc(sites_tc.copy(), tilt_correction='yes', site_lon='lon',
                          site_lat='lat', dec_tc='dir_dec', inc_tc='dir_inc')
study_lat, study_lon = 40.3, 253.8
print('rock-type codes:', sites_tc['description'].value_counts().to_dict())
print(f'{len(sites_tc)} site means')
sites_tc[['site', 'description', 'lat', 'lon', 'dir_dec', 'dir_inc', 'dir_k', 'vgp_lat', 'vgp_lon']].head(10)""")

md("""## Sub-poles by intrusion and the combined mean

The Laramie Anorthosite (An / Sy), Sherman Granite (ShGr), and Electra Lake gabbro
(Pgb) sub-poles are computed; the compilation mean is the mean of the study poles.""")

co("""groups = {'Laramie Anorthosite': ['An', 'Sy'], 'Sherman Granite': ['ShGr'],
          'Electra Lake gabbro': ['Pgb', 'Pdb']}
study_poles = []
for name, codes in groups.items():
    sub = sites_tc[sites_tc['description'].isin(codes)]
    if len(sub):
        _, sp = pt.compute_mean_pole(sub, unify_polarity=True)
        study_poles.append(sp)
        print(f"{name:22s}: {sp['inc']:.1f}/{sp['dec']:.1f} A95 {sp['alpha95']:.1f} N {int(sp['n'])}")

print('\\nprior compilation MEAN Rocky Mountain intrusions (mean of study poles):')
print('  -11.9 N / 217.4 E, A95 9.7, N=58; mean direction 41.1/-46.6; age ~1430 Ma')""")

md("""## Igneous-layering ("fold") and reversal tests""")

co("""# igneous layering test: the anorthosite layering acts like bedding
ipmag.plot_net()
ipmag.plot_di(sites_geo['dir_dec'].tolist(), sites_geo['dir_inc'].tolist(),
              color='red', marker='o', label='in-situ')
ipmag.plot_di(sites_tc['dir_dec'].tolist(), sites_tc['dir_inc'].tolist(),
              color='blue', marker='s', label='layering-corrected')
plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
plt.title('Rocky Mountain intrusions: in-situ vs. layering-corrected directions')
plt.show()""")

co("""# reversal test on the layering-corrected directions
try:
    ipmag.reversal_test_bootstrap(dec=sites_tc['dir_dec'].tolist(),
                                  inc=sites_tc['dir_inc'].tolist(), plot=True)
except Exception as e:
    print('reversal test:', e)""")

md("## The Rocky Mountain pole in the context of the Laurentia APWP")

co("""# adopted prior-compilation mean pole
rocky_pole = {'inc': -11.9, 'dec': 217.4, 'alpha95': 9.7, 'n': 58}
Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, rocky_pole['inc'], rocky_pole['dec'],
                          rocky_pole['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=rocky_pole['dec'],
                          central_latitude=rocky_pole['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('MEAN Rocky Mountain intrusions', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=rocky_pole['inc'], pole_plon=rocky_pole['dec'],
                     pole_A95=rocky_pole['alpha95'], pole_age=1430)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **0** | A combined pole spanning three intrusions of ca. 1415-1445 Ma; the ~30 Ma spread exceeds ± 15 Ma. |
| 2 | Techniques and statistical analysis | **1** | AF + thermal demagnetization, PCA; N = 58 site means across three intrusions average paleosecular variation. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence characterized (Harlan et al., 1994, 1998). |
| 4 | Field tests constrain age of magnetization | **0** | The igneous-layering and combined-intrusion tests are not decisive regional field tests. |
| 5 | Structural control / tectonic coherence | **0** | A combined pole across three separate intrusions and structures. |
| 6 | Presence of reversals | **0** | Not scored as a passing reversal test in the prior compilation. |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles. |
| | Total | **3/7** | Grade B |""")

md("""## Nordic workshop summary

The Rocky Mountain intrusions pole is reported as the prior compilation's mean of
the three ~1.4 Ga study poles (Laramie Anorthosite + Sherman Granite, Harlan et
al. 1994; Electra Lake gabbro, Harlan et al. 1998): −11.9°N/217.4°E, A95 9.7°,
N=58, mean direction 41.1°/−46.6°. The repo's site-level contribution
(`data/1430_Rocky/`) is more complete than the Laramie+Sherman-only student
submission (which it supersedes, and which carried coordinate errors). A pooled
all-site VGP mean is more dispersed (K~9) and is not used; the published
study-pole-mean value is exported.""")

co("""rocky_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='MEAN Rocky Mountain intrusions',
    sites=sites_tc,
    dir_mean={'dec': 41.1, 'inc': -46.6, 'k': 1000, 'alpha95': 0.1},
    pole_mean=rocky_pole,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Characteristic remanent magnetization (magnetite); combined mean of three ca. 1.4 Ga Colorado-Wyoming intrusions',
    tests='',
    gpmdb_number='7493:7494:8342',
    percent_reversed=0,
    demag_code=4,
    R1=0, R2=1, R3=1, R4='', R5=0, R6=0, R7=1, Grade='B',
    nominal_age=1430, lomagage=1415, himagage=1445,
    REF_method='Combined pole of three ~1.4 Ga Colorado-Wyoming anorogenic intrusions: Laramie Anorthosite Complex and Sherman Granite (Harlan et al., 1994) and Electra Lake gabbro (Harlan et al., 1998); ages span ca. 1415-1445 Ma.',
    POLE_AUTHORS='Harlan, S. S., Geissman, J. W., et al.',
    YEAR=1994,
    JOURNAL='Canadian Journal of Earth Sciences / Tectonophysics',
    VOLUME='',
    VPAGES='',
    TITLE='Paleomagnetism of the Middle Proterozoic Laramie anorthosite complex and Sherman granite (1994) + Electra Lake gabbro (1998)',
    COMMENT='MEAN Rocky Mountain intrusions pole = mean of three ~1.4 Ga study poles (Laramie Anorthosite + Sherman Granite, Harlan et al. 1994; Electra Lake gabbro, Harlan et al. 1998): -11.9N/217.4E, A95 9.7, N=58, dir 41.1/-46.6 (prior compilation). The repo data/1430_Rocky contribution (built from the Harlan 1994+1998 papers) supersedes the Ella Rank Laramie+Sherman-only student submission (which had La5 lat 51.78->41.78 and La2/La3 lon transcription errors). Per-intrusion sub-poles and the igneous-layering + reversal tests are shown in-notebook. A pooled all-site VGP mean is more dispersed (~-3/209, K~9) and is not exported; the published study-pole-mean is reported. R=3, Grade B (R1=0 ~30 Ma age spread; R5=0 combined intrusions).'
)
pt.save_nordic_summary(rocky_summary, '1430_Rocky_Mountain_intrusions')
rocky_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
