"""Generate pole_notebooks/1108_Nipigon.ipynb from the site-level compilation."""
import json
from pathlib import Path

import nbformat as nbf

md = nbf.v4.new_markdown_cell
code = nbf.v4.new_code_cell

cells = []

cells.append(md(r"""# Nipigon (Logan) sills ca. 1108 Ma paleomagnetic pole

## Geologic context

The Nipigon sills (long known as the Logan sills) are reversely magnetized mafic sills
and related intrusions emplaced into the Mesoproterozoic Sibley Group and older
Animikie and Archean rocks of the Nipigon Embayment and the Thunder Bay district on
the north shore of Lake Superior. They were intruded during the earliest magmatic
stage of the Midcontinent Rift, contemporaneous with the lower Osler Volcanic Group,
and are among the most widely sampled units of the "Keweenawan" reverse-polarity
interval. Sill thicknesses reach several hundred metres and the sills are subhorizontal
to gently dipping toward the rift axis.

## Paleomagnetic data

Earlier versions of this compilation carried the Lulea Working Group (2009) eight-unit
grand mean as reported in Evans et al. (2021), for which the underlying per-unit data
were not digitized. This entry instead recomputes the pole from a site-level compilation
of the published results (E. J. Iloranta, 2026), drawing on DuBois (1962), Robertson and
Fahrig (1971), Pesonen (1979), Middleton et al. (2004), Borradaile and Middleton (2006),
and the Thunder Intrusion result of Piispa et al. (in preparation). VGPs were recomputed
from the published site mean directions and site coordinates rather than adopted from the
original papers, and reverse-polarity VGPs are reported as the antipodal north pole.

Because many of the same sills were sampled by more than one study, and because
Middleton et al. (2004) list the same locality repeatedly under different demagnetization
treatments, the published results are averaged in two stages into independent cooling
units before the pole is computed. Two further same-sill correlations established from the
source publications are applied, and three that the publications do not resolve are
reported below as a sensitivity test. The construction of the site table, the merging
rules, and the exclusions are documented in
`../data/1108_Nipigon/contribution_build/build_nipigon_contribution.py`.

## Age

The Nipigon (Logan) sills cluster in age near 1105-1109 Ma. The Nipigon sills are dated
at 1108.2 +/- 0.9 Ma (zircon upper intercept, recalculated from the data of Davis and
Sutcliffe, 1985), and Bleeker et al. (2020) report new high-precision U-Pb dates of
1106.3 +/- 2.0 Ma for the main Logan Sill at Mount McKay and 1105.5 +/- 3.0 Ma for the
Inspiration sill. The pole is assigned an age of 1108 +/- 2 Ma. The older Heaman et al.
(2007) baddeleyite dates are imprecise and superseded and are not used here.
"""))

cells.append(code("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline"""))

cells.append(md(r"""## Site-level data

`../data/1108_Nipigon/sites.txt` holds one row per independent cooling unit. Each row is
either a single published site mean or the Fisher mean of several determinations of the
same sill (the `description` field records which). The sills are intrusive and are treated
in geographic coordinates (`dir_tilt_correction = 0`), so the geographic frame is used
throughout. Results excluded from the pole are retained in the file with
`result_quality = 'b'` and are dropped on load.
"""))

cells.append(code("""units, _ = pt.load_magic_sites('../data/1108_Nipigon/sites.txt')
all_rows, _ = pt.load_magic_sites('../data/1108_Nipigon/sites.txt', drop_bad=False)

print(f'{len(units)} independent units entering the pole '
      f'({units["dir_n_samples"].sum()} samples)')
print(units['location'].value_counts().to_string())
print()
print('excluded and retained in the contribution:')
print(all_rows[all_rows['result_quality'] == 'b'][['site', 'description']].to_string(index=False))"""))

cells.append(md(r"""The published results behind each unit, including the excluded ones and the exclusion
reason, are kept alongside the build script.
"""))

cells.append(code("""results = pd.read_csv('../data/1108_Nipigon/contribution_build/accepted_results.csv')
merged = results[results['accepted']].groupby('unit').filter(lambda g: len(g) > 1)
print('units built by merging more than one published determination:')
print(merged[['unit', 'site', 'Reference', 'dir_dec', 'dir_inc',
              'dir_alpha95', 'dir_n_samples']].to_string(index=False))"""))

cells.append(code("""pt.plot_site_map(units, zoom_start=7)"""))

cells.append(md(r"""## The Nipigon sills pole

The pole is the Fisher mean of the unit VGPs.
"""))

cells.append(code("""STUDY_LAT, STUDY_LON = units['lat'].mean(), units['lon'].mean()

vgp_block, pole_mean = pt.compute_mean_pole(units, unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(5, 5))"""))

cells.append(code("""dir_block, dir_mean = pt.compute_mean_direction(units, unify_polarity=False)
ipmag.print_direction_mean(dir_mean)
ipmag.plot_net()
ipmag.plot_di(di_block=dir_block, color='blue', marker='o')
ipmag.plot_di_mean(dir_mean['dec'], dir_mean['inc'], dir_mean['alpha95'], color='red', marker='s')
plt.title('Nipigon sills — unit mean directions (geographic)')
plt.show()"""))

cells.append(code("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
_ = pt.plot_Deenen_test(pole_mean)"""))

cells.append(code("""fishqq_result = pt.fishqq_vgps(units, unify_polarity=True)
fishqq_result"""))

cells.append(code("""svei_result = pt.svei_test_vgps(units, STUDY_LON, STUDY_LAT, model='TK03_GAD', plot=True)
print(f"paleolatitude = {svei_result['lat']:.1f}\\u00b0; elongation E = {svei_result['E']:.2f} "
      f"({'consistent' if svei_result['E_result'] else 'inconsistent'} with TK03.GAD)")"""))

cells.append(md(r"""The VGP distribution is consistent with a Fisher distribution. The elongation is
higher than the TK03.GAD prediction, which is expected for a compilation of sills of
heterogeneous quality drawn from six studies over a ~200 km wide region rather than a
single stratigraphically controlled sequence; the elongation direction is not
interpreted here.
"""))

cells.append(md(r"""## Comparison between regions

The compiled sills fall into two geographic groups: those of the Nipigon Embayment
proper, and those of the Thunder Bay-Logan area to the southwest (including the Sibley
Peninsula and the Mount McKay sill). Both were emplaced during the same short-lived
magmatic episode, so the two sub-poles provide an internal consistency check.
"""))

cells.append(code("""region_results = {}
for region in ['Nipigon Embayment', 'Thunder Bay-Logan']:
    sub = units[units['location'] == region].reset_index(drop=True)
    _, p = pt.compute_mean_pole(sub, unify_polarity=True)
    _, d = pt.compute_mean_direction(sub, unify_polarity=False)
    region_results[region] = dict(sub=sub, pole=p, dir=d)
    print(f'{region:20s} N={p["n"]:3d}  plat={p["inc"]:5.1f}  plon={p["dec"]:6.1f}  '
          f'A95={p["alpha95"]:5.1f}  K={p["k"]:5.1f}')"""))

cells.append(code("""ax = pt.plot_apwp_context(pt.get_Laurentia_poles(), pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=1050, age_max=1200,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
for region, color in [('Nipigon Embayment', 'darkorange'), ('Thunder Bay-Logan', 'seagreen')]:
    p = region_results[region]['pole']
    ipmag.plot_pole(ax, p['dec'], p['inc'], p['alpha95'], color=color,
                    markersize=40, label=f'{region} (N={p["n"]})')
ax.legend(loc='lower left')
plt.title('Regional sub-poles within the Nipigon sills compilation')
plt.show()"""))

cells.append(md(r"""The two regional means overlap within their confidence circles, so the compilation is
treated as a single population. The Thunder Bay-Logan mean is displaced to lower pole
longitude, in the direction of the younger part of the Keweenawan track; whether that
reflects a real age difference between the Mount McKay-area sills (1106.3 +/- 2.0 Ma;
Bleeker et al., 2020) and the Nipigon Embayment sills, or simply the lower precision of
that smaller group, is not resolved by these data.
"""))

cells.append(md(r"""## Comparison between studies

Agreement between studies that used very different laboratory methods over four decades
is a further check that the compiled directions record a single magnetization.
"""))

cells.append(code("""STUDIES = {'10.4095/100589': 'DuBois (1962)',
           '10.1139/e71-125': 'Robertson & Fahrig (1971)',
           '10.17741/bgsf/51.1-2.004': 'Pesonen (1979)',
           '10.1029/2003JB002581': 'Middleton et al. (2004)',
           '10.1016/j.precamres.2005.10.007': 'Borradaile & Middleton (2006)',
           'Piispa et al., in preparation': 'Piispa et al. (in prep.)'}

for doi, name in STUDIES.items():
    sub = units[units['citations'].str.contains(doi, regex=False)]
    if len(sub) < 3:
        print(f'{name:32s} N={len(sub):3d}  (too few units for a mean)')
        continue
    _, p = pt.compute_mean_pole(sub, unify_polarity=True)
    print(f'{name:32s} N={p["n"]:3d}  plat={p["inc"]:5.1f}  plon={p["dec"]:6.1f}  '
          f'A95={p["alpha95"]:5.1f}  K={p["k"]:5.1f}')"""))

cells.append(md(r"""## Are any sills still counted twice?

Each sill should contribute a single VGP. Two same-sill correlations beyond those
identified in the source compilation are established by the source publications and are
applied in the site table:

* **DuBois N29 (the sill at Red Rock) and N27-N28 (the sediment it baked).** DuBois's
  Table XIII places all three at one cliff above the CPR and CNR tracks at Red Rock, with
  N27-N28 "from sediment immediately below diabase sill" and N29 "from lower part of the
  sill". A sill and the sediment it baked record the same cooling event, so both join the
  Red Rock North unit.
* **Robertson and Fahrig S7 and DuBois's Mount Mackay sites.** S7 plots 0.33 km from the
  lower Mount McKay sill at Fort William, the exposure dated by Bleeker et al. (2020) at
  1106.3 +/- 2.0 Ma. DuBois obtained badly scattered directions there, attributed them to
  a superimposed random component and resorted to great-circle fits (his Table XV); S7
  (k = 92) is much the better determination of that sill.

A third candidate is **refuted**. Pesonen R47 and R49 carry identical coordinates and
identical notes ("North shore of Nipigon Bay, same sill as R&F S14&15?") in the source
compilation, which suggested they might be one sill. The Geological Survey of Finland
report version of Pesonen (1979) (Report 1942, Q 20/27.2/79/1) states that the 270 hand
samples came "from 18 reversed and 22 normal dikes and from 10 reversed sills", and its
Figure 1 plots R47 and R49 as two separate sill sites on the north shore of Nipigon Bay
with R49 slightly north of R47. Sites R40-R49 are ten distinct sills by Pesonen's own
count, so they are kept separate; the shared coordinates are a digitizing artifact to
correct in the source compilation.

Two candidates remain unresolved. Robertson and Fahrig (1971) describe none of their
sites, so their site positions here are digitized from that paper's Figure 1.

| Candidate | Basis | Status |
|---|---|---|
| Robertson and Fahrig S8 with DuBois N19-21 | 3.5 km apart on the Highway 11 scarp south of Orient Bay | needs the map |
| DuBois N16-18 with N9-11 (Doghead Mountain) | DuBois locates N16 "2 miles southwest of Ozone Station" and N9-N11 "just south of Ozone Station" | needs the map |

A structural caveat limits how far this can be pushed: the Logan sills are laterally
extensive sheets that transgress from one bedding plane to another (Robertson and Fahrig,
1971), so exposures tens of kilometres apart may belong to one cooling unit. The source
compilation already merges Pesonen R45 into the Terry Fox unit on that basis: Pesonen's
Figure 1 places R45 near Pass Lake at the neck of the Sibley Peninsula, ~27 km northeast
of the Current River / Terry Fox exposure sampled by DuBois, Robertson and Fahrig and
Middleton et al., so the merge asserts that one sheet is continuous over that distance.
Applied consistently, such map-based correlation would likely reduce the unit count
further; N = 37 is therefore an upper bound on the number of independent cooling units.
"""))

cells.append(code("""CANDIDATE_MERGES = {'S8': 'N19-21', 'N16-18': 'Doghead Mountain'}


def pole_with_merges(merges=(), drop=()):
    \"\"\"Recomputes the pole after further merging or dropping named units.

    Units named in ``merges`` are folded into their target unit by taking the
    Fisher mean of the two units' VGPs before the pole is computed, so that a
    sill counted twice contributes a single VGP.

    Args:
        merges (dict or tuple): Mapping of unit name to the unit it merges into.
        drop (iterable): Unit names to exclude entirely.

    Returns:
        dict: Fisher mean pole from ``ipmag.fisher_mean``.
    \"\"\"
    u = units[~units['site'].isin(drop)].copy()
    u['merged_site'] = u['site'].map(dict(merges)).fillna(u['site'])
    rows = []
    for _, g in u.groupby('merged_site'):
        m = ipmag.fisher_mean(dec=list(g['vgp_lon']), inc=list(g['vgp_lat']))
        rows.append((m['dec'], m['inc']))
    return ipmag.fisher_mean(dec=[r[0] for r in rows], inc=[r[1] for r in rows])


def row(label, **kw):
    p = pole_with_merges(**kw)
    return {'selection': label, 'N': p['n'], 'plat': round(p['inc'], 1),
            'plon': round(p['dec'], 1), 'A95': round(p['alpha95'], 1),
            'K': round(p['k'], 1)}

same_sill = pd.DataFrame([
    row('preferred (applied merges only)'),
    row('+ S8 into N19-21', merges={'S8': 'N19-21'}),
    row('+ N16-18 into Doghead Mountain', merges={'N16-18': 'Doghead Mountain'}),
    row('+ both candidate merges', merges=CANDIDATE_MERGES),
    row('+ both merges, refuted R49 merge applied anyway', merges=dict(CANDIDATE_MERGES, **{'R49': 'Kama Hill'})),
])
same_sill"""))

cells.append(md(r"""The pole moves by less than half a degree under any combination of the
unresolved merges, so the outcome of the map check will not change the result; what it
changes is the honest value of N and hence A95. K rises as duplicated determinations are
removed, which is the expected signature of collapsing non-independent readings of the
same cooling event.
"""))

cells.append(md(r"""## Sensitivity of the pole to compilation choices

Several units carry caveats from the source studies. Robertson and Fahrig's S5 and S6
("Colville Lake") may belong to the younger *ca.* 1100.8 Ma McIntyre diabase; the Seagull
Intrusion is dated at 1112 +/- 2.4 Ma (Hart and Whaley, 2005, as reported by Borradaile
and Middleton, 2006) and so may predate the sills; and several units have very poorly
determined mean directions (DuBois's Mount McKay and P8m results, for which he invoked
great-circle fits and noted possible lightning remagnetization, and Middleton's Hele 2 and
Havoc sites, whose thermal and low-temperature treatments give discordant directions).
"""))

cells.append(code("""def summarize(sub, label):
    _, p = pt.compute_mean_pole(sub, unify_polarity=True)
    return {'selection': label, 'N': p['n'], 'plat': round(p['inc'], 1),
            'plon': round(p['dec'], 1), 'A95': round(p['alpha95'], 1),
            'K': round(p['k'], 1)}

sensitivity = pd.DataFrame([
    summarize(units, 'all units (preferred)'),
    summarize(units[units['site'] != 'Colville Lake'], 'drop Colville Lake (possible McIntyre diabase)'),
    summarize(units[units['site'] != 'Seagull Intrusion'], 'drop Seagull Intrusion (ca. 1112 Ma)'),
    summarize(units[units['dir_alpha95'] <= 30], 'drop units with mean-direction alpha95 > 30 deg'),
    summarize(units[units['dir_alpha95'] <= 20], 'drop units with mean-direction alpha95 > 20 deg'),
    summarize(units[units['site'] != 'P8m'], 'drop P8m (2 discordant samples, alpha95 138 deg)'),
    summarize(units[~units['site'].isin(['P8m', 'N31m'])], 'drop P8m and N31m (N31m is a dyke cutting the Sleeping Giant sill)'),
    summarize(units[units['location'] == 'Nipigon Embayment'], 'Nipigon Embayment only'),
    summarize(units[units['location'] == 'Thunder Bay-Logan'], 'Thunder Bay-Logan only'),
])
sensitivity"""))

cells.append(md(r"""The pole position is stable at the few-degree level across all of these choices, well
within the A95 of the preferred mean, so no unit is excluded on these grounds. Note that
for units built by merging determinations from several studies the tabulated alpha95
describes the scatter *between* those determinations rather than within a site, so the
alpha95 filters preferentially remove well-replicated sills and are a conservative rather
than a sharpening test.
"""))

cells.append(md(r"""## Statistical assessment (Meert et al., 2020; R2)"""))

cells.append(code("""_ = pt.assess_R2(units, pole_mean)"""))

cells.append(md(r"""## Position within the Laurentia apparent polar wander path"""))

cells.append(code("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1300,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()"""))

cells.append(code("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
# 'MEAN Nipigon sills and lavas' is this unit's ROCKNAME in Laurentia_poles.csv,
# so naming it here excludes the superseded entry from the younger-pole comparison
_ = pt.plot_pole_overlap('MEAN Nipigon sills and lavas', Laurentia_stricto_poles,
                         pt.Torsvik2012_Laurentia,
                         pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                         pole_A95=pole_mean['alpha95'], pole_age=1108)
plt.show()"""))

cells.append(md(r"""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | Nipigon sills 1108.2 ± 0.9 Ma (recalculated from Davis & Sutcliffe, 1985); Logan Sill at Mount McKay 1106.3 ± 2.0 Ma and Inspiration sill 1105.5 ± 3.0 Ma (Bleeker et al., 2020). |
| 2 | Techniques and statistical analysis | **1** | See the R2 assessment above: N = 37 independent units, K = 21.8 and A95 = 5.2° within the Deenen et al. (2011) envelope, and the unit VGPs are consistent with a Fisher distribution. Stepwise AF and thermal demagnetization with PCA in Pesonen (1979), Middleton et al. (2004) and Borradaile & Middleton (2006); DuBois (1962) and Robertson & Fahrig (1971) used blanket cleaning only. Individual determinations are of mixed quality, but the pole is stable to the exclusion of the poorest units and the four study means agree within uncertainty. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite (with ilmenite intergrowths and minor pyrrhotite) identified by Robertson & Fahrig (1971); rock-magnetic and opaque-mineral characterization in Middleton et al. (2004) and Borradaile & Middleton (2006). |
| 4 | Field tests constrain age of magnetization | **c** | Positive inverse baked-contact tests: baked Sibley Group and Rove Formation adjacent to the reversed sills carry the sill direction while unbaked host rock does not (Pesonen, 1979, Table 4). |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous Superior craton; subhorizontal sills, geographic coordinates. |
| 6 | Presence of reversals | **0** | Single (reverse) polarity. |
| 7 | No resemblance to younger poles | **1** | On the Keweenawan track and distinct from the younger normal-polarity Midcontinent Rift poles. |
| | Total | **6/7** | Grade A |
"""))

cells.append(md(r"""## Nordic workshop summary

The pole is exported in the Nordic compilation format. In contrast to the previous entry,
which carried the Lulea Working Group grand mean with sentinel values for the site and
sample counts and the direction precision, the exported statistics are now the recomputed
site-level quantities.

| | This study (unit VGPs) | Previous Nordic compilation |
|---|---|---|
| Component | Nipigon (Logan) sills, ChRM | MEAN Nipigon sills and lavas |
| N (units / sites) | 37 | 86 |
| Pole lat (°N) | see below | 47.2 |
| Pole lon (°E) | see below | 217.8 |
| A95 (°) | see below | 4.0 |
"""))

cells.append(code("""nipigon_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Nipigon sills',
    sites=units,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=STUDY_LON,
    study_lat=STUDY_LAT,
    component_comment=('Single reverse-polarity ChRM of the Nipigon (Logan) sills and '
                       'related intrusions; Fisher mean of independent unit VGPs '
                       'recompiled from DuBois (1962), Robertson & Fahrig (1971), '
                       'Pesonen (1979), Middleton et al. (2004), Borradaile & Middleton '
                       '(2006) and Piispa et al. (in prep.)'),
    tests='c+ (positive inverse baked-contact test on Sibley Group and Rove Formation; Pesonen, 1979)',
    gpmdb_number='',
    magic_id='',
    percent_reversed=100,
    demag_code=3,
    R1=1, R2=1, R3=1, R4='c', R5=1, R6=0, R7=1, Grade='A', Grade_E21='A',
    nominal_age=1108, lomagage=1106, himagage=1110,
    REF_method=('The Nipigon (Logan) sills cluster in age near 1105-1109 Ma. The Nipigon '
                'sills are dated at 1108.2 +/- 0.9 Ma (zircon upper intercept, '
                'recalculated from the data of Davis and Sutcliffe, 1985), and Bleeker et '
                'al. (2020) report new high-precision U-Pb dates of 1106.3 +/- 2.0 Ma for '
                'the main Logan Sill at Mount McKay and 1105.5 +/- 3.0 Ma for the '
                'Inspiration sill. The pole is assigned an age of 1108 +/- 2 Ma. The older '
                'Heaman et al. (2007) baddeleyite dates are imprecise and superseded and '
                'are not used here.'),
    POLE_AUTHORS=('Iloranta, E. J. compilation of DuBois (1962), Robertson & Fahrig (1971), '
                  'Pesonen (1979), Middleton et al. (2004), Borradaile & Middleton (2006) '
                  'and Piispa et al. (in prep.)'),
    YEAR=2026,
    JOURNAL='',
    VOLUME='',
    VPAGES='',
    TITLE='Nipigon (Logan) sills ca. 1108 Ma pole, recomputed from site-level data',
    COMMENT=('Replaces the Lulea Working Group (2009) eight-unit grand mean carried in '
             'Evans et al. (2021) (47.2 N / 217.8 E, A95 4.0, B = 86) with a pole '
             'recomputed from published site means. Repeat determinations of the same '
             'sill are averaged into independent cooling units in two stages: the 15 '
             'entries of Middleton et al. (2004) Table 1 collapse to their 9 sampling '
             'localities, and sills sampled by more than one study are merged. Pillar '
             'Lake Lava, three normal-polarity Middleton et al. (2004) entries '
             'interpreted as incompletely removed overprints, and two results flagged as '
             'excluded in the source compilation are omitted. Geographic coordinates '
             '(intrusive, no tilt correction).')
)
pt.save_nordic_summary(nipigon_summary, '1108_Nipigon')
nipigon_summary"""))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata = {
    'kernelspec': {'display_name': 'Python [conda env:miniforge3-ess-jbook]',
                   'language': 'python',
                   'name': 'conda-env-miniforge3-ess-jbook-py'},
    'language_info': {'name': 'python'},
}
out = Path('/Users/hematite/Documents/GitHub/2026_Laurentia_Precambrian_poles/pole_notebooks/1108_Nipigon.ipynb')
nbf.write(nb, str(out))
print('wrote', out, len(cells), 'cells')
