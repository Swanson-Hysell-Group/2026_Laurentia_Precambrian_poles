"""Build the central North Greenland (Victoria Fjord) dolerite dyke pole (1382_Victoria).

Source: a student MagIC contribution (spreadsheet) for Abrahamsen & Van der Voo
(1987), "Palaeomagnetism of middle Proterozoic (c. 1.25 Ga) dykes from central
North Greenland", Geophys. J. R. astr. Soc. 91, 597-611
(doi:10.1111/j.1365-246X.1987.tb01660.x), audited against the paper. (New student,
no prior review; the xlsx carried several errors.)

Audited fixes (instructor review): method code LT-AF-I -> LP-DIR-AF:LP-DIR-T:
DE-BFL:DE-FM; site longitude 44.7 (E) -> 315.3 (= 44.7 W, the nunatak at
81.5 N / 44.7 W); D6 dir_dec 267 -> 275; D8 dir_n_samples 9 -> 5; D10 dir_alpha95
16.5 -> 18.5 (and D10 excluded from the mean, result_quality=b); lithology
"Dibase" -> "Diabase"; per-site VGPs recomputed with pmag.dia_vgp at the corrected
locality.

Abrahamsen & Van der Voo (1987): 9 dyke sites + 3 baked gneiss sites give D=265,
I=21.5, N=12, a95=5.6, confirmed by a positive baked-contact test. Observed pole
(before the Bullard rotation) ~12.5 N / 232 E. Single polarity (opposite to the
Zig-Zag Dal / Midsommerso units). Age: the paper infers ~1.25 Ga (Rb-Sr on
related intrusives); the compilation adopts 1382 +/- 2 Ma by correlation with the
U-Pb-dated Midsommerso Dolerites (Upton et al., 2005) -- a genuine age
uncertainty.
"""
import os
import openpyxl
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
XLSX = ('/Users/penokean/Dropbox/Teaching/ESCI_4204_8204/pole_project/MagIC_revised/'
        'ueckergabriel_LATE_1642971_62184407_contributespreaadsheetGabrieul uecker - .xlsx')

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
rows = list(wb['Sheet1'].iter_rows(values_only=True))
for i, r in enumerate(rows):
    if r and r[0] == 'site':
        hdr = [c for c in r if c is not None]
        start = i + 1
        break
recs = [list(r) for r in rows[start:] if r and r[0]]
sites = pd.DataFrame([rec[:len(hdr)] for rec in recs], columns=hdr)

# ---- audited fixes -----------------------------------------------------------
sites['lon'] = 315.3                                   # 44.7 W -> 0-360 E
sites['lithologies'] = 'Diabase'
sites['method_codes'] = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM'
sites['location'] = 'central North Greenland'
sites.loc[sites['site'] == 'D6', 'dir_dec'] = 275.0
sites.loc[sites['site'] == 'D8', 'dir_n_samples'] = 5
sites.loc[sites['site'] == 'D10', 'dir_alpha95'] = 18.5
sites.loc[sites['site'] == 'D10', 'result_quality'] = 'b'   # excluded from the mean
sites['dir_polarity'] = 'r'

def vgp_row(r):
    plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], r['dir_alpha95'], r['lat'], r['lon'])
    return pd.Series({'vgp_lon': round(plon, 1), 'vgp_lat': round(plat, 1),
                      'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1)})
sites = pd.concat([sites.drop(columns=['vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm']),
                   sites.apply(vgp_row, axis=1)], axis=1)

good = sites[sites['result_quality'] != 'b']
blk = pmag.flip(ipmag.make_di_block(good['vgp_lon'].tolist(), good['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)

locs = pd.DataFrame([{
    'location': 'central North Greenland', 'location_type': 'Region',
    'result_name': 'Central North Greenland dolerite dykes (Victoria Fjord) ca. 1382 Ma pole',
    'result_type': 'a', 'sites': ':'.join(good['site'].tolist()),
    'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP:ST-C',
    'citations': '10.1111/j.1365-246X.1987.tb01660.x:10.1007/s00410-004-0634-7',
    'geologic_classes': 'Intrusive', 'lithologies': 'Diabase',
    'lat_s': 81.5, 'lat_n': 81.5, 'lon_w': 315.3, 'lon_e': 315.3,
    'age': 1382, 'age_low': 1380, 'age_high': 1384, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1), 'pole_alpha95': round(p['alpha95'], 1),
    'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'description': 'Central North Greenland dolerite dyke pole (Abrahamsen & Van der Voo, 1987), nunatak at 81.5 N / 44.7 W. Observed (unrotated) pole; positive baked-contact test (3 baked gneiss sites carry the dyke direction). Single polarity, opposite to the Zig-Zag Dal / Midsommerso units. Age 1382+/-2 Ma adopted by correlation (Upton et al., 2005); the paper infers ~1.25 Ga (Rb-Sr).'}])

def write_magic(df, kind, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])}')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/1382_Victoria.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Central North Greenland (Victoria Fjord) dolerite dykes ca. 1382 Ma paleomagnetic pole

## Geologic context

These middle-Proterozoic basic dykes cut Archaean basement exposed as nunataks at
the head of Victoria Fjord in central North Greenland (81.5°N, 44.7°W;
Abrahamsen & Van der Voo, 1987). Their characteristic remanence is of a single
polarity that is antiparallel to that of the comagmatic Zig-Zag Dal Basalts and
Midsommerso Dolerites of eastern North Greenland, and confirms that Greenland was
part of the North American craton in the Mesoproterozoic. The dyke magnetization
is shown to be primary by a detailed positive baked-contact test.

## Pole

This notebook recreates the central North Greenland dolerite pole at the dyke
level from the 9 dyke sites of Abrahamsen & Van der Voo (1987), with the site
longitude and per-site VGPs corrected. The age is adopted as 1382 ± 2 Ma by
correlation with the U-Pb-dated Midsommerso Dolerites (Upton et al., 2005); the
paper itself infers ~1.25 Ga from Rb-Sr — a genuine age uncertainty.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data and the dyke pole

Loaded from `../data/1382_Victoria_Fjord/`. The pole uses the 9 accepted dyke
sites (D1-D9); D10 (a discordant easterly direction, α95 18.5°) is excluded.""")

co("""sites_geo, _ = pt.load_magic_sites('../data/1382_Victoria_Fjord/sites.txt')
study_lat, study_lon = 81.5, 315.3
print(f'{len(sites_geo)} accepted dyke sites')
sites_geo[['site', 'lat', 'lon', 'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k', 'dir_n_samples', 'vgp_lat', 'vgp_lon']]""")

co("""vgp_block, pole_mean = pt.compute_mean_pole(sites_geo, unify_polarity=True)
if pole_mean['inc'] < 0:
    vgp_block, pole_mean = pt.compute_mean_pole(sites_geo, unify_polarity=True, flip=True)
ipmag.print_pole_mean(pole_mean)
print('\\nprior compilation (Victoria Fjord, GPMDB 489): 10.3 N / 231.7 E, A95 5.9, N=12')
print('Abrahamsen & Van der Voo (1987) observed pole (before rotation): ~12.5 N / 232 E')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('Central North Greenland dolerite dyke pole')
plt.show()""")

co("""dir_block, dir_mean = pt.compute_mean_direction(sites_geo, unify_polarity=True)
ipmag.print_direction_mean(dir_mean)
ipmag.plot_net()
ipmag.plot_di(di_block=dir_block, color='blue', marker='o')
ipmag.plot_di_mean(dir_mean['dec'], dir_mean['inc'], dir_mean['alpha95'], color='red', marker='s')
plt.title('Central N Greenland dyke directions')
plt.show()""")

md("""## Field test: baked contact

Abrahamsen & Van der Voo (1987) sampled three baked gneiss sites adjacent to the
dykes; the baked Archaean gneiss carries the dyke direction, a positive
baked-contact test confirming the dyke magnetization is primary (the baked-gneiss
sites are not in this site-level table but are documented in the paper).""")

md("""## Paleosecular variation""")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

md("## The pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Victoria Fjord dolerite dykes', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=1382)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | Age adopted as 1382 ± 2 Ma by correlation with the U-Pb-dated Midsommerso Dolerites (Upton et al., 2005). (The paper infers ~1.25 Ga from Rb-Sr — see note.) |
| 2 | Techniques and statistical analysis | **0** | AF + thermal demagnetization; the VGPs are over-concentrated (K ≈ 104 ≫ 70, A95 below the Deenen et al. (2011) lower bound), indicating under-averaged paleosecular variation. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence characterized by rock-magnetic measurements (Abrahamsen & Van der Voo, 1987). |
| 4 | Field tests constrain age of magnetization | **C** | Positive baked-contact test (three baked gneiss sites carry the dyke direction). |
| 5 | Structural control / tectonic coherence | **0** | Isolated nunatak exposure; possible local rotation and the unresolved age leave the structural/tectonic coherence weak. |
| 6 | Presence of reversals | **0** | Single polarity within the swarm (antiparallel to the eastern North Greenland units, but no reversal within this collection). |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles. |
| | Total | **4/7** | Grade B |""")

md("""## Nordic workshop summary

The central North Greenland (Victoria Fjord) dolerite pole is recreated from the
9 dyke sites of Abrahamsen & Van der Voo (1987), after correcting the site
longitude (44.7°E → 315.3°E) and recomputing VGPs. The recreation
(≈12.7°N/231.3°E, A95 5.1°) reproduces the prior compilation Victoria Fjord pole
(10.3°N/231.7°E, GPMDB 489). Positive baked-contact test (R4 = C); over-
concentrated VGPs (R2 = 0).

**Age flag for review:** the paper infers ~1.25 Ga (Rb-Sr on related intrusives);
the compilation adopts 1382 ± 2 Ma by correlation with the U-Pb-dated Midsommerso
Dolerites — a genuine age uncertainty (R1 depends on this choice).""")

co("""victoria_summary = pt.make_nordic_summary(
    terrane='Laurentia-Greenland',
    rockname='Victoria Fjord dolerite dykes',
    sites=sites_geo,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Primary characteristic remanent magnetization (magnetite); single polarity (antiparallel to eastern North Greenland units)',
    tests='C+ (positive baked-contact test; three baked gneiss sites carry the dyke direction)',
    gpmdb_number='489',
    percent_reversed=0,
    demag_code=3,
    R1=1, R2=0, R3=1, R4='C', R5=0, R6=0, R7=1, Grade='B',
    nominal_age=1382, lomagage=1380, himagage=1384,
    REF_method='Age adopted as 1382 +/- 2 Ma by correlation with the U-Pb baddeleyite-dated Midsommerso Dolerites (Upton et al., 2005); Abrahamsen & Van der Voo (1987) infer ~1.25 Ga from Rb-Sr on related intrusives.',
    POLE_AUTHORS='Abrahamsen, N., & Van der Voo, R.',
    YEAR=1987,
    JOURNAL='Geophysical Journal of the Royal Astronomical Society',
    VOLUME='91',
    VPAGES='597-611',
    TITLE='Palaeomagnetism of middle Proterozoic (c. 1.25 Ga) dykes from central North Greenland',
    COMMENT='Central North Greenland (Victoria Fjord) dolerite pole recreated from the 9 dyke sites of Abrahamsen & Van der Voo (1987) (audited spreadsheet: site longitude 44.7E->315.3E, D6 dec 267->275, D8 n 9->5, D10 a95 16.5->18.5 + excluded, method code + lithology Dibase->Diabase fixed; VGPs recomputed). Recreated 12.7N/231.3E, A95 5.1, N=9 -- reproduces prior compilation Victoria Fjord pole (10.3/231.7, GPMDB 489). Positive baked-contact test (R4=C; the paper describes a normal baked-contact test, so C rather than the compilation R4=c). K~104>>70 -> R2=0; single polarity R6=0. AGE FLAG: paper infers ~1.25 Ga (Rb-Sr); compilation adopts 1382+/-2 Ma via correlation with Upton 2005 -- genuine uncertainty. R=4, Grade B.'
)
pt.save_nordic_summary(victoria_summary, '1382_Victoria_Fjord')
victoria_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
