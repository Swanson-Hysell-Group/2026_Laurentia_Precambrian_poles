"""Build the Michikamau Intrusion MagIC contribution + notebook (1460_Michikamau).

Source: a student MagIC contribution (id 20668) for the Michikamau anorthosite
pole, combining Emslie, Irving & Park (1976), CJES 13, 1052-1057
(doi:10.1139/e76-108) with Murthy, Fahrig & Jones (1968), CJES 5, 1139-1144
(doi:10.1139/e68-111), audited against the papers.

Audited fixes (instructor review): Michikamau location citations corrected to
10.1139/e76-108:10.1139/e68-111 (was the wrong e65-030 = Emslie 1965); the
Petscapiskau (site 22 meta-andesite) VGP corrected to -12/254 (the source carried
the polarity-flipped antipode 11.9/73.9); ST-C added to the Michikamau location
method codes (the baked-contact test supports the Michikamau primary
magnetization); DE-DI removed in favor of DE-VGP.

Emslie et al. (1976) Table 2: combined 12 localities / 54 cores give direction
259.5/+11 (k=43, a95=6.5) and pole 1.5 S / 142 W (=-1.5/218), A95 4.5, dp 3.5/dm 7.
Remanence carried by titanium-free magnetite + hematite; positive baked-contact
(Bruhnes) test (site 13 inside the aureole matches the intrusion; site 22
Petscapiskau outside differs). Age of magnetization ~1450-1460 Ma (U-Pb zircon
1460 Ma on adamellite; Krogh & Davis, 1973).
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Emslie1976_Murthy1968_Michikamau_magic_20668_source.txt')

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

# fix Petscapiskau site VGP antipode -> -12/254
pet = sites['location'].str.contains('Petscapiskau', na=False)
sites.loc[pet, 'vgp_lat'] = -12.0
sites.loc[pet, 'vgp_lon'] = 254.0

# Michikamau location: correct citations, add ST-C, prefer DE-VGP
mich = locs['location'].str.contains('Michikamau', na=False)
locs.loc[mich, 'citations'] = '10.1139/e76-108:10.1139/e68-111'
locs.loc[mich, 'method_codes'] = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP:ST-C'

# recompute Michikamau pole from the 12 intrusion site VGPs
mi = sites[sites['location'] == 'Michikamau Intrusion']
blk = ipmag.make_di_block(mi['vgp_lon'].tolist(), mi['vgp_lat'].tolist())
p = pmag.fisher_mean(blk)
locs.loc[mich, 'pole_lat'] = round(p['inc'], 1)
locs.loc[mich, 'pole_lon'] = round(p['dec'], 1)
locs.loc[mich, 'pole_alpha95'] = round(p['alpha95'], 1)
locs.loc[mich, 'pole_k'] = round(p['k'], 1)
locs.loc[mich, 'pole_n_sites'] = int(p['n'])
locs.loc[mich, 'age'] = 1460

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; Michikamau pole '
      f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])}')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/1460_Michikamau.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Michikamau Intrusion ca. 1460 Ma paleomagnetic pole

## Geologic context

The Michikamau Intrusion of central Labrador is a large, weakly to strongly
layered leucotroctolite-anorthosite massif with minor late adamellite, surrounded
by a thermal contact aureole up to several kilometres wide developed in pelitic
rocks of the older Petscapiskau Group (Emslie, Irving & Park, 1976). The intrusion
has not been regionally deformed or metamorphosed since consolidation. A U-Pb
zircon age of ~1460 Ma (Krogh & Davis, 1973) on the youngest adamellite dates
emplacement; the magnetization age is ~1450-1460 Ma. The remanence is carried by
fine-grained, titanium-free magnetite and hematite.

## Pole

This notebook recreates the Michikamau pole at the site level from the 12
collecting localities of Murthy et al. (1968) and Emslie et al. (1976). The
intrusion's westerly, shallow magnetization is shown to be primary by a positive
baked-contact (Bruhnes) test: country rock within the thermal aureole (Emslie's
site 13) carries the intrusion direction, while the older Petscapiskau Group
outside the aureole (site 22) carries a distinct northeast direction.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data

Loaded from `../data/1460_Michikamau/`. The pole uses the 12 Michikamau Intrusion
localities; the Petscapiskau Group meta-andesite (the baked-contact reference
outside the aureole) is held separately. Directions are in-situ (intrusive rocks).""")

co("""sites_geo, _ = pt.load_magic_sites('../data/1460_Michikamau/sites.txt')
michikamau = sites_geo[sites_geo['location'] == 'Michikamau Intrusion'].reset_index(drop=True)
petscapiskau = sites_geo[sites_geo['location'].str.contains('Petscapiskau', na=False)].reset_index(drop=True)
study_lat, study_lon = 54.5, 296.0
print(f'{len(michikamau)} Michikamau localities; {len(petscapiskau)} Petscapiskau reference')
michikamau[['site', 'lat', 'lon', 'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k',
            'dir_n_samples', 'vgp_lat', 'vgp_lon']]""")

md("""## Mean pole from site VGPs

The intrusion carries a single westerly, shallow magnetization; the Fisher mean
of the 12 locality VGPs gives the pole.""")

co("""vgp_block, pole_mean = pt.compute_mean_pole(michikamau, unify_polarity=False)
ipmag.print_pole_mean(pole_mean)
print('\\nEmslie et al. (1976): pole 1.5 S / 142 W (= -1.5 N, 218 E), A95 4.5, N=12')
print('prior compilation (GPMDB 2274): -1.5 N / 217.5 E, A95 4.7')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('Michikamau Intrusion pole')
plt.show()""")

co("""dir_block, dir_mean = pt.compute_mean_direction(michikamau, unify_polarity=False)
ipmag.print_direction_mean(dir_mean)""")

md("""## Field test: baked contact (Bruhnes test)

The country rock within the thermal aureole (Emslie's site 13) carries the
westerly intrusion direction, whereas the older Petscapiskau Group meta-andesite
outside the aureole (site 22) carries a distinct northeast-and-up direction — a
positive baked-contact test indicating the intrusion's magnetization is primary
(acquired on cooling). The intrusion mean and the Petscapiskau reference
direction are compared below.""")

co("""ipmag.plot_net()
ipmag.plot_di(michikamau['dir_dec'].tolist(), michikamau['dir_inc'].tolist(),
              color='royalblue', marker='o', label='Michikamau intrusion')
if len(petscapiskau):
    ipmag.plot_di(petscapiskau['dir_dec'].tolist(), petscapiskau['dir_inc'].tolist(),
                  color='darkorange', marker='^', markersize=100,
                  label='Petscapiskau (outside aureole)')
ipmag.plot_di_mean(dir_mean['dec'], dir_mean['inc'], dir_mean['alpha95'],
                   color='red', marker='s')
plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
plt.title('Baked-contact test: intrusion vs. Petscapiskau country rock')
plt.show()""")

md("""## Paleosecular variation and VGP-shape diagnostics

The anorthosite VGPs are tightly clustered (K ≫ 70), so the Deenen et al. (2011)
test flags possible under-sampling of paleosecular variation, consistent with the
prior compilation scoring R2 = 0.""")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

co("""fishqq_result = pt.fishqq_vgps(michikamau, unify_polarity=False)
fishqq_result""")

co('''try:
    svei_result = pt.svei_test_vgps(michikamau, study_lon, study_lat, model='TK03_GAD', plot=True)
except TypeError:
    svei_result = pt.svei_test_vgps(michikamau, study_lon, study_lat, model='TK03_GAD', plot=False)
    print('(SVEI elongation plot skipped: E below the TK03.GAD model minimum)')
print(f"paleolatitude = {svei_result['lat']:.1f} deg; elongation E = {svei_result['E']:.2f} "
      f"({'consistent' if svei_result['E_result'] else 'inconsistent'} with TK03.GAD)")''')

md("## The Michikamau pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Michikamau Intrusion Combined', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=1460)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | U-Pb zircon ~1460 Ma on the youngest adamellite (Krogh & Davis, 1973); magnetization age ~1450-1460 Ma. |
| 2 | Techniques and statistical analysis | **0** | AF + thermal demagnetization with stable end points; however the VGPs are over-concentrated (K ≈ 104 ≫ 70, A95 below the Deenen et al. (2011) lower bound), indicating paleosecular variation is probably under-averaged. |
| 3 | Magnetic mineralogy characterized | **1** | Titanium-free magnetite and hematite identified by thermal/AF demagnetization and polished-section study (Emslie et al., 1976; Murthy et al., 1970). |
| 4 | Field tests constrain age of magnetization | **C** | Positive baked-contact (Bruhnes) test: country rock within the aureole (site 13) carries the intrusion direction, the Petscapiskau Group outside (site 22) differs. |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous Labrador craton, north of the Grenville Front; not regionally deformed since consolidation. |
| 6 | Presence of reversals | **0** | The 12 localities are of a single (westerly) polarity. A reversed (easterly) hematite component appears within some specimens antiparallel to the magnetite, but there are no reversed site means, so no site-level reversal test applies. |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles; baked-contact test independently brackets the age. |
| | Total | **5/7** | Grade A |""")

md("""## Nordic workshop summary

The Michikamau pole is recreated at the site level from the 12 localities of
Murthy et al. (1968) + Emslie et al. (1976). The recreation reproduces the prior
compilation pole (−1.5°N/217.5°E, GPMDB 2274) and Emslie's published pole
(1.5°S/142°W, A95 4.5). The positive baked-contact test scores R4 = C.

**Note (flag for review):** the prior compilation scored R6 = 1, but the 12
localities are single (westerly) polarity — the reversal is only a within-specimen
magnetite(normal)/hematite(reversed) antiparallelism, not a site-level reversal —
so R6 is scored 0 here (R total 6 → 5; Grade A unchanged).""")

co("""michikamau_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Michikamau Intrusion Combined',
    sites=michikamau,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Primary westerly, shallow characteristic remanent magnetization carried by titanium-free magnetite and hematite; single polarity',
    tests='C+ (positive baked-contact / Bruhnes test: aureole country rock matches the intrusion, the older Petscapiskau Group outside the aureole differs)',
    gpmdb_number='2274',
    percent_reversed=0,
    demag_code=3,
    R1=1, R2=0, R3=1, R4='C', R5=1, R6=0, R7=1, Grade='A',
    nominal_age=1460, lomagage=1455, himagage=1465,
    REF_method='U-Pb zircon ~1460 Ma on the youngest Michikamau adamellite (Krogh & Davis, 1973); magnetization age ~1450-1460 Ma (Emslie, Irving & Park, 1976).',
    POLE_AUTHORS='Emslie, R. F., Irving, E., & Park, J. K.; Murthy, G. S., Fahrig, W. F., & Jones, D. L.',
    YEAR=1976,
    JOURNAL='Canadian Journal of Earth Sciences',
    VOLUME='13',
    VPAGES='1052-1057',
    TITLE='Further paleomagnetic results from the Michikamau Intrusion, Labrador',
    COMMENT='Michikamau pole recreated at the site level from the 12 localities of Murthy et al. (1968) + Emslie et al. (1976) (audited MagIC contribution 20668: citations corrected to e76-108:e68-111, Petscapiskau VGP de-antipoded to -12/254, ST-C added). Reproduces Emslie published pole (1.5S/142W, A95 4.5) and prior compilation (-1.5N/217.5E, GPMDB 2274). Positive baked-contact (Bruhnes) test -> R4=C. K~104>>70 with A95 below the Deenen bound -> R2=0. R6 set 0 (single-polarity site means; the reversal is only a within-specimen magnetite/hematite antiparallelism) -> R total 6->5 vs prior compilation; Grade A. The older Petscapiskau Group (1559 Ma, N=1) is the baked-contact reference, not exported as a pole.'
)
pt.save_nordic_summary(michikamau_summary, '1460_Michikamau')
michikamau_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
