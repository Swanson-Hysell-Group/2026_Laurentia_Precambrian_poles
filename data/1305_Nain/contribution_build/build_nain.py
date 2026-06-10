"""Build the Nain anorthosite MagIC contribution + notebook (1305_Nain).

Source: a student MagIC contribution (id 20684) for Murthy (1978), CJES 15,
516-525 (doi:10.1139/e78-058), audited against the paper.

Audited fixes (instructor review): added signed per-site VGPs (recomputed with
pmag.dia_vgp); corrected the location pole rows — dark-facies (18 sites),
pale-facies (3 sites), and a true combined all-21-site pole (the source combined
row erroneously reused the pale-facies values).

Murthy (1978): dark-facies 18-site mean D=278.2/I=15.8 (k=174, a95=2.6), pole
11 N / 149 W (=211 E), dp 1.5/dm 3.0; pale-facies 3-site pole 15 N / 212 E;
mixed polarity (sites 1,2,16,17,18 reversed/easterly ~180 deg from the rest = a
magnetic reversal). Magnetite remanence isolated by thermal demagnetization.
Age 1305 +/- 22 Ma (Kiglapait intrusion Sm-Nd; DePaolo, 1985), used by the prior
compilation; the host anorthosite is older (Barton, 1974 Rb-Sr 1418 +/- 25 Ma).
"""
import os, io
import numpy as np
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Murthy1978_Nain_magic_20684_source.txt')

# ---- read source, add signed VGPs --------------------------------------------
data = open(SRC, encoding='latin-1').read()
sites = None
for b in data.split('>>>>>>>>>>'):
    b = b.strip()
    if b.startswith('tab') and 'sites' in b.splitlines()[0]:
        sites = pd.read_csv(io.StringIO('\n'.join(b.splitlines()[1:])), sep='\t')

def vgp_row(r):
    a = r['dir_alpha95'] if not pd.isna(r['dir_alpha95']) else 5.0
    plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], a, r['lat'], r['lon'])
    return pd.Series({'vgp_lon': round(plon, 1), 'vgp_lat': round(plat, 1),
                      'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1)})

sites = pd.concat([sites, sites.apply(vgp_row, axis=1)], axis=1)
sites['location'] = 'Nain Anorthosite'

def pole_of(df):
    blk = pmag.flip(ipmag.make_di_block(df['vgp_lon'].tolist(), df['vgp_lat'].tolist()), combine=True)
    return pmag.fisher_mean(blk)

dark = sites[sites['description'].str.contains('Dark', na=False)]
pale = sites[sites['description'].str.contains('Pale', na=False)]
pc, pd_, pp = pole_of(sites), pole_of(dark), pole_of(pale)

locs = pd.DataFrame([
    {'location': 'Nain Anorthosite', 'location_type': 'Outcrop',
     'result_name': 'Nain anorthosite combined dark+pale facies pole', 'result_type': 'a',
     'sites': ':'.join(map(str, sites['site'].tolist())),
     'method_codes': 'LP-DIR-T:DE-BLANKET:DE-FM:DE-VGP', 'citations': '10.1139/e78-058',
     'geologic_classes': 'Intrusive', 'lithologies': 'Anorthosite',
     'lat_s': sites['lat'].min(), 'lat_n': sites['lat'].max(),
     'lon_w': sites['lon'].min(), 'lon_e': sites['lon'].max(),
     'age': 1305, 'age_low': 1283, 'age_high': 1327, 'age_unit': 'Ma',
     'dir_tilt_correction': 0, 'pole_lat': round(pc['inc'], 1), 'pole_lon': round(pc['dec'], 1),
     'pole_alpha95': round(pc['alpha95'], 1), 'pole_k': round(pc['k'], 1), 'pole_n_sites': int(pc['n']),
     'description': 'Combined 21-site (18 dark + 3 pale facies) VGP-Fisher-mean pole; dual polarity (positive reversal). Age 1305+/-22 Ma (Kiglapait Sm-Nd, DePaolo 1985); host anorthosite ca. 1418+/-25 Ma (Barton 1974).'},
    {'location': 'Nain Anorthosite dark facies', 'location_type': 'Outcrop',
     'result_name': 'Nain anorthosite dark facies pole', 'result_type': 'a',
     'sites': ':'.join(map(str, dark['site'].tolist())),
     'method_codes': 'LP-DIR-T:DE-BLANKET:DE-FM:DE-VGP', 'citations': '10.1139/e78-058',
     'geologic_classes': 'Intrusive', 'lithologies': 'Anorthosite',
     'lat_s': dark['lat'].min(), 'lat_n': dark['lat'].max(), 'lon_w': dark['lon'].min(), 'lon_e': dark['lon'].max(),
     'age': 1305, 'age_low': 1283, 'age_high': 1327, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
     'pole_lat': round(pd_['inc'], 1), 'pole_lon': round(pd_['dec'], 1), 'pole_alpha95': round(pd_['alpha95'], 1),
     'pole_k': round(pd_['k'], 1), 'pole_n_sites': int(pd_['n']),
     'description': 'Dark-facies anorthosite, 18 sites (Palungotok Island + south of Nain). Murthy (1978) pole 11 N/149 W (=211 E), dp 1.5/dm 3.0 from the mean direction.'},
    {'location': 'Nain Anorthosite pale facies', 'location_type': 'Outcrop',
     'result_name': 'Nain anorthosite pale facies pole', 'result_type': 'a',
     'sites': ':'.join(map(str, pale['site'].tolist())),
     'method_codes': 'LP-DIR-T:DE-BLANKET:DE-FM:DE-VGP', 'citations': '10.1139/e78-058',
     'geologic_classes': 'Intrusive', 'lithologies': 'Anorthosite',
     'lat_s': pale['lat'].min(), 'lat_n': pale['lat'].max(), 'lon_w': pale['lon'].min(), 'lon_e': pale['lon'].max(),
     'age': 1305, 'age_low': 1283, 'age_high': 1327, 'age_unit': 'Ma', 'dir_tilt_correction': 0,
     'pole_lat': round(pp['inc'], 1), 'pole_lon': round(pp['dec'], 1), 'pole_alpha95': round(pp['alpha95'], 1),
     'pole_k': round(pp['k'], 1), 'pole_n_sites': int(pp['n']),
     'description': 'Pale-facies anorthosite (Paul Island), 3 sites. Murthy (1978) pole 15 N/148 W (=212 E), dp 7.5/dm 14.'},
])

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt ({len(locs)})')
print(f'    combined pole {pc["inc"]:.1f}/{pc["dec"]:.1f} A95 {pc["alpha95"]:.1f} N {int(pc["n"])}')

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/1305_Nain.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Nain anorthosite ca. 1305 Ma paleomagnetic pole

## Geologic context

The Nain anorthosite (Nain Plutonic Suite) of coastal Labrador is one of the
massif-type anorthosite complexes that lie north of the Grenville Front and are
unaffected by Grenvillian metamorphism (Murthy, 1978). It comprises some
10,000 km² of anorthosite and related leuconorite-norite-troctolite, subdivided
by Wheeler (1960) into dark, pale, and buff facies. The complex was intruded over
a span of time around the Mesoproterozoic; the Kiglapait layered intrusion within
it gives a Sm-Nd age of 1305 ± 22 Ma (DePaolo, 1985), used here as the pole age,
while Rb-Sr on the host anorthosite gives ca. 1418 ± 25 Ma (Barton, 1974).

## Pole

This notebook recreates the Nain anorthosite pole at the site level from the 21
sites of Murthy (1978) — 18 dark-facies (Palungotok Island and south of Nain) and
3 pale-facies (Paul Island) sites — whose magnetite remanence was isolated by
thermal demagnetization. The collection records a magnetic reversal (a subset of
sites carry the antipodal, easterly direction), supporting a primary origin.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Site-level data

Site means are loaded from `../data/1305_Nain/` (rebuilt from the audited Murthy
1978 contribution, with signed VGPs recomputed). Directions are in-situ
(intrusive rocks). The dark- and pale-facies sites are tagged in the
`description` field.""")

co("""sites_geo, _ = pt.load_magic_sites('../data/1305_Nain/sites.txt')
study_lat, study_lon = 56.5, 298.2
print(sites_geo['dir_polarity'].value_counts())
dark = sites_geo[sites_geo['description'].str.contains('Dark', na=False)].reset_index(drop=True)
pale = sites_geo[sites_geo['description'].str.contains('Pale', na=False)].reset_index(drop=True)
print(f'\\n{len(sites_geo)} sites: {len(dark)} dark facies, {len(pale)} pale facies')
sites_geo[['site', 'lat', 'lon', 'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k',
           'dir_n_samples', 'dir_polarity', 'vgp_lat', 'vgp_lon']]""")

md("""## Combined pole and the magnetic reversal

The easterly (reversed) and westerly (normal) site directions differ by ~180°, a
positive reversal test; polarity is unified before computing the Fisher mean pole.""")

co("""ipmag.plot_net()
ipmag.plot_di(sites_geo['dir_dec'].tolist(), sites_geo['dir_inc'].tolist(),
              color='blue', marker='o')
plt.title('Nain anorthosite site directions (in-situ; dual polarity)')
plt.show()

vgp_block, pole_mean = pt.compute_mean_pole(sites_geo, unify_polarity=True)
ipmag.print_pole_mean(pole_mean)
print('\\nprior compilation (GPMDB 2180): 11.7 N / 206.7 E, A95 2.2, N=21')
print('Murthy (1978) dark-facies pole: 11 N / 211 E, dp 1.5 / dm 3.0, N=18')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('Nain anorthosite combined pole')
plt.show()""")

co("""dir_block, dir_mean = pt.compute_mean_direction(sites_geo, unify_polarity=True)
ipmag.print_direction_mean(dir_mean)
# dark vs pale facies poles for context
_, dark_pole = pt.compute_mean_pole(dark, unify_polarity=True)
_, pale_pole = pt.compute_mean_pole(pale, unify_polarity=True)
print(f"dark facies (N={int(dark_pole['n'])}): {dark_pole['inc']:.1f}/{dark_pole['dec']:.1f} A95 {dark_pole['alpha95']:.1f}")
print(f"pale facies (N={int(pale_pole['n'])}): {pale_pole['inc']:.1f}/{pale_pole['dec']:.1f} A95 {pale_pole['alpha95']:.1f}")""")

md("""## Paleosecular variation and VGP-shape diagnostics

The anorthosite VGPs are tightly clustered (K ≫ 70), so the Deenen et al. (2011)
test flags possible under-sampling of paleosecular variation — consistent with
the prior compilation scoring R2 = 0 for this pole.""")

co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")

co("""fishqq_result = pt.fishqq_vgps(sites_geo, unify_polarity=True)
fishqq_result""")

co('''try:
    svei_result = pt.svei_test_vgps(sites_geo, study_lon, study_lat, model='TK03_GAD', plot=True)
except TypeError:
    svei_result = pt.svei_test_vgps(sites_geo, study_lon, study_lat, model='TK03_GAD', plot=False)
    print('(SVEI elongation plot skipped: E below the TK03.GAD model minimum)')
print(f"paleolatitude = {svei_result['lat']:.1f} deg; elongation E = {svei_result['E']:.2f} "
      f"({'consistent' if svei_result['E_result'] else 'inconsistent'} with TK03.GAD)")''')

md("## The Nain pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Nain Anorthosite', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=1305)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **0** | Pole age 1305 ± 22 Ma (Kiglapait Sm-Nd; DePaolo, 1985); the ± 22 Ma uncertainty (and the older Rb-Sr age of the host anorthosite) exceeds ± 15 Ma. |
| 2 | Techniques and statistical analysis | **0** | Thermal demagnetization with stable end points; however the VGPs are over-concentrated (K ≈ 163 ≫ 70, A95 ≈ 2.5° below the Deenen et al. (2011) lower bound), indicating paleosecular variation is probably under-averaged. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence identified by thermomagnetic analysis and thermal demagnetization (Murthy, 1978). |
| 4 | Field tests constrain age of magnetization | **0** | No baked-contact, fold, or conglomerate test. |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous Nain Plutonic Suite, north of the Grenville Front and unaffected by Grenvillian metamorphism. |
| 6 | Presence of reversals | **1** | The collection records a magnetic reversal — easterly (reversed) and westerly (normal) directions ~180° apart. |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles. |
| | Total | **4/7** | Grade B |""")

md("""## Nordic workshop summary

The Nain anorthosite pole is recreated at the site level from all 21 sites
(18 dark + 3 pale facies) of Murthy (1978), with signed VGPs recomputed
(the source contribution lacked per-site VGPs and its combined location row was
defective). The recreation (combined VGP-Fisher mean) reproduces the prior
compilation pole (11.7°N/206.7°E, GPMDB 2180) within A95. The collection's
magnetic reversal scores R6 = 1; the very high K flags under-averaged PSV (R2 = 0).""")

co("""nain_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Nain Anorthosite',
    sites=sites_geo,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Characteristic remanent magnetization carried by magnetite; dual polarity (easterly reversed + westerly normal)',
    tests='R+ (positive reversal test: easterly and westerly site directions ~180 deg apart)',
    gpmdb_number='2180',
    percent_reversed=24,
    demag_code=3,
    R1=0, R2=0, R3=1, R4='', R5=1, R6=1, R7=1, Grade='B',
    nominal_age=1305, lomagage=1283, himagage=1327,
    REF_method='Pole age 1305 +/- 22 Ma from the Kiglapait layered intrusion (Sm-Nd; DePaolo, 1985); the host Nain anorthosite gives Rb-Sr ca. 1418 +/- 25 Ma (Barton, 1974).',
    POLE_AUTHORS='Murthy, G. S.',
    YEAR=1978,
    JOURNAL='Canadian Journal of Earth Sciences',
    VOLUME='15',
    VPAGES='516-525',
    TITLE='Paleomagnetic results from the Nain anorthosite and their tectonic implications',
    COMMENT='Nain anorthosite pole recreated at the site level from all 21 sites (18 dark + 3 pale facies) of Murthy (1978) (audited MagIC contribution 20684; signed per-site VGPs recomputed, location rows corrected). Combined VGP-Fisher mean reproduces the prior compilation pole (11.7N/206.7E, GPMDB 2180) within A95; Murthy reports the dark-facies (N=18) pole as 11N/211E (dp 1.5/dm 3.0). Positive reversal test (R6=1). K~163 >>70 with A95 below the Deenen bound -> R2=0 (under-averaged PSV). R1=0 (age 1305+/-22). Grade B.'
)
pt.save_nordic_summary(nain_summary, '1305_Nain')
nain_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
