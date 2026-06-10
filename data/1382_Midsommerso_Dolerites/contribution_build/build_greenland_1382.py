"""Build the Midsommerso Dolerites + Zig-Zag Dal Basalt notebooks (both 1382).

Source: a student MagIC contribution for Marcussen & Abrahamsen (1983),
"Palaeomagnetism of the Proterozoic Zig-Zag Dal Basalt and the Midsommerso
Dolerites, eastern North Greenland" (doi:10.1111/j.1365-246X.1983.tb03328.x),
audited against the paper. One contribution carries both units; two notebooks are
built (1382_Midsommersoe and 1382_Zigzag).

Audited fixes (instructor review): corrupted Upton et al. (2005) DOI corrected
(...0634-7); dir_tilt_correction set to 100 on all site rows (Marcussen &
Abrahamsen report bedding-corrected directions; the revision had wrongly set 0,
inconsistent with the location rows).

Both units are ca. 1382 Ma (U-Pb baddeleyite on the intrusive Midsommerso
Dolerites, 1382 +/- 2 Ma; Upton et al., 2005). Prior compilation poles
(VGP-Fisher mean convention here, reported at the positive-latitude antipode):
Midsommerso Dolerites GPMDB 99 (10.0 N / 242.0 E, N=10); Zig-Zag Dal Basalts
GPMDB 98 (12.2 N / 242.8 E, N=19). Both single polarity; over-concentrated VGPs
(K > 70) -> R2 = 0; no field test; R = 4, Grade B.
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
MID_OUT = os.path.dirname(HERE)
ZZ_OUT = os.path.abspath(os.path.join(HERE, '../../1382_ZigZag_Dal_Basalt'))
SRC = os.path.join(HERE, 'Marcussen1983_Greenland_magic_source.txt')

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

sites['dir_tilt_correction'] = 100   # bedding-corrected directions (issues fix)

def write_magic(df, kind, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

# split by formation, write a contribution into each data dir
for loc, out in [('Midsommerso Dolerites', MID_OUT), ('Zig-Zag Dal Basalt Formation', ZZ_OUT)]:
    s = sites[sites['location'] == loc].copy()
    lo = locs[locs['location'] == loc].copy()
    write_magic(s, 'sites', os.path.join(out, 'sites.txt'))
    write_magic(lo, 'locations', os.path.join(out, 'locations.txt'))

print('-I- wrote sites/locations for both units')


def build_notebook(title_md, rockname, gpmdb, datadir, nb_name, pole_target, dir_target,
                   N, slat, slon, context_md):
    NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/' + nb_name))
    cells = []
    md = lambda s: cells.append(new_markdown_cell(s))
    co = lambda s: cells.append(new_code_cell(s))
    md(title_md)
    co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")
    md(context_md)
    co(f"""sites_geo, sites_tc = pt.load_magic_sites('../data/{datadir}/sites.txt')
sub = sites_tc if len(sites_tc) else sites_geo
study_lat, study_lon = {slat}, {slon}
print(f'{{len(sub)}} sites (bedding-corrected)')
sub[['site', 'lat', 'lon', 'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k', 'dir_n_samples', 'vgp_lat', 'vgp_lon']]""")
    co(f"""# the VGPs are southern in the student convention; report at the positive-latitude
# antipode to match the prior compilation (GPMDB {gpmdb})
vgp_block, pm = pt.compute_mean_pole(sub, unify_polarity=True)
if pm['inc'] < 0:
    vgp_block, pm = pt.compute_mean_pole(sub, unify_polarity=True, flip=True)
pole_mean = pm
ipmag.print_pole_mean(pole_mean)
print('\\nprior compilation (GPMDB {gpmdb}): {pole_target}, N={N}')
_ = pt.plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'], figsize=(6, 6))
plt.title('{rockname} pole')
plt.show()""")
    co("""dir_block, dir_mean = pt.compute_mean_direction(sub, unify_polarity=True)
ipmag.print_direction_mean(dir_mean)
ipmag.plot_net()
ipmag.plot_di(di_block=dir_block, color='blue', marker='o')
ipmag.plot_di_mean(dir_mean['dec'], dir_mean['inc'], dir_mean['alpha95'], color='red', marker='s')
plt.title('Site directions (bedding-corrected)')
plt.show()""")
    md("""## Paleosecular variation

The VGPs are over-concentrated (K > 70), so the Deenen et al. (2011) test flags
under-sampling of paleosecular variation — consistent with the prior compilation
scoring R2 = 0.""")
    co("""pt.Deenen_test(pole_mean['n'], pole_mean['alpha95'])
pt.plot_Deenen_test(pole_mean)""")
    co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, pole_mean['inc'], pole_mean['dec'],
                          pole_mean['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=pole_mean['dec'],
                          central_latitude=pole_mean['inc'])
plt.show()""")
    co(f"""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('{rockname}', Laurentia_stricto_poles, pt.Torsvik2012_Laurentia,
                     pole_plat=pole_mean['inc'], pole_plon=pole_mean['dec'],
                     pole_A95=pole_mean['alpha95'], pole_age=1382)""")
    md(f"""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | U-Pb baddeleyite 1382 ± 2 Ma on the intrusive Midsommerso Dolerites (Upton et al., 2005). |
| 2 | Techniques and statistical analysis | **0** | AF demagnetization; the VGPs are over-concentrated (K > 70, A95 below the Deenen et al. (2011) lower bound), indicating under-averaged paleosecular variation. |
| 3 | Magnetic mineralogy characterized | **1** | Magnetite remanence (Marcussen & Abrahamsen, 1983). |
| 4 | Field tests constrain age of magnetization | **0** | No baked-contact or fold test. |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous eastern North Greenland (Laurentia-Greenland); bedding-corrected. |
| 6 | Presence of reversals | **0** | Single polarity. |
| 7 | No resemblance to younger poles | **1** | Distinct from younger Laurentia poles. |
| | Total | **4/7** | Grade B |""")
    co(f"""summary = pt.make_nordic_summary(
    terrane='Laurentia-Greenland',
    rockname='{rockname}',
    sites=sub,
    dir_mean=dir_mean,
    pole_mean=pole_mean,
    study_lon=study_lon,
    study_lat=study_lat,
    component_comment='Characteristic remanent magnetization (magnetite); single polarity',
    tests='',
    gpmdb_number='{gpmdb}',
    percent_reversed=0,
    demag_code=3,
    R1=1, R2=0, R3=1, R4='', R5=1, R6=0, R7=1, Grade='B',
    nominal_age=1382, lomagage=1380, himagage=1384,
    REF_method='U-Pb baddeleyite age 1382 +/- 2 Ma on the intrusive Midsommerso Dolerites (Upton et al., 2005); the Zig-Zag Dal Basalts are correlated.',
    POLE_AUTHORS='Marcussen, C., & Abrahamsen, N.',
    YEAR=1983,
    JOURNAL='Geophysical Journal of the Royal Astronomical Society',
    VOLUME='73',
    VPAGES='',
    TITLE='Palaeomagnetism of the Proterozoic Zig-Zag Dal Basalt and the Midsommerso Dolerites, eastern North Greenland',
    COMMENT='{rockname} pole recreated at the site level from Marcussen & Abrahamsen (1983) (audited contribution: dir_tilt_correction set to 100 = bedding-corrected; corrupted Upton 2005 DOI fixed to ...0634-7). Reported at the positive-latitude antipode to match the prior compilation (GPMDB {gpmdb}; the student computed the southern antipode). Single polarity, over-concentrated VGPs (K>70) -> R2=0; no field test; U-Pb 1382+/-2 Ma -> R1=1. R=4, Grade B.'
)
pt.save_nordic_summary(summary, '{datadir}')
summary""")
    nb = new_notebook(cells=cells, metadata={
        'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
        'language_info': {'name': 'python'}})
    nbformat.write(nb, NB)
    print('wrote', NB, len(cells), 'cells')


build_notebook(
    title_md="""# Midsommerso Dolerites ca. 1382 Ma paleomagnetic pole

## Geologic context

The Midsommerso Dolerites are mafic sills and dykes intruding the Mesoproterozoic
sedimentary succession and the Zig-Zag Dal Basalt of eastern North Greenland
(Marcussen & Abrahamsen, 1983). U-Pb baddeleyite dating gives 1382 ± 2 Ma (Upton
et al., 2005). The remanence is a single-polarity, magnetite-carried
magnetization.

## Pole

This notebook recreates the Midsommerso Dolerites pole at the site level from the
10 sites of Marcussen & Abrahamsen (1983).""",
    rockname='Midsommerso Dolerites', gpmdb='99',
    datadir='1382_Midsommerso_Dolerites', nb_name='1382_Midsommersoe.ipynb',
    pole_target='10.0 N / 242.0 E', dir_target='85.6/-14.5', N=10, slat=81.6, slon=333.4,
    context_md="""## Site-level data

Site means are loaded from `../data/1382_Midsommerso_Dolerites/` (the bedding-
corrected directions of Marcussen & Abrahamsen, 1983).""")

build_notebook(
    title_md="""# Zig-Zag Dal Basalt Formation ca. 1382 Ma paleomagnetic pole

## Geologic context

The Zig-Zag Dal Basalt Formation is a thick Mesoproterozoic continental flood-
basalt succession of eastern North Greenland, intruded by the comagmatic
Midsommerso Dolerites (Marcussen & Abrahamsen, 1983). U-Pb baddeleyite on the
dolerites gives 1382 ± 2 Ma (Upton et al., 2005), dating the magmatic event. The
remanence is single-polarity, magnetite-carried.

## Pole

This notebook recreates the Zig-Zag Dal Basalt pole at the site level from the
19 sites of Marcussen & Abrahamsen (1983).""",
    rockname='Zig-Zag Dal Basalts', gpmdb='98',
    datadir='1382_ZigZag_Dal_Basalt', nb_name='1382_Zigzag.ipynb',
    pole_target='12.2 N / 242.8 E', dir_target='93.6/-22', N=19, slat=81.2, slon=334.8,
    context_md="""## Site-level data

Site means are loaded from `../data/1382_ZigZag_Dal_Basalt/` (the bedding-
corrected directions of Marcussen & Abrahamsen, 1983).""")
