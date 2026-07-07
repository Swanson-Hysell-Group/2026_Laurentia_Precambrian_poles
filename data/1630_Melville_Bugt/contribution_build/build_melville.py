"""Build the Melville Bugt dyke swarm MagIC contribution (1630_Melville).

Source: a student MagIC contribution for Halls, Hamilton & Denyszyn (2011),
"The Melville Bugt Dyke Swarm of Greenland" (Dyke Swarms book ch. 27,
doi:10.1007/978-3-642-12496-9_27), audited against the paper.

This script writes the MagIC tables (sites, locations, ages) and the combined
upload file, then validates it. The **notebook is maintained separately** (edited
via the IDE / JupyterLab), so this script no longer generates it.

Audited fixes (instructor review):
- Pole reported with an internally consistent k/A95 pair: a true Fisher mean of
  the 9 dyke VGPs gives k = 31.3, A95 = 9.3 (the source mixed the paper A95 of
  8.7 with a recomputed k). The chilled-margin site "OF Margin" is the margin of
  the OF dyke and is not counted as an independent dyke in the 9-dyke pole.
  Mixed-polarity pole may be reported at the -4 antipode.
- QT dir_n_samples stays 8: Table 27.1 lists QT as N = 8 (samples providing data,
  i.e. in the mean) and NT = 9 (collected); dir_n_samples is the in-mean count
  (a prior edit that set it to the collected 9 is reverted here).
- Nominal age 1630 Ma (mean of the dated dykes; Klausen & Nilsson 2019 also name
  it the "1630 Ma MBDS"). The four Halls et al. (2011) pole dykes are dated (2σ)
  OF 1622.1 ± 3.2, MB2 1629.4 ± 0.8, QT 1632.0 ± 1.1, MB9 1635.0 ± 2.7 Ma; the
  location bracket is their 2σ envelope, 1619-1638 Ma. (All four dated dykes are
  in the 9-dyke pole, including the youngest, OF, so the pole time-averages the
  ~13 Myr span and an old-end nominal such as 1635 is not appropriate.)
- Source read as UTF-8 (was latin-1, which corrupted the ± in the MB2 date).
- Ages table added for the four U-Pb-dated dykes (MB9, QT, MB2, OF). The paper
  reports 2sigma uncertainties; age_sigma is the 1sigma (half of that).
"""
import io
import os
import sys
from datetime import date

import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Halls2011_MelvilleBugt_magic_source.txt')
DOI = '10.1007/978-3-642-12496-9_27'

# read UTF-8 (reading as latin-1 corrupted the ± in the MB2 U-Pb date)
data = open(SRC, encoding='utf-8').read()
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

# undated dykes carry the swarm's age bracket = the 2sigma envelope of the four
# dated dykes (1619-1638 Ma); the dated dykes (OF, MB2, QT, MB9) keep their own
# point ages (age_low == age_high)
_bracket = (sites['age_low'] == 1622) & (sites['age_high'] == 1635)
sites.loc[_bracket, 'age_low'] = 1619
sites.loc[_bracket, 'age_high'] = 1638

# 9-dyke pole (primary component, VGP present, exclude the OF chilled margin)
nine = sites[(sites['dir_comp_name'] == 'primary') & (sites['vgp_lat'].notna()) &
             (sites['site'] != 'OF Margin')]
blk = pmag.flip(ipmag.make_di_block(nine['vgp_lon'].tolist(),
                                    nine['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)
locs['pole_lat'] = round(p['inc'], 1)
locs['pole_lon'] = round(p['dec'], 1)
locs['pole_alpha95'] = round(p['alpha95'], 1)
locs['pole_k'] = round(p['k'], 1)
locs['pole_n_sites'] = int(p['n'])

# rich location metadata: the polarity-unified mean direction (underlies the
# pole) and geographic descriptors. No reversal_test column: the swarm is of
# mixed polarity, but the reversal test is equivocal (fails the Watson V common-
# mean test though the bootstrap passes), so it is not recorded as a positive
# field test. (Dual polarity is noted in the description; the equivocal reversal
# test means R6 = 0 in the notebook.)
dblk = pmag.flip(ipmag.make_di_block(nine['dir_dec'].tolist(),
                                     nine['dir_inc'].tolist()), combine=True)
dm = pmag.fisher_mean(dblk)
locs['dir_dec'] = round(dm['dec'], 1)
locs['dir_inc'] = round(dm['inc'], 1)
locs['dir_alpha95'] = round(dm['alpha95'], 1)
locs['dir_k'] = round(dm['k'], 1)
locs['dir_n_sites'] = int(dm['n'])
# dir_n_samples is intentionally not set at the location level: the sample count
# is inherent in the sites table (sum of per-site dir_n_samples)
locs['continent_ocean'] = 'Greenland'
locs['country'] = 'Greenland'

# adopted nominal 1630 Ma (mean of the four dated pole dykes) with the 2sigma
# envelope of those dates as the bracket (OF 1622.1-3.2 to MB9 1635.0+2.7)
locs['age'] = 1630
locs['age_low'] = 1619
locs['age_high'] = 1638

locs['result_name'] = 'Melville Bugt dyke swarm ca. 1630 Ma pole'
desc = str(locs['description'].iloc[0]).replace(
    # drop a superfluous source sentence (k is already given below with A95)
    ' pole_k of 31.35 from Fisher recalculation of site VGPs.', '').replace(
    'nominal pole age ~1628 Ma (midpoint of dated sites).',
    'nominal pole age 1630 Ma (mean of the four dated pole dykes: OF 1622.1 +/- '
    '3.2, MB2 1629.4 +/- 0.8, QT 1632.0 +/- 1.1, MB9 1635.0 +/- 2.7 Ma, all 2sigma; '
    'Halls et al., 2011); the location age bracket 1619-1638 Ma is the 2sigma '
    'envelope of those dates. A fifth dyke ~1000 km SE (Snehatten) is 1629.7 +/- '
    '3.9 Ma (Klausen & Nilsson, 2019).')
locs['description'] = (
    desc +
    ' Pole is the Fisher mean of the 9 dyke VGPs (k=%.1f, A95=%.1f); the OF '
    'chilled margin is the margin of the OF dyke and is not counted separately.'
    % (p['k'], p['alpha95']))

# ---- ages table: the four U-Pb-dated dykes ---------------------------------
# The paper reports 2sigma (95% conf.); age_sigma is the 1sigma (half of that).
AGE_COLS = ['location', 'site', 'age', 'age_sigma', 'age_unit', 'method_codes',
            'citations', 'timescale_eon', 'timescale_era', 'timescale_period',
            'description']
_DATED = [  # site, age, 2sigma, published-uncertainty note
    ('MB9', 1635.0, 2.7),
    ('QT',  1632.0, 1.1),
    ('MB2', 1629.4, 0.8),
    ('OF',  1622.1, 3.2),
]
ages = pd.DataFrame([{
    'location': 'Melville Bugt', 'site': s, 'age': a,
    'age_sigma': round(two_sigma / 2, 2), 'age_unit': 'Ma',
    'method_codes': 'GM-UPB', 'citations': DOI,
    'timescale_eon': 'Proterozoic', 'timescale_era': 'Paleoproterozoic',
    'timescale_period': 'Statherian',
    'description': 'U-Pb baddeleyite %.1f ± %.1f Ma (2σ; Halls et al., 2011)'
                   % (a, two_sigma),
} for s, a, two_sigma in _DATED], columns=AGE_COLS)


def write_magic(df, kind, path, cols=None):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', columns=cols,
              encoding='utf-8')


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
write_magic(ages, 'ages', os.path.join(OUT, 'ages.txt'), cols=AGE_COLS)
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt, ages.txt ({len(ages)}); '
      f'pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} K {p["k"]:.1f} '
      f'N {int(p["n"])}')

# ---- combined upload file (locations + sites + ages) -----------------------
stamp = date.today().strftime('%d.%b.%Y')
out_upload = os.path.join(OUT, f'Halls2011_MelvilleBugt_{stamp}.txt')
with open(out_upload, 'w', encoding='utf-8') as f:
    for kind, df, cols in [('locations', locs, None), ('sites', sites, None),
                           ('ages', ages, AGE_COLS)]:
        if kind != 'locations':
            f.write('>>>>>>>>>>\n')
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False, columns=cols)
print(f'-I- wrote {out_upload} (locations + sites + ages, canonical markers)')

# ---- validate before submission --------------------------------------------
repo_root = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(repo_root, 'scripts'))
try:
    from validate_magic_contribution import validate_upload_file
    if not validate_upload_file(out_upload, tables=['locations', 'sites', 'ages']):
        print('-W- validation reported problems above; review before uploading')
except Exception as exc:  # pragma: no cover
    print(f'-W- could not run validator ({exc}); validate on upload instead')
