"""Build the Michikamau Intrusion MagIC contribution (1469_Michikamau).

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
Petscapiskau outside differs). Age 1469 +/- 1 Ma (U-Pb baddeleyite; Kerr &
McNicoll, 2010), superseding the older ~1460 Ma U-Pb zircon of Krogh & Davis
(1973).
"""
import os, io
from datetime import date
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Emslie1976_Murthy1968_Michikamau_magic_20668_source.txt')

# read UTF-8 (the source uses a real degree glyph; latin-1 corrupted it to 'Âº')
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

# --- ages: intrusion U-Pb baddeleyite 1469 +/- 1 Ma (Kerr & McNicoll, 2010;
#     ID-TIMS concordia upper-intercept on pegmatitic leuconorite, Juno zone),
#     recorded as a 1468-1470 Ma bracket. This supersedes the older U-Pb zircon
#     ~1460 Ma of Krogh & Davis (1973) and matches the notebook/Nordic/compilation
#     nominal (1469). The Petscapiskau Group baked-contact reference is the older
#     host rock (Rb-Sr 1559 +/- 60 Ma; Emslie et al. 1976), carried as 1499-1619
#     on the Petscapiskau rows only. ---
for df in (sites, locs):
    if 'age_low' not in df.columns:
        pos = df.columns.get_loc('age_unit')
        df.insert(pos, 'age_low', pd.NA)
        df.insert(pos + 1, 'age_high', pd.NA)

pet_s = sites['location'].str.contains('Petscapiskau', na=False)
sites.loc[~pet_s, ['age', 'age_low', 'age_high']] = [1469, 1468, 1470]
sites.loc[pet_s, ['age', 'age_low', 'age_high']] = [1559, 1499, 1619]
if 'age_sigma' in sites.columns:
    sites['age_sigma'] = ''

pet_l = locs['location'].str.contains('Petscapiskau', na=False)
locs.loc[mich, ['age', 'age_low', 'age_high']] = [1469, 1468, 1470]
locs.loc[pet_l, ['age', 'age_low', 'age_high']] = [1559, 1499, 1619]
if 'age_sigma' in locs.columns:
    locs['age_sigma'] = ''

# structured positive baked-contact test (normal; Emslie 1976 Bruhnes test) on the
# Michikamau pole location only
locs['contact_test'] = ''
locs.loc[mich, 'contact_test'] = 'C+'
locs.loc[mich, 'pole_reversed_perc'] = 0   # single (westerly) polarity pole

# --- lithology audit vs the papers: site 7 is Murthy et al. (1968) site 1, the
#     "olivine gabbro of the border group" (not anorthosite as the source had);
#     sites 1-2 (Emslie 18/19) are also olivine gabbro, the rest anorthosite. ---
sites.loc[sites['site'].astype(str) == '7', 'lithologies'] = 'Olivine Gabbro'
# specific geologic_classes: the intrusion (anorthosite/gabbro) is Intrusive;
# the baked biotite gneiss (site 6) and Petscapiskau meta-andesite (site 13)
# stay Metamorphic
sites.loc[sites['geologic_classes'] == 'Igneous', 'geologic_classes'] = 'Intrusive'

# wording: "No tilt correction needed" -> "No tilt correction applied"
sites['description'] = sites['description'].str.replace(
    'No tilt correction needed', 'No tilt correction applied', regex=False)

# --- richer location metadata (geographic / tectonic) on both locations ---
for col, val in [('continent_ocean', 'North America'), ('country', 'Canada'),
                 ('state_province', 'Newfoundland and Labrador'),
                 ('terranes', 'Laurentia')]:
    locs[col] = val

# --- evergreen location descriptions (no pole-result restatement, no project
#     workflow references) ---
locs.loc[mich, 'description'] = (
    'Paleomagnetic pole for the Michikamau Intrusion, a large layered '
    'leucotroctolite-anorthosite massif with minor late adamellite in central '
    'Labrador, combining the results of Murthy et al. (1968) and Emslie et al. '
    '(1976). Directions are in geographic (in-situ) coordinates: Murthy et al. '
    '(1968) report a negative fold test using the igneous layering, so no tilt '
    'correction is applied. The westerly, shallow characteristic magnetization is '
    'carried by titanium-free magnetite and hematite. A positive baked-contact '
    'test (Emslie et al., 1976) supports a primary, emplacement-age magnetization: '
    'baked biotite gneiss within the thermal aureole carries the intrusion '
    'direction, whereas the older Petscapiskau Group meta-andesite outside the '
    'aureole does not. Age from U-Pb baddeleyite dating of Michikamau anorthosite '
    '(1469 +/- 1 Ma; Kerr and McNicoll, 2010), superseding the U-Pb zircon ~1460 '
    'Ma of Krogh and Davis (1973). Original Murthy/Emslie site numbers are given '
    'in the site descriptions.')
locs.loc[pet_l, 'description'] = (
    'Petscapiskau Group meta-andesite (Emslie et al., 1976 site 22), the '
    'unbaked country rock outside the Michikamau thermal aureole (Rb-Sr ~1559 '
    'Ma). Its distinct northeast direction is the reference arm of the positive '
    'baked-contact test for the Michikamau Intrusion and is not a paleopole in '
    'its own right.')

# The Petscapiskau Group location is the baked-contact reference (a single site's
# direction), NOT a paleopole: clear the pole-result columns, mark it an
# individual (result_type 'i') result, and rename it away from "... pole".
_pole_cols = ['pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k', 'pole_n_sites']
for _pc in _pole_cols:
    locs[_pc] = locs[_pc].astype(object)   # allow blanks alongside the Michikamau floats
locs.loc[pet_l, _pole_cols] = ''
locs.loc[pet_l, 'result_type'] = 'i'
locs.loc[pet_l, 'result_name'] = 'Petscapiskau Group baked-contact reference'


def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


def write_combined(path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('tab delimited\tlocations\n'); locs.to_csv(f, sep='\t', index=False)
        f.write('>>>>>>>>>>\n')
        f.write('tab delimited\tsites\n'); sites.to_csv(f, sep='\t', index=False)


def validate(path):
    import sys
    root = os.path.dirname(os.path.dirname(OUT))
    sys.path.insert(0, os.path.join(root, 'scripts'))
    from validate_magic_contribution import validate_upload_file
    return validate_upload_file(path, tables=['locations', 'sites'])


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Michikamau-Intrusion_{stamp}.txt')
write_combined(combined)
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt, combined '
      f'{os.path.basename(combined)}; Michikamau pole '
      f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])}')
validate(combined)
