"""Build the Nain anorthosite MagIC contribution (1305_Nain).

Source: a student MagIC contribution (id 20684) for Murthy (1978), Can. J.
Earth Sci. 15, 516-525 (doi:10.1139/e78-058), audited against the paper.

Audited fixes:
  - signed per-site VGPs recomputed with pmag.dia_vgp (source lacked them);
  - location rows rebuilt as VGP-Fisher-mean poles for the combined 21-site set
    plus the dark (18) and pale (3) facies (the source combined row erroneously
    reused the pale-facies values);
  - age made consistent with the notebook / Nordic export: 1305 Ma, bracket
    1283-1328 (the source carried the stale 1320 +/- 30 label);
  - rich location metadata added (continent_ocean / country / state_province /
    terranes / pole_reversed_perc); descriptions rewritten as evergreen data.

Directions verified against Murthy (1978): thermal demagnetization to stable end
points (AF was tried but not useful); no PCA; remanence carried by magnetite with
subordinate hematite (thermomagnetic curves). Five sites (1, 2, 16, 17, 18, all
dark facies) carry the antipodal easterly direction, recording a magnetic
reversal; no formal (statistical) reversal test or baked-contact/fold/
conglomerate test was performed, so no structured field-test CV column is set.
Age: the pole age of 1305 +/- 22 Ma is the Sm-Nd age of the Kiglapait layered
intrusion within the complex (DePaolo, 1985); the host anorthosite gives an older
Rb-Sr age of ca. 1418 +/- 25 Ma (Barton, 1974).
"""
import os, io
from datetime import date
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Murthy1978_Nain_magic_20684_source.txt')

AGE, AGE_LOW, AGE_HIGH = 1305, 1283, 1328   # matches notebook / Nordic export

# ---- read source (UTF-8), add signed VGPs ------------------------------------
data = open(SRC, encoding='utf-8').read()
sites = None
for b in data.split('>>>>>>>>>>'):
    b = b.strip()
    if b.startswith('tab') and 'sites' in b.splitlines()[0]:
        sites = pd.read_csv(io.StringIO('\n'.join(b.splitlines()[1:])), sep='\t')


def vgp_row(r):
    a = r['dir_alpha95'] if not pd.isna(r['dir_alpha95']) else 5.0
    plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], a, r['lat'], r['lon'])
    return pd.Series({'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
                      'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1)})


sites = pd.concat([sites, sites.apply(vgp_row, axis=1)], axis=1)
sites['location'] = 'Nain Anorthosite'

# age consistent across the whole contribution (notebook / Nordic / MagIC)
sites['age'] = AGE
sites['age_low'] = AGE_LOW
sites['age_high'] = AGE_HIGH
sites = sites.drop(columns=[c for c in ['age_sigma'] if c in sites.columns])


def pole_of(df):
    blk = pmag.flip(ipmag.make_di_block(df['vgp_lon'].tolist(), df['vgp_lat'].tolist()), combine=True)
    return pmag.fisher_mean(blk)


def pct_reversed(df):
    return round(100.0 * (df['dir_polarity'] == 'r').sum() / len(df))


dark = sites[sites['description'].str.contains('Dark', na=False)]
pale = sites[sites['description'].str.contains('Pale', na=False)]
pc, pd_, pp = pole_of(sites), pole_of(dark), pole_of(pale)

GEO = dict(continent_ocean='North America', country='Canada',
           state_province='Newfoundland and Labrador', terranes='Laurentia')


def loc_row(name, result_name, subset, pole, description):
    r = {'location': name, 'location_type': 'Region', 'result_name': result_name,
         'result_type': 'a', 'result_quality': 'g',
         'sites': ':'.join(map(str, subset['site'].tolist())),
         'method_codes': 'LP-DIR-T:DE-BLANKET:DE-FM:DE-VGP', 'citations': '10.1139/e78-058',
         'geologic_classes': 'Intrusive', 'lithologies': 'Anorthosite',
         'lat_s': subset['lat'].min(), 'lat_n': subset['lat'].max(),
         'lon_w': subset['lon'].min(), 'lon_e': subset['lon'].max(),
         'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
         'dir_tilt_correction': 0,
         'pole_lat': round(pole['inc'], 1), 'pole_lon': round(pole['dec'], 1),
         'pole_alpha95': round(pole['alpha95'], 1), 'pole_k': round(pole['k'], 1),
         'pole_n_sites': int(pole['n']), 'pole_reversed_perc': pct_reversed(subset),
         'description': description}
    r.update(GEO)
    return r


locs = pd.DataFrame([
    loc_row('Nain Anorthosite', 'Nain anorthosite ca. 1305 Ma pole', sites, pc,
            "Paleomagnetic pole for the Nain anorthosite (Nain Plutonic Suite), a "
            "massif-type anorthosite complex of coastal Labrador that lies north of "
            "the Grenville Front and is unaffected by Grenvillian metamorphism "
            "(Murthy, 1978). The pole is the Fisher mean of the in-situ site virtual "
            "geomagnetic poles from all 21 sites, combining the dark facies "
            "(olivine-bearing anorthosite; Palungotok Island and the region south of "
            "Nain) and the pale facies (hypersthene-bearing anorthosite; Paul "
            "Island). The characteristic remanence, carried by magnetite with "
            "subordinate hematite, was isolated by thermal demagnetization to stable "
            "end points; five sites carry the antipodal easterly (reversed) "
            "direction, recording a magnetic reversal consistent with a primary "
            "origin. The pole age is the Sm-Nd age of the Kiglapait layered intrusion "
            "within the complex (DePaolo, 1985); the host anorthosite yields an older "
            "Rb-Sr age (Barton, 1974). Original Murthy (1978) site numbers are "
            "retained in the site descriptions."),
    loc_row('Nain Anorthosite dark facies', 'Nain anorthosite dark-facies pole', dark, pd_,
            "Dark-facies anorthosite (olivine-bearing) of the Nain complex, sampled "
            "on Palungotok Island and in the region south of Nain (Murthy, 1978, "
            "sites 1-18). The pole is the Fisher mean of the in-situ site virtual "
            "geomagnetic poles; both normal (westerly) and reversed (easterly) "
            "polarities are present."),
    loc_row('Nain Anorthosite pale facies', 'Nain anorthosite pale-facies pole', pale, pp,
            "Pale-facies anorthosite (hypersthene-bearing) of the Nain complex, "
            "sampled on Paul Island (Murthy, 1978, sites 19-21). The pole is the "
            "Fisher mean of the in-situ site virtual geomagnetic poles. Murthy (1978) "
            "found the pale-facies mean direction to differ significantly from the "
            "dark facies."),
])


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
combined = os.path.join(OUT, f'Murthy1978_Nain_{stamp}.txt')
write_combined(combined)
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt ({len(locs)}), combined '
      f'{os.path.basename(combined)}')
print(f'    combined pole {pc["inc"]:.1f}/{pc["dec"]:.1f} A95 {pc["alpha95"]:.1f} '
      f'K {pc["k"]:.1f} N {int(pc["n"])}')
print(f'    dark  {pd_["inc"]:.1f}/{pd_["dec"]:.1f} A95 {pd_["alpha95"]:.1f} N {int(pd_["n"])}')
print(f'    pale  {pp["inc"]:.1f}/{pp["dec"]:.1f} A95 {pp["alpha95"]:.1f} N {int(pp["n"])}')
validate(combined)
