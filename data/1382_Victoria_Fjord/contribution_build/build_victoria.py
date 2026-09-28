"""Build the Victoria Fjord dolerite dyke pole contribution (1382_Victoria).

Source: Abrahamsen & Van der Voo (1987), "Palaeomagnetism of middle Proterozoic
(c. 1.25 Ga) dykes from central North Greenland," Geophys. J. R. astr. Soc. 91,
597-611 (doi:10.1111/j.1365-246X.1987.tb01660.x). The site table is transcribed
directly from the paper's Table 1 (p. 601, read from page images) below; an
earlier student transcription (Abrahamsen1987_Victoria_sites_source.txt) is kept
for reference only.

Audit against the paper (2026-09-28):
  - D6 (site 269) D = 267 (Table 1); a prior build carried 275.
  - D8 (site 272) N1 = 9 (Table 1: N 9, N1 9, k 56, R 8.858 -> (N1-1)/(N1-R) = 56);
    the student source had 5. The dolerite mean then uses 73 samples (Table 2).
  - Treatment per site from Table 1: D4, D5 and the three baked sites are AF only;
    the other dykes AF + thermal. Directions are stable end points read from
    Zijderveld plots (no PCA) -> LP-DIR-AF[:LP-DIR-T]:DE-FM.
  - Orientation: solar compass for nine dyke sites, Brunton magnetic compass for
    D10 because of overcast conditions (p. 600) -> SO-SUN / SO-MAG.
  - D10 (site 279) gave an anomalous NE direction, 'probably of a different
    (younger?) age' (p. 601), and is excluded from the pole (result_quality 'b').
  - The baked-contact test is on the 23 m wide dyke D1 (site 260; Fig. 5): site 261
    in gneiss 0.1-0.3 m from the south contact, site 262 single cores in gneiss at
    1, 2, 4, 8, 16 and 32 m, site 268 in the small quartz-dioritic body within 3 m
    of the north contact; an unbaked gneiss site (263) 60 m S has no tabulated
    statistics. Resetting is complete to 16 m; the 32 m core is scattered and the
    60 m site shows no resetting (Fig. 6). The three baked sites are included as
    rows but excluded from the pole ('b'), as they record the cooling of D1.
  - Paper p. 600: '10 dolerite sites from eight dykes' -- which sites share a
    dyke is not stated.
  - Directions are in situ (near-vertical dykes in basement unconformably overlain
    by flat-lying sediments; no tilt correction in the paper).
  - Precambrian single-polarity dykes: no dir_polarity column; pole_reversed_perc 0
    (no sites of opposite polarity).

Age: the dykes are undated. The compilation assigns ca. 1382 Ma by correlation
(petrographic similarity and antiparallel directions) with the Midsommerso
Dolerites, one sill of which gave a U-Pb baddeleyite age of 1382 +/- 2 Ma
(Upton et al., 2005, doi:10.1007/s00410-004-0634-7). The paper itself infers
c. 1.25 Ga from Rb-Sr isochrons on the correlative eastern North Greenland
intrusives and notes that the opposite polarities may indicate a slight
difference in intrusive age (p. 609).
"""
import os
from datetime import date
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)

CIT = '10.1111/j.1365-246X.1987.tb01660.x'
CIT_AGE = '10.1111/j.1365-246X.1987.tb01660.x:10.1007/s00410-004-0634-7'
AGE, AGE_LOW, AGE_HIGH = 1382, 1380, 1384
LAT, LON = 81.5, 315.3            # 81.5N, 44.7W (p. 600); one coordinate for all sites
LOCATION = 'Victoria Fjord dolerite dykes'

# Table 1 (p. 601): site, paper site no., N (collected), N1 (measured/used),
# treatment, D, I, k, alpha95, R
TABLE1 = [
    ('D1', '260', 10, 10, 'AF+TH', 265, 21, 41, 7.6, 9.782),
    ('D2', '264', 10, 10, 'AF+TH', 264, 29, 38, 7.9, 9.763),
    ('D3', '265', 10, 10, 'AF+TH', 273, 20, 121, 4.4, 9.926),
    ('D4', '266', 8, 6, 'AF', 268, 21, 150, 5.5, 5.967),
    ('D5', '267', 6, 6, 'AF', 274, 14, 1807, 1.6, 5.9972),
    ('D6', '269', 8, 8, 'AF+TH', 267, 28, 149, 4.6, 7.953),
    ('D7', '271', 10, 9, 'AF+TH', 258, 29, 72, 6.1, 8.889),
    ('D8', '272', 9, 9, 'AF+TH', 260, 37, 56, 6.9, 8.858),
    ('D9', '273', 5, 5, 'AF+TH', 256, 27, 152, 6.2, 4.974),
    ('D10', '279', 9, 8, 'AF+TH', 46, 30, 10, 18.5, 7.294),
    ('261', '', 10, 9, 'AF', 265, 12, 92, 5.4, 8.913),
    ('262', '', 6, 5, 'AF', 266, 4, 164, 6.0, 4.976),
    ('268', '', 6, 5, 'AF', 264, 14, 148, 6.3, 4.973),
]

AGE_NOTE = ('Undated; age assigned by correlation with the U-Pb baddeleyite-dated '
            'Midsommerso Dolerites (Upton et al., 2005).')
SITE_INFO = {
    'D1': ('dyke', 'Dolerite dyke D1 (paper site 260), c. 23 m wide; the dyke of the '
                   'baked-contact test. ' + AGE_NOTE),
    'D10': ('dyke', 'Dolerite dyke D10 (paper site 279). Oriented by magnetic compass '
                    'under overcast conditions; anomalous NE direction, probably of a '
                    'different age, excluded from the mean by the authors. [excluded]'),
    '261': ('gneiss', 'Baked Archaean gneiss 0.1-0.3 m from the south contact of dyke D1. '
                      '[excluded: records the cooling of dyke D1]'),
    '262': ('gneiss', 'Baked Archaean gneiss, single cores at 1, 2, 4, 8, 16 and 32 m from '
                      'the south contact of dyke D1. [excluded: records the cooling of '
                      'dyke D1]'),
    '268': ('qdiorite', 'Baked quartz-diorite (small quartz-dioritic body in the gneiss) '
                        'within 3 m of the north contact of dyke D1. [excluded: records '
                        'the cooling of dyke D1]'),
}


def site_row(t):
    site, alt, n, n1, treat, dec, inc, k, a95, r = t
    kind, desc = SITE_INFO.get(site, ('dyke', None))
    if desc is None:
        desc = f'Dolerite dyke {site} (paper site {alt}). ' + AGE_NOTE
    codes = ['LP-DIR-AF'] + (['LP-DIR-T'] if 'TH' in treat else []) + ['DE-FM']
    if kind == 'dyke':
        codes = (['SO-MAG'] if site == 'D10' else ['SO-SUN']) + codes
    geo = {'dyke': ('Intrusive', 'Volcanic Dike', 'Diabase'),
           'gneiss': ('Metamorphic', 'Baked Contact', 'Gneiss'),
           'qdiorite': ('Intrusive', 'Baked Contact', 'Quartz Diorite')}[kind]
    plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95, LAT, LON)
    return {
        'site': site, 'site_alternatives': alt, 'location': LOCATION,
        'result_type': 'i',
        'result_quality': 'b' if site in ('D10', '261', '262', '268') else 'g',
        'method_codes': ':'.join(codes), 'citations': CIT_AGE if kind == 'dyke' else CIT,
        'geologic_classes': geo[0], 'geologic_types': geo[1], 'lithologies': geo[2],
        'lat': LAT, 'lon': LON,
        'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': 0,
        'dir_dec': float(dec), 'dir_inc': float(inc), 'dir_alpha95': a95, 'dir_k': float(k),
        'dir_r': r, 'dir_n_samples': n1, 'dir_n_total_samples': n,
        'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
        'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1),
        'description': desc,
    }


sites = pd.DataFrame([site_row(t) for t in TABLE1])
sites.loc[sites['site_alternatives'] == '', 'site_alternatives'] = None

good = sites[sites['result_quality'] == 'g']
assert len(good) == 9 and good['dir_n_samples'].sum() == 73   # Table 2: 9 sites, 73 samples

# pole: Fisher mean of the site VGPs, reported at the positive-latitude antipode
p = pmag.fisher_mean(ipmag.make_di_block(good['vgp_lon'].tolist(), good['vgp_lat'].tolist()))
if p['inc'] < 0:
    p['dec'], p['inc'] = (p['dec'] + 180) % 360, -p['inc']
# mean direction in the observed polarity of the sites (all W/down)
d = pmag.fisher_mean(ipmag.make_di_block(good['dir_dec'].tolist(), good['dir_inc'].tolist()))

locs = pd.DataFrame([{
    'location': LOCATION, 'location_alternatives': 'central North Greenland dykes:NDL',
    'location_type': 'Region',
    'result_name': 'Victoria Fjord dolerite dykes ca. 1382 Ma pole',
    'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(good['site'].tolist()),
    'method_codes': 'SO-SUN:LP-DIR-AF:LP-DIR-T:DE-FM:DE-VGP:ST-C', 'citations': CIT_AGE,
    'geologic_classes': 'Intrusive', 'lithologies': 'Diabase',
    'lat_s': LAT, 'lat_n': LAT, 'lon_w': LON, 'lon_e': LON,
    'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
    'dir_tilt_correction': 0,
    'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
    'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1), 'dir_n_sites': int(d['n']),
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
    'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'pole_reversed_perc': 0, 'contact_test': 'C+',
    'continent_ocean': 'Greenland', 'country': 'Greenland',
    'description': (
        "Paleomagnetic pole for dolerite dykes cutting Archaean gneiss on the major "
        "(eastern) nunatak at the head of Victoria Fjord, central North Greenland "
        "(Abrahamsen & Van der Voo, 1987). The paper reports 10 dolerite sites from "
        "eight dykes; which sites share a dyke is not stated. The pole is the Fisher "
        "mean of the in-situ virtual geomagnetic poles of the nine accepted dyke "
        "sites, reported at the positive-latitude antipode in present-day Greenland "
        "coordinates (not restored to North America). Dyke D10, oriented by magnetic "
        "compass, gave an anomalous NE direction and is excluded, as are the three "
        "baked host-rock sites. The characteristic remanence, carried by Ti-poor "
        "titanomagnetite, was isolated by alternating-field and, at most sites, "
        "thermal demagnetization as stable end points on orthogonal plots. The "
        "baked-contact test, described by the authors as detailed and positive, is on "
        "the 23 m wide dyke D1: baked gneiss and quartz-diorite within 16 m of the "
        "dyke carry the dyke direction, a core at 32 m is scattered, and gneiss 60 m "
        "from the dyke shows no resetting. The dykes are undated. Their age is "
        "assigned by correlation with the Midsommerso Dolerites of eastern North "
        "Greenland, to which they are petrographically similar and antiparallel in "
        "direction, and one sill of which gave a U-Pb baddeleyite age (Upton et al., "
        "2005). The authors note that the opposite polarities may indicate a slight "
        "difference in intrusive age, and they inferred c. 1.25 Ga from Rb-Sr "
        "isochrons on correlative intrusives.")}])


def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


def validate(path):
    import sys
    root = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
    sys.path.insert(0, os.path.join(root, 'scripts'))
    from validate_magic_contribution import validate_upload_file
    return validate_upload_file(path, tables=['locations', 'sites'])


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Abrahamsen1987_Victoria_Fjord_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    f.write('tab delimited\tlocations\n'); locs.to_csv(f, sep='\t', index=False)
    f.write('>>>>>>>>>>\n')
    f.write('tab delimited\tsites\n'); sites.to_csv(f, sep='\t', index=False)
print(f'-I- Victoria: sites {len(sites)} ({len(good)} accepted), pole '
      f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} K {p["k"]:.1f} N {int(p["n"])}; '
      f'dir {d["dec"]:.1f}/{d["inc"]:.1f} a95 {d["alpha95"]:.1f} k {d["k"]:.1f}')
validate(combined)
