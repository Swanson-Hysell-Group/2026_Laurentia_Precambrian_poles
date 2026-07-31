"""Build the Mistastin Pluton ca. 1441 Ma MagIC contribution (1441_Mistastin).

Combined pole (user decision): the primary quartz-monzonite sites of Fahrig &
Jones (1976) plus the two basement sites of Herve et al. (2015) that escaped
overprinting by the ca. 36 Ma Mistastin Lake impact.

Sources (audited against the papers):
  - Fahrig, W.F. & Jones, D.L. (1976), Can. J. Earth Sci. 13, 832-837
    (doi:10.1139/e76-086): 11 sites, 7 primary. AF (blanket, 25 mT) demagnetization,
    no PCA; quartz-monzonite (adamellite) pluton; six of seven accepted sites
    reversed, one normal.
  - Herve, G. et al. (2015), Earth Planet. Sci. Lett. 417, 151-163
    (doi:10.1016/j.epsl.2015.02.011): basement sites 8 & 16 retain the Mesoproterozoic
    direction; sites 13 & 15, cored within a metre of the impact-melt contact, were
    thermally remagnetized to the ~36 Ma impact direction. Thermal + AF
    demagnetization with principal-component analysis; Ti-poor titanomagnetite.
  - Age: four LA-ICP-MS U-Pb zircon Concordia Ages on the Mistastin batholith target
    basement (Marion & Sylvester, 2010, Planet. Space Sci. 58, 552-573;
    doi:10.1016/j.pss.2009.09.018), verified against the paper: mangerite 1451 +/- 12,
    granodiorite gneiss 1440 +/- 13, anorthosite 1438.7 +/- 8.9, granodiorite
    1429 +/- 10 Ma (all 2 sigma, decay-constant errors included). The four are within
    2 sigma of one another; their 2 sigma envelope is 1419-1463 Ma. Pole age 1441 Ma
    (bracket 1419-1463) is the envelope midpoint, consistent with the notebook / Nordic
    export, and supersedes the Rb-Sr ca. 1346 Ma of Fahrig & Jones (1976).

The pole is the Fisher mean of the in-situ site VGPs from the nine accepted sites
(FJ76 36,37,38,41,42,43,44 + Herve 8,16), reproducing the notebook / Nordic
export (-1.4 / 205.1, A95 7.3, N9). The impact-baking relationship of Herve et al.
(2015) is a qualitative (inverse) argument, not a formal positive baked-contact
test, and the combined set fails a bootstrap reversal test, so no structured
field-test CV column is set (R4 / reversal indeterminate). The six impact- and
secondary-overprinted sites are retained with result_quality 'b'.

Audited fixes: canonical table markers; method_codes separators corrected (';' ->
':', 'LP-AF' -> 'LP-DIR-AF') and Herve PCA line-fitting recorded (DE-BFL);
geologic_classes Igneous -> Intrusive, geologic_types Intrusion -> Pluton;
FJ76 longitude normalized to 0-360 degE; age made consistent at 1441 (1419-1463).
"""
import os, io
from datetime import date
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
SRC = os.path.join(HERE, 'Mistastin_FJ76_Herve_sites_source.txt')

FJ76, HERVE, MS10 = '10.1139/e76-086', '10.1016/j.epsl.2015.02.011', '10.1016/j.pss.2009.09.018'
AGE, AGE_LOW, AGE_HIGH = 1441, 1419, 1463

raw = open(SRC, encoding='utf-8').read().splitlines()
sites = pd.read_csv(io.StringIO('\n'.join(raw[1:])), sep='\t')

# ---- canonical fixes ---------------------------------------------------------
sites['lon'] = sites['lon'] % 360                     # FJ76 -63.59 -> 296.41
sites['geologic_classes'] = 'Intrusive'
sites['geologic_types'] = 'Pluton'
sites['lithologies'] = 'Quartz Monzonite'
sites['method_codes'] = sites['citations'].map(
    lambda c: 'LP-DIR-AF:DE-BLANKET:DE-FM' if c == FJ76
    else 'LP-DIR-T:LP-DIR-AF:DE-BFL:DE-FM')
sites['result_type'] = 'i'
sites['age'] = AGE
sites['age_low'] = AGE_LOW
sites['age_high'] = AGE_HIGH
sites = sites.drop(columns=[c for c in ['age_sigma', 'age_high.1'] if c in sites.columns])

good = sites[sites['result_quality'] == 'g']


def flip_north(m):
    """Report a Fisher mean in the easterly-down (positive-inclination) convention used by the Nordic export."""
    if m['inc'] < 0:
        m['dec'] = (m['dec'] + 180) % 360
        m['inc'] = -m['inc']
    return m


# pole: Fisher mean of site VGPs, reported at the southern-hemisphere pole (Nordic convention: -1.4/205.1)
blk = ipmag.do_flip(di_block=pmag.flip(ipmag.make_di_block(good['vgp_lon'].tolist(),
                                                           good['vgp_lat'].tolist()), combine=True))
p = pmag.fisher_mean(blk)
# mean direction in the easterly-down convention (matches Nordic DEC 90.4 / INC 3.8)
dm = flip_north(pmag.fisher_mean(pmag.flip(ipmag.make_di_block(good['dir_dec'].tolist(),
                                                              good['dir_inc'].tolist()), combine=True)))

# Fahrig & Jones (1976) sites alone, for the combined-vs-FJ76-only note in the description
fj_good = good[good['citations'] == FJ76]
pfj = pmag.fisher_mean(ipmag.do_flip(di_block=pmag.flip(
    ipmag.make_di_block(fj_good['vgp_lon'].tolist(), fj_good['vgp_lat'].tolist()), combine=True)))
sep = float(pmag.angle([p['dec'], p['inc']], [pfj['dec'], pfj['inc']])[0])

locs = pd.DataFrame([{
    'location': 'Mistastin Pluton', 'location_type': 'Region',
    'result_name': 'Mistastin Pluton ca. 1441 Ma pole', 'result_type': 'a',
    'result_quality': 'g', 'sites': ':'.join(map(str, good['site'].tolist())),
    'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BLANKET:DE-BFL:DE-FM:DE-VGP',
    'citations': f'{FJ76}:{HERVE}:{MS10}',
    'geologic_classes': 'Intrusive', 'lithologies': 'Quartz Monzonite',
    'lat_s': round(good['lat'].min(), 3), 'lat_n': round(good['lat'].max(), 3),
    'lon_w': round(good['lon'].min(), 3), 'lon_e': round(good['lon'].max(), 3),
    'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
    'dir_tilt_correction': 0,
    'dir_dec': round(dm['dec'], 1), 'dir_inc': round(dm['inc'], 1),
    'dir_alpha95': round(dm['alpha95'], 1), 'dir_k': round(dm['k'], 1), 'dir_n_sites': int(dm['n']),
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
    'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'continent_ocean': 'North America', 'country': 'Canada',
    'state_province': 'Newfoundland and Labrador', 'terranes': 'Laurentia',
    'description': (
        "Paleomagnetic pole for the ca. 1441 Ma Mistastin Pluton, a quartz-monzonite "
        "(adamellite) intrusion of the anorthosite-quartz-monzonite suite in "
        "Labrador, north of the Grenville Front. The pole is the Fisher mean of the "
        "in-situ site virtual geomagnetic poles from nine sites, combining the "
        "primary, shallow characteristic remanence isolated from seven "
        "quartz-monzonite sites by Fahrig & Jones (1976) using alternating-field "
        "demagnetization with two basement sites of Herve et al. (2015) that escaped "
        "thermal overprinting by the ca. 36 Ma Mistastin Lake impact and retain the "
        "Mesoproterozoic direction (isolated by thermal and alternating-field "
        "demagnetization with principal-component analysis). Computed from the seven "
        f"Fahrig & Jones (1976) sites alone, the pole lies at {abs(pfj['inc']):.1f} S, "
        f"{pfj['dec']:.1f} E (A95 {pfj['alpha95']:.1f} deg), reproducing the previous "
        "Nordic compilation pole (GPMDB 2271); adding the two Herve et al. (2015) "
        f"basement sites shifts the combined pole reported here about {sep:.1f} deg "
        "with little change in A95, and the two poles are statistically "
        "indistinguishable, although the combined nine-site set does not pass a "
        "bootstrap common-mean test between the two studies or a bootstrap reversal "
        "test. Basement sites cored "
        "within a metre of the impact-melt contact (retained here with quality 'b') "
        "were remagnetized to the impact direction, whereas sites farther from the "
        "contact preserve the pre-impact magnetization -- a qualitative impact-baking "
        "relationship rather than a formal baked-contact test. Age from four LA-ICP-MS "
        "U-Pb zircon Concordia Ages spanning 1429-1451 Ma on the target-basement phases "
        "of the Mistastin batholith (mangerite, two granodiorites and anorthosite; "
        "Marion & Sylvester, 2010), whose 2 sigma envelope (1419-1463 Ma) brackets the "
        "adopted pole age of 1441 Ma and supersedes the Rb-Sr ca. 1346 Ma of Fahrig & "
        "Jones (1976).")}])


# ages table: the four LA-ICP-MS U-Pb zircon Concordia Ages of the Mistastin batholith
# target basement (Marion & Sylvester, 2010), replacing the superseded Rb-Sr of legacy 18165.
# Published uncertainties are 2 sigma (decay-constant errors included); age_sigma below is 1 sigma
# (half the published value), per the contribution convention. The four ages are within 2 sigma of
# one another; their 2 sigma envelope (1419-1463 Ma) brackets the adopted pole age of 1441 Ma.
def age_row(nom, sig, sample, rock, extra=''):
    return {'location': 'Mistastin Pluton', 'age': nom, 'age_sigma': sig,
            'age_unit': 'Ma', 'timescale_eon': 'Proterozoic', 'timescale_era': 'Mesoproterozoic',
            'timescale_period': 'Calymmian', 'method_codes': 'GM-UPB',
            'timescale_citations': 'Gradstein et al. 2004', 'citations': MS10,
            'description': ('LA-ICP-MS U-Pb zircon Concordia Age of the ' + rock + ' (sample '
                            + sample + ') of the Mistastin batholith target basement; '
                            + nom.__format__('.1f').rstrip('0').rstrip('.') + ' +/- '
                            + f'{2 * sig:g}' + ' Ma (2 sigma), Marion & Sylvester (2010).' + extra)}


ages = pd.DataFrame([
    age_row(1451, 6.0, 'W05-45', 'mangerite'),
    age_row(1440, 6.5, 'CM001', 'granodiorite gneiss'),
    age_row(1438.7, 4.45, 'CM032', 'anorthosite'),
    age_row(1429, 5.0, 'CM003', 'granodiorite',
            ' The four target-basement Concordia Ages (1429-1451 Ma) are within 2 sigma of '
            'each other; their 2 sigma envelope (1419-1463 Ma) brackets the adopted pole age '
            'and supersedes the Rb-Sr age of Fahrig & Jones (1976).'),
])


def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
write_magic(ages, 'ages', os.path.join(OUT, 'ages.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Mistastin_FJ76_Herve_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    f.write('tab delimited\tlocations\n'); locs.to_csv(f, sep='\t', index=False)
    f.write('>>>>>>>>>>\n')
    f.write('tab delimited\tsites\n'); sites.to_csv(f, sep='\t', index=False)
    f.write('>>>>>>>>>>\n')
    f.write('tab delimited\tages\n'); ages.to_csv(f, sep='\t', index=False)
print(f'-I- Mistastin: sites {len(sites)} ({len(good)} in pole), pole '
      f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} K {p["k"]:.1f} N {int(p["n"])}; '
      f'dir {dm["dec"]:.1f}/{dm["inc"]:.1f} a95 {dm["alpha95"]:.1f}')

import sys
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from validate_magic_contribution import validate_upload_file
validate_upload_file(combined, tables=['locations', 'sites', 'ages'])
