"""Build the Western Channel Diabase (WCD) MagIC contribution.

Source data
-----------
``Irving1972_WCD_sites_source.txt`` — the site-level directions of the Western
Channel Diabase transcribed (by the authoring group, consumed by the notebook)
from Table 1 of Irving, Donaldson & Park (1972), *Canadian Journal of Earth
Sciences* 9, 960-971 (doi:10.1139/e72-080). This is an old (pre-MagIC) study
with no specimen- or measurement-level data available, so the contribution is
built at the site + location level, mirroring the Cleaver Dykes and Lake Shore
Traps contributions. Each site ``dir_k`` was recomputed from the published
``n`` and resultant ``R`` (Table 1), and the site VGPs are the unified-polarity
(northern-hemisphere) VGPs that reproduce the recalculated pole.

What this script does
---------------------
It reads the authored site table and rebuilds a canonical, MagIC-valid
``sites.txt``, then derives the ``locations.txt`` and ``ages.txt`` and the
combined upload file:

- geologic_classes / geologic_types / lithologies are mapped to controlled-
  vocabulary terms (the draft used the non-vocabulary terms ``Volcanic igneous``
  and ``Sedimentary``): diabase sheets -> ``Intrusive`` / ``Sill`` / ``Diabase``,
  diabase dikes -> ``Intrusive`` / ``Volcanic Dike`` / ``Diabase``, baked country
  rock -> ``Sedimentary`` / ``Burnt Sediment`` / ``Baked Sediment``, and the
  combined site means carry both (colon-delimited). Sites 26A/26B, which the
  draft mis-classed ``Sedimentary`` but which sit under Irving's diabase-sheet
  heading with no contact (``c``) suffix, are corrected to diabase sheets;
- method_codes are made canonical (the draft wrote ``LP-DIR-AF: DE-FM`` with a
  stray space); thermal demagnetization (``LP-DIR-T``) is added on the baked
  contact rocks and the Port Radium sill, the units on which Irving et al. ran
  thermal treatment;
- the stray site label ``14 mean`` is normalized to ``14mean`` to match the
  other combined-site means and the location ``sites`` list;
- the recalculated pole is the Fisher mean of the 35 good-site VGPs
  (result_quality ``g``); the 16 subsumed sub-site rows and the reversed Port
  Radium sill (site 144) are flagged ``result_quality='b'`` so
  ``pole_tools.load_magic_sites`` drops them and the pole uses exactly the 35
  sites of the published Irving et al. (1972) mean.

Age control
-----------
The nominal 1592 Ma pole age is anchored by three U-Pb baddeleyite dates on WCD
sheets, carried in ``ages.txt`` (measured dates only): 1590.2 +/- 3.8 Ma
(FA66-42) and 1591.9 +/- 2.9 Ma (FA66-44) of Hamilton & Buchan (2010,
doi:10.1016/j.precamres.2010.06.009), and 1592.4 +/- 2.5 Ma (FA66-45/6) of
Rogers et al. (2018, doi:10.1016/j.lithos.2018.06.002). All three uncertainties
are stated at the 95% (2 sigma) level, so ``age_sigma`` is the 1 sigma (half of
that). Their inverse-variance weighted mean is 1591.8 +/- 1.7 Ma (2 sigma),
i.e. 1592 to the nearest Ma - the two most precise and concordant determinations
(1591.9 and 1592.4) both cluster at ~1592, so 1592 is adopted as the nominal
pole age (updated from the historical 1590 Ma label). The location and site rows
carry the 1589-1594.8 Ma bracket spanning the individual 207Pb/206Pb fraction
ages (Hamilton & Buchan, 2010) and the Rogers et al. (2018) determination.

Field test
----------
Irving et al. (1972) report a positive baked-contact test: baked country rock
within ~0.5 m of WCD contacts carries the diabase direction, while distal rock
(>~8 m) retains a different, pre-intrusion direction (``ST-C``, contact_test
``C+``). The complementary inverse baked-contact test - that the younger WCD
cuts and does not overprint the ca. 1740 Ma Cleaver Dykes - is recorded on the
Cleaver Dykes contribution, not here.

Outputs
-------
``../sites.txt``      — 51 site rows (35 pole sites + 16 excluded), MagIC v3.0.
``../locations.txt``  — one location row: the WCD ca. 1592 Ma pole.
``../ages.txt``       — the three U-Pb dates.
plus the combined upload file ``Western-Channel-Diabase_<DD.Mon.YYYY>.txt``.

Validation: the recreated Fisher mean of the 35 good-site VGPs reproduces the
notebook ``pt.compute_mean_pole`` result (plon 245.5, plat 9.2, K 32.3,
A95 4.3) and the Irving et al. (1972) published pole (9 deg N, 245 deg E,
A95 6.0) to within mutual uncertainty.
"""

import os
import re
import sys
from datetime import date

import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.dirname(HERE)

IRVING_DOI = '10.1139/e72-080'                     # Irving, Donaldson & Park (1972)
HAMILTON_DOI = '10.1016/j.precamres.2010.06.009'   # Hamilton & Buchan (2010)
ROGERS_DOI = '10.1016/j.lithos.2018.06.002'        # Rogers et al. (2018)

LOCATION = 'Western Channel Diabase, Great Bear Lake-west side, NWT'

SITE_COLS = ['site', 'location', 'result_type', 'result_quality', 'method_codes',
             'citations', 'geologic_classes', 'geologic_types', 'lithologies',
             'lat', 'lon', 'age', 'age_sigma', 'age_low', 'age_high', 'age_unit',
             'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95',
             'dir_n_samples', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm',
             'description']


# Irving et al. (1972) distinguish two baked-contact lithologies: baked
# porphyritic lava (sites 20c and 12c -- 'Porphyritic lavas near diabase were
# sampled at 2 sites, one (site 20c) ... and another (site 12c) from a thin
# flow') vs. baked fine-grained sediment (14c, 31c, and the site-29 sediment).
# Their combined site means 12mean (sheet 12 + lava 12c) and 20mean (dike 20 +
# lava 20c) therefore pair diabase with baked porphyritic lava, not sediment.
LAVA_CONTACTS = {'12c', '20c'}
LAVA_MEANS = {'12mean', '20mean'}


def classify(site, description):
    """Return the geologic role of a site row from its label and description."""
    if site == '144':
        return 'port_radium'
    if site in LAVA_CONTACTS:
        return 'lava_contact'
    if re.match(r'^\d+c\d*[a-z]?$', site):    # 14c, 29c1a/b/c, 29c2a/b/c, 31c
        return 'contact'
    if 'mean' in site:
        base = 'dike' if 'dike' in description.lower() else 'sheet'
        suffix = '_lava_mean' if site in LAVA_MEANS else '_mean'
        return base + suffix
    if 'dike' in description.lower():
        return 'dike'
    return 'sheet'


# controlled-vocabulary geology (class, type, lithology) per role
GEOLOGY = {
    'sheet':            ('Intrusive', 'Sill', 'Diabase'),
    'dike':             ('Intrusive', 'Volcanic Dike', 'Diabase'),
    'port_radium':      ('Intrusive', 'Sill', 'Diabase'),
    'contact':          ('Sedimentary', 'Burnt Sediment', 'Baked Sediment'),
    'lava_contact':     ('Extrusive', 'Lava Flow', 'Porphyritic Basalt'),
    'sheet_mean':       ('Intrusive:Sedimentary', 'Sill:Burnt Sediment',
                         'Diabase:Baked Sediment'),
    'dike_mean':        ('Intrusive:Sedimentary', 'Volcanic Dike:Burnt Sediment',
                         'Diabase:Baked Sediment'),
    'sheet_lava_mean':  ('Intrusive:Extrusive', 'Sill:Lava Flow',
                         'Diabase:Porphyritic Basalt'),
    'dike_lava_mean':   ('Intrusive:Extrusive', 'Volcanic Dike:Lava Flow',
                         'Diabase:Porphyritic Basalt'),
}

# thermal demagnetization was applied to the baked contacts and the Port Radium
# sill; the combined means include a thermally treated contact sample
THERMAL_ROLES = {'contact', 'lava_contact', 'port_radium', 'sheet_mean',
                 'dike_mean', 'sheet_lava_mean', 'dike_lava_mean'}


def method_codes_for(role):
    codes = ['LP-DIR-AF']
    if role in THERMAL_ROLES:
        codes.append('LP-DIR-T')
    codes.append('DE-FM')
    return ':'.join(codes)


def build_sites(src):
    rows = []
    for _, r in src.iterrows():
        site = '14mean' if r['site'] == '14 mean' else r['site']
        role = classify(site, r['description'])
        gc, gt, lith = GEOLOGY[role]
        row = {c: r.get(c, '') for c in SITE_COLS}
        row.update({
            'site': site,
            'location': LOCATION,
            'result_type': 'i',
            'method_codes': method_codes_for(role),
            'citations': IRVING_DOI,
            'geologic_classes': gc,
            'geologic_types': gt,
            'lithologies': lith,
            'age': '', 'age_sigma': '', 'age_low': 1589, 'age_high': 1594.8,
            'age_unit': 'Ma', 'dir_tilt_correction': 0,
        })
        rows.append(row)
    return pd.DataFrame(rows, columns=SITE_COLS).fillna('')


def compute_pole(good):
    vgps = good[['vgp_lon', 'vgp_lat']].astype(float).values.tolist()
    return pmag.fisher_mean(vgps)


def build_locations(sites_df):
    good = sites_df[sites_df['result_quality'] == 'g'].copy()
    pole = compute_pole(good)
    meandir = ipmag.fisher_mean(dec=good['dir_dec'].astype(float).tolist(),
                                inc=good['dir_inc'].astype(float).tolist())
    lat = good['lat'].astype(float)
    lon = good['lon'].astype(float)

    print('Recalculated Western Channel Diabase pole '
          '(Fisher mean of %d site VGPs):' % pole['n'])
    print('  plon %.1f  plat %.1f  K %.1f  A95 %.1f'
          % (pole['dec'], pole['inc'], pole['k'], pole['alpha95']))
    print('  mean dir: dec %.1f inc %.1f a95 %.1f k %.1f'
          % (meandir['dec'], meandir['inc'], meandir['alpha95'], meandir['k']))
    print('  notebook pt.compute_mean_pole: plon 245.5 plat 9.2 K 32.3 A95 4.3')
    print('  published Irving et al. (1972): plat 9  plon 245 (115 W)  A95 6.0')

    loc = {
        'location': LOCATION, 'location_type': 'Region',
        'result_name': 'Western Channel Diabase ca. 1592 Ma pole',
        'result_type': 'a', 'result_quality': 'g',
        'method_codes': ('LP-DIR-AF:LP-DIR-T:DE-FM:DE-VGP:DA-DIR-GEO:'
                         'GM-UPB-CC:ST-C'),
        'citations': ':'.join([IRVING_DOI, HAMILTON_DOI, ROGERS_DOI]),
        'geologic_classes': 'Intrusive', 'lithologies': 'Diabase',
        'lat_s': round(lat.min(), 2), 'lat_n': round(lat.max(), 2),
        'lon_w': round(lon.min(), 2), 'lon_e': round(lon.max(), 2),
        'continent_ocean': 'North America', 'country': 'Canada',
        'age': 1592, 'age_low': 1589, 'age_high': 1594.8, 'age_unit': 'Ma',
        'dir_tilt_correction': 0,
        'dir_dec': round(meandir['dec'], 1), 'dir_inc': round(meandir['inc'], 1),
        'dir_alpha95': round(meandir['alpha95'], 1),
        'dir_k': round(meandir['k'], 1), 'dir_n_sites': int(meandir['n']),
        'pole_lat': round(pole['inc'], 1), 'pole_lon': round(pole['dec'], 1),
        'pole_alpha95': round(pole['alpha95'], 1), 'pole_k': round(pole['k'], 1),
        'pole_n_sites': int(pole['n']),
        'contact_test': 'C+',
        'sites': ':'.join(good['site'].tolist()),
        'description': (
            'Paleomagnetic pole for the ca. 1592 Ma Western Channel Diabase, a '
            'suite of diabase sheets (sills) and NNE-trending dikes of the Great '
            'Bear magmatic zone (Wopmay Orogen) that intrude the Cameron Bay and '
            'Echo Bay groups and the overlying Hornby Bay Group on the west side '
            'of Great Bear Lake, NWT, Canada (Irving et al., 1972). Pole is the '
            'Fisher mean of the VGPs of the 35 site means retained by Irving et '
            'al. (1972; sheets plus dikes); the reversed-polarity Port Radium '
            'sill (site 144) and the 16 subsumed sub-site rows are excluded from '
            'the mean. Directions are in geographic coordinates; the sheets and '
            'dikes are interpreted to not require tilt correction. A positive '
            'baked-contact test '
            '(baked country rock within ~0.5 m of contacts carries the diabase '
            'direction while distal rock does not) supports a primary, '
            'emplacement-age magnetization (ST-C). Age from U-Pb baddeleyite '
            'dating of WCD sheets (Hamilton and Buchan, 2010; Rogers et al., '
            '2018).'),
    }
    return pd.DataFrame([loc]).fillna('')


AGE_COLS = ['location', 'age', 'age_sigma', 'age_low', 'age_high', 'age_unit',
            'method_codes', 'citations', 'timescale_eon', 'timescale_era',
            'timescale_period', 'description']


def build_ages():
    """The three measured U-Pb baddeleyite dates on Western Channel Diabase
    sheets. Uncertainties in the sources are 95% (2 sigma); age_sigma is the
    1 sigma (half). Only specific determinations go here (not inferred ranges).
    """
    rows = [
        {'location': LOCATION, 'age': '1590.2', 'age_sigma': '1.9',
         'age_low': '', 'age_high': '', 'age_unit': 'Ma',
         'method_codes': 'GM-UPB-CC', 'citations': HAMILTON_DOI,
         'timescale_eon': 'Proterozoic', 'timescale_era': 'Mesoproterozoic',
         'timescale_period': 'Calymmian',
         'description': ('WCD sheet FA66-42 (MacAlpine Channel), weighted-mean '
                         '207Pb/206Pb U-Pb single-grain baddeleyite 1590.2 +/- '
                         '3.8 Ma (95% conf.; Hamilton and Buchan, 2010)')},
        {'location': LOCATION, 'age': '1591.9', 'age_sigma': '1.45',
         'age_low': '', 'age_high': '', 'age_unit': 'Ma',
         'method_codes': 'GM-UPB-CC', 'citations': HAMILTON_DOI,
         'timescale_eon': 'Proterozoic', 'timescale_era': 'Mesoproterozoic',
         'timescale_period': 'Calymmian',
         'description': ('WCD sheet FA66-44 (Hornby Bay), weighted-mean '
                         '207Pb/206Pb U-Pb baddeleyite 1591.9 +/- 2.9 Ma '
                         '(95% conf.; Hamilton and Buchan, 2010)')},
        {'location': LOCATION, 'age': '1592.4', 'age_sigma': '1.25',
         'age_low': '', 'age_high': '', 'age_unit': 'Ma',
         'method_codes': 'GM-UPB', 'citations': ROGERS_DOI,
         'timescale_eon': 'Proterozoic', 'timescale_era': 'Mesoproterozoic',
         'timescale_period': 'Calymmian',
         'description': ('WCD sample FA66-045, U-Pb baddeleyite upper-intercept '
                         '1592.4 +/- 2.5 Ma (95% conf.; Rogers et al., 2018)')},
    ]
    return pd.DataFrame(rows, columns=AGE_COLS)


def write_table(df, path, table_name):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{table_name}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')
    print(f'-I- wrote {path} ({len(df)} rows)')


def write_combined_upload(loc_df, sites_df, ages_df):
    stamp = date.today().strftime('%d.%b.%Y')
    out = os.path.join(OUT_DIR, f'Western-Channel-Diabase_{stamp}.txt')
    with open(out, 'w', encoding='utf-8') as f:
        for table, df in [('locations', loc_df), ('sites', sites_df),
                          ('ages', ages_df)]:
            if table != 'locations':
                f.write('>>>>>>>>>>\n')
            f.write(f'tab delimited\t{table}\n')
            df.to_csv(f, sep='\t', index=False)
    print(f'-I- wrote {out} (locations + sites + ages, canonical markers)')
    return out


def _validate(upload_path, tables):
    repo_root = os.path.dirname(os.path.dirname(OUT_DIR))
    sys.path.insert(0, os.path.join(repo_root, 'scripts'))
    try:
        from validate_magic_contribution import validate_upload_file
    except Exception as exc:  # pragma: no cover
        print(f'-W- could not import validator ({exc}); skipping validation')
        return
    if not validate_upload_file(upload_path, tables=tables):
        print('-W- validation reported problems above; review before uploading')


def main():
    src = pd.read_csv(os.path.join(HERE, 'Irving1972_WCD_sites_source.txt'),
                      sep='\t', skiprows=1, dtype=str,
                      encoding='utf-8').fillna('')
    sites_df = build_sites(src)
    write_table(sites_df, os.path.join(OUT_DIR, 'sites.txt'), 'sites')

    loc_df = build_locations(sites_df)
    write_table(loc_df, os.path.join(OUT_DIR, 'locations.txt'), 'locations')

    ages_df = build_ages()
    write_table(ages_df, os.path.join(OUT_DIR, 'ages.txt'), 'ages')

    out_upload = write_combined_upload(loc_df, sites_df, ages_df)
    _validate(out_upload, tables=['locations', 'sites', 'ages'])


if __name__ == '__main__':
    main()
