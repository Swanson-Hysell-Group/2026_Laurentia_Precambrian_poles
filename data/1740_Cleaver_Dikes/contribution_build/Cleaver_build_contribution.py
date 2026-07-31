"""Build the Cleaver Dykes MagIC contribution (sites.txt + locations.txt).

Source data
-----------
``Irving2004_Cleaver_Table1.csv`` — the site-level directions of the Cleaver
Dykes and associated rocks transcribed from Table 1 of Irving, Baker, Hamilton
& Wynne (2004), *Precambrian Research* 129, 251-270
(doi:10.1016/j.precamres.2003.10.002). This is an old (pre-MagIC) study with no
specimen- or measurement-level data available, so the contribution is built at
the site + location level, mirroring the legacy Lake Shore Traps contributions.

Provenance
----------
``student_contribution_11.May.2026.txt`` is the draft MagIC
contribution assembled by an undergraduate student from the same Table 1.
This script rebuilds the contribution from the published table and improves it
for submission to the MagIC database:

- geologic_classes set to the controlled-vocabulary terms ``Intrusive`` /
  ``Extrusive`` (the draft used the non-vocabulary term ``Igneous``);
- the excluded Echo Bay re-sampling (site 172) is flagged ``result_quality='b'``
  so ``pole_tools.load_magic_sites`` drops it and the pole uses exactly the 17
  dykes of the published mean (B = 17);
- the location pole result carries the positive baked-contact (``ST-C``) and
  inverse baked-contact (``ST-C-I``) test codes that match the GPMDB 9139 /
  Nordic compilation field-test record (``C+, C*+``); the conglomerate test
  (``ST-G``) is recorded on the Labine Group location it tests, not on the
  Cleaver pole;
- per-site VGPs (``vgp_lat``/``vgp_lon``/``vgp_dp``/``vgp_dm``) are recomputed
  from the published directions with ``pmag.dia_vgp``;
- the Labine Group baked-contact/conglomerate controls and the ca. 780 Ma Hottah
  (Gunbarrel) Gabbro are kept in their own location rows to document the field
  tests, with the high-precision Labine Group ages of Ootes et al. (2015,
  doi:10.1139/cjes-2015-0026).

Outputs
-------
``../sites.txt``     — all 24 site rows (17 pole dykes + excluded site 172 +
                       5 Labine controls + Hottah Gabbro), MagIC v3.0.
``../locations.txt`` — three location rows: the Cleaver Dykes pole result, the
                       Labine Group test container, and the Hottah Gabbro.

Validation: the recreated Fisher mean of the 17 site VGPs reproduces the
published Irving et al. (2004) Cleaver Dykes paleopole (19.4 deg N, 276.7 deg E,
A95 6.1) to within ~1.5 deg, i.e. within mutual A95. The recreated mean direction
is slightly more scattered than the published grand mean (k 50 vs 64): Irving et
al. computed their mean from the underlying specimen-level line fits (213
specimens), not from the 17 rounded site means tabulated in Table 1, and that
difference -- not rounding of the table, which a +/-0.05 deg Monte Carlo
perturbation shows shifts the pole by only ~0.01 deg and leaves k unchanged --
accounts for the small offset and the slightly larger A95. The recreated value is
the reproducible result of computing the pole strictly from the published table.
"""

import os
import sys
from datetime import date

import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.dirname(HERE)

IRVING_DOI = '10.1016/j.precamres.2003.10.002'
OOTES_DOI = '10.1139/cjes-2015-0026'  # Ootes et al. (2015) CJES 52, 1062-1092
HARLAN_DOI = '10.1130/g19944.1'  # Harlan et al. (2003) Geology 31, 1053-1056

# demagnetization / analysis method codes common to every site (AF + thermal
# stepwise demagnetization, PCA best-fit lines, Fisher site means, geographic
# coordinates)
SITE_METHODS = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DA-DIR-GEO'

# per-role geologic classification (controlled-vocabulary terms)
GEOLOGY = {
    'pole':              ('Intrusive', 'Volcanic Dike', 'Diabase'),
    'excluded':          ('Intrusive', 'Volcanic Dike', 'Diabase'),
    'conglomerate':      ('Extrusive', 'Volcanic Other', 'Volcaniclastic Breccia'),
    'baked_contact':     None,   # set per lithology below
    'baked_distal':      None,
    'younger_intrusion': ('Intrusive', 'Sill', 'Gabbro'),
}


def geology_for(row):
    """Return (geologic_classes, geologic_types, lithologies) for a CSV row."""
    role = row['role']
    if role in ('baked_contact', 'baked_distal'):
        # Labine Group controls: tuff sites 32/33, andesite (lava flow) 35/36
        if 'tuff' in row['note'].lower():
            return ('Extrusive', 'Tuff', 'Tuff')
        return ('Extrusive', 'Lava Flow', 'Andesite')
    return GEOLOGY[role]


def age_fields_for(row):
    """Return dict of age columns for a CSV row (age | age_low/high | sigma)."""
    role, note = row['role'], row['note'].lower()
    if row['location'] == 'Cleaver Dykes':
        return {'age': '', 'age_sigma': '', 'age_low': 1736, 'age_high': 1745,
                'age_unit': 'Ma'}
    if row['location'] == 'Hottah Gabbro':
        # Gunbarrel gabbro (named for Gunbarrel Inlet, SE Great Bear Lake -- the
        # 'Gunbarrel Gabbro' of Park et al. 1995 = this Hottah sheet): U-Pb
        # baddeleyite 779.6 +/- 1.4 Ma (95% conf.; Harlan et al., 2003), so
        # age_sigma = 0.7 Ma (1sigma). NOT the Calder gabbro (779.5 +/- 1.8), a
        # different Gunbarrel-event sheet.
        return {'age': 779.6, 'age_sigma': 0.7, 'age_low': '', 'age_high': '',
                'age_unit': 'Ma'}
    # Labine Group controls (Ootes et al., 2015 U-Pb); the cited weighted-mean
    # ages are 95% confidence, so age_sigma is the 1sigma (half of that)
    if 'tuff' in note:
        return {'age': 1876.4, 'age_sigma': 1.2, 'age_low': '', 'age_high': '',
                'age_unit': 'Ma'}
    if 'andesite' in note and row['role'] != 'conglomerate':
        return {'age': 1872.6, 'age_sigma': 0.45, 'age_low': '', 'age_high': '',
                'age_unit': 'Ma'}
    # volcanic breccia (conglomerate test) — broad Labine Group age
    return {'age': '', 'age_sigma': '', 'age_low': 1860, 'age_high': 1880,
            'age_unit': 'Ma'}


def citations_for(row):
    if row['location'] == 'Labine Group' and row['role'] != 'conglomerate':
        return f'{IRVING_DOI}:{OOTES_DOI}'
    if row['role'] == 'conglomerate':
        return f'{IRVING_DOI}:{OOTES_DOI}'
    if row['location'] == 'Hottah Gabbro':
        # Irving et al. (2004) direction + Harlan et al. (2003) U-Pb age
        return f'{IRVING_DOI}:{HARLAN_DOI}'
    return IRVING_DOI


def build_sites(df):
    rows = []
    for _, r in df.iterrows():
        gc, gt, lith = geology_for(r)
        elon = round(360.0 - float(r['lon_W']), 2)
        dec, inc, a95 = float(r['dec']), float(r['inc']), float(r['a95'])
        plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95, float(r['lat']), elon)
        ages = age_fields_for(r)
        rows.append({
            'age': ages['age'], 'age_high': ages['age_high'],
            'age_low': ages['age_low'], 'age_sigma': ages['age_sigma'],
            'age_unit': ages['age_unit'],
            'citations': citations_for(r),
            'description': r['note'],
            'dir_alpha95': a95, 'dir_dec': dec, 'dir_inc': inc,
            'dir_k': float(r['k']),
            'dir_n_samples': int(r['n_cores']),
            'dir_n_specimens': int(r['n_spec']),
            'dir_polarity': 'n' if r['location'] == 'Cleaver Dykes' else '',
            'dir_tilt_correction': 0,
            'geologic_classes': gc, 'geologic_types': gt,
            'lat': float(r['lat']), 'lithologies': lith,
            'location': r['location'], 'lon': elon,
            'method_codes': SITE_METHODS,
            'result_quality': r['quality'], 'result_type': 'i',
            'site': str(r['site']),
            'vgp_dm': round(dm, 1), 'vgp_dp': round(dp, 1),
            'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
        })
    return pd.DataFrame(rows)


def compute_pole(df_pole):
    vgps = df_pole[['vgp_lon', 'vgp_lat']].values.tolist()
    return pmag.fisher_mean(vgps)


def build_locations(df, sites_df):
    pole_sites = df[(df['location'] == 'Cleaver Dykes') &
                    (df['quality'] == 'g')]
    pole_site_rows = sites_df[sites_df['site'].isin(pole_sites['site'].astype(str))]
    pole = compute_pole(pole_site_rows)
    meandir = ipmag.fisher_mean(dec=pole_site_rows['dir_dec'].astype(float).tolist(),
                                inc=pole_site_rows['dir_inc'].astype(float).tolist())
    n_samples = int(pole_site_rows['dir_n_samples'].sum())
    n_specimens = int(pole_site_rows['dir_n_specimens'].sum())
    print('Recreated Cleaver Dykes pole (Fisher mean of %d site VGPs):' % pole['n'])
    print('  plon %.1f  plat %.1f  K %.1f  A95 %.1f  N_samples %d'
          % (pole['dec'], pole['inc'], pole['k'], pole['alpha95'], n_samples))
    print('  mean dir: dec %.1f inc %.1f a95 %.1f k %.1f'
          % (meandir['dec'], meandir['inc'], meandir['alpha95'], meandir['k']))
    print('  published Irving et al. (2004): plon 276.7  plat 19.4  A95 6.1  B 17')

    def bbox(sites):
        rows = sites_df[sites_df['site'].isin([str(s) for s in sites])]
        return (rows['lat'].max(), rows['lat'].min(),
                rows['lon'].max(), rows['lon'].min())

    cleaver_sites = list(pole_sites['site'].astype(str))
    lat_n, lat_s, lon_e, lon_w = bbox(cleaver_sites)

    cleaver = {
        'age': '', 'age_high': 1745, 'age_low': 1736, 'age_sigma': '',
        'age_unit': 'Ma', 'citations': IRVING_DOI,
        'geologic_classes': 'Intrusive', 'lat_n': lat_n, 'lat_s': lat_s,
        'lithologies': 'Diabase', 'location': 'Cleaver Dykes',
        'location_type': 'Region', 'lon_e': lon_e, 'lon_w': lon_w,
        'continent_ocean': 'North America', 'country': 'Canada',
        'description': (
            'Paleomagnetic pole for the ca. 1740 Ma Cleaver Dykes, a swarm of '
            'subvertical, post-orogenic diabase dykes of the Great Bear magmatic '
            'zone (Hottah terrane), Wopmay Orogen, northwest Canadian Shield '
            '(Irving et al., 2004). U-Pb baddeleyite age 1740 +5/-4 Ma. Pole is '
            'the Fisher mean of the VGPs of the 17 Cleaver dyke site means of the '
            "published 'Cleaver magnetization' (single, SE-and-down 'normal' "
            'polarity carried by fine-grained magnetite, thermal unblocking '
            '500-600 C). The Echo Bay re-sampling (site 172, the resampling of '
            'site 144 of Irving et al. 1972a) is excluded from the mean. Directions '
            'are in geographic coordinates; the subvertical dykes require no tilt '
            'correction. Positive baked-contact test (dyke 34/22 bake adjacent '
            'Labine Group andesite/tuff) and positive inverse baked-contact test '
            '(the younger ca. 1590 Ma Western Channel Diabase and ca. 780 Ma Hottah '
            'Gabbro that cut the dykes carry different directions) support that '
            'the magnetization is primary.'),
        'dir_tilt_correction': 0,
        'method_codes': ('LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP:DA-DIR-GEO:'
                         'GM-UPB-CC:ST-C:ST-C-I'),
        # location-level mean direction (geographic; underlies the pole)
        'dir_dec': round(meandir['dec'], 1), 'dir_inc': round(meandir['inc'], 1),
        'dir_alpha95': round(meandir['alpha95'], 1), 'dir_k': round(meandir['k'], 1),
        # stored as str so the count columns (absent on the container rows) are
        # not upcast to float (99.0) when the location rows are combined
        'dir_n_sites': str(int(meandir['n'])), 'dir_n_samples': str(n_samples),
        'dir_n_specimens': str(n_specimens),
        'pole_alpha95': round(pole['alpha95'], 1), 'pole_k': round(pole['k'], 1),
        'pole_lat': round(pole['inc'], 1), 'pole_lon': round(pole['dec'], 1),
        'pole_n_sites': int(pole['n']),
        # structured field-test result (controlled vocabulary): positive
        # baked-contact test (clean at dyke 34). method_codes also carries the
        # positive inverse contact test (ST-C-I).
        'contact_test': 'C+',
        'result_name': 'Cleaver Dykes ca. 1740 Ma pole',
        'result_quality': 'g', 'result_type': 'a',
        'sites': ':'.join(cleaver_sites),
    }

    labine_sites = ['24', '32', '33', '35', '36']
    lat_n2, lat_s2, lon_e2, lon_w2 = bbox(labine_sites)
    labine = {
        'age': '', 'age_high': 1880, 'age_low': 1860, 'age_sigma': '',
        'age_unit': 'Ma', 'citations': f'{IRVING_DOI}:{OOTES_DOI}',
        'geologic_classes': 'Extrusive', 'lat_n': lat_n2, 'lat_s': lat_s2,
        'lithologies': 'Tuff:Andesite:Volcaniclastic Breccia',
        'location': 'Labine Group', 'location_type': 'Region',
        'lon_e': lon_e2, 'lon_w': lon_w2,
        'continent_ocean': 'North America', 'country': 'Canada',
        'description': (
            'Labine Group (Great Bear magmatic zone, ca. 1875 Ma; Ootes et al., '
            '2015) host rocks used for the Cleaver Dykes field tests (Irving et '
            'al., 2004). The clearest positive baked-contact test is at dyke 34: '
            'baked andesite at the contact (site 35) carries the Cleaver '
            'magnetization while andesite 50 m away (site 36) retains a distinct '
            'pre-Cleaver direction (ST-C). A companion test is recorded at dyke 22 '
            'against Labine tuff (contact site 32, distal site 33); note that the '
            'Table 1 directions of Irving et al. (2004) for the dyke-22 tuff pair '
            'are inconsistent with the accompanying text (which reports the '
            'contact tuff in agreement with the dyke), so only the dyke-34 test is '
            'treated as unambiguous here. Site directions are archived exactly as '
            'published in Table 1. Andesite boulders in a Labine Group volcanic '
            'breccia (site 24) show dispersed directions (conglomerate test, '
            'ST-G), indicating the host has not been regionally remagnetized.'),
        'dir_tilt_correction': 0,
        'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:ST-C:ST-G',
        'pole_alpha95': '', 'pole_k': '', 'pole_lat': '', 'pole_lon': '',
        'pole_n_sites': '',
        'result_name': 'Labine Group baked-contact and conglomerate tests',
        'result_quality': 'g', 'result_type': 'a',
        'sites': ':'.join(labine_sites),
    }

    lat_n3, lat_s3, lon_e3, lon_w3 = bbox(['45'])
    hottah = {
        'age': 779.6, 'age_high': '', 'age_low': '', 'age_sigma': 0.7,
        'age_unit': 'Ma', 'citations': f'{IRVING_DOI}:{HARLAN_DOI}',
        'geologic_classes': 'Intrusive', 'lat_n': lat_n3, 'lat_s': lat_s3,
        'lithologies': 'Gabbro', 'location': 'Hottah Gabbro',
        'location_type': 'Region', 'lon_e': lon_e3, 'lon_w': lon_w3,
        'continent_ocean': 'North America', 'country': 'Canada',
        'description': (
            "Hottah Gabbro sheet (the 'Gunbarrel Gabbro' of Park et al., 1995), a "
            'ca. 780 Ma intrusion dated at Gunbarrel Inlet (SE Great Bear Lake) '
            'by U-Pb baddeleyite at 779.6 +/- 1.4 Ma (95% conf.; the Gunbarrel '
            'gabbro of Harlan et al., 2003) that post-dates and '
            'cuts the Cleaver Dykes. Its magnetization direction differs from the '
            'Cleaver magnetization, contributing to the positive inverse '
            'baked-contact test (the dykes were not overprinted by this younger '
            'event).'),
        'dir_tilt_correction': 0,
        'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM',
        'pole_alpha95': '', 'pole_k': '', 'pole_lat': '', 'pole_lon': '',
        'pole_n_sites': '',
        'result_name': 'Hottah (Gunbarrel) Gabbro ca. 780 Ma',
        'result_quality': 'g', 'result_type': 'a',
        'sites': '45',
    }
    # the pole row carries dir_*/contact_test columns the container rows lack;
    # union the columns and leave them empty for Labine/Hottah
    return pd.DataFrame([cleaver, labine, hottah]).fillna('')


AGE_COLS = ['location', 'age', 'age_sigma', 'age_low', 'age_high', 'age_unit',
            'method_codes', 'citations', 'timescale_eon', 'timescale_era',
            'timescale_period', 'description']


def build_ages():
    """The measured radiometric ages that anchor the three locations.

    Only specific determinations go here (not inferred brackets): the Cleaver
    dyke U-Pb date, the two Labine Group host-rock dates, and the Gunbarrel
    (Hottah) gabbro date. Uncertainties from the sources are 95%-confidence
    weighted means; age_sigma is the 1sigma (half of that), while the published
    95%-conf value is stated in the description. The Cleaver date is an
    asymmetric upper-intercept (+5/-4), which age_sigma cannot express, so it is
    carried as age_low/age_high -- measurement uncertainty on a specific date,
    not an inferred range.
    """
    rows = [
        {'location': 'Cleaver Dykes', 'age': '1740', 'age_sigma': '',
         'age_low': '1736', 'age_high': '1745', 'age_unit': 'Ma',
         'method_codes': 'GM-UPB-CC-T0', 'citations': IRVING_DOI,
         'timescale_eon': 'Proterozoic', 'timescale_era': 'Paleoproterozoic',
         'timescale_period': 'Statherian',
         'description': ('U-Pb baddeleyite, modified York regression upper '
                         'intercept 1740 +5/-4 Ma (Irving et al., 2004)')},
        {'location': 'Labine Group', 'age': '1876.4', 'age_sigma': '1.2',
         'age_low': '', 'age_high': '', 'age_unit': 'Ma',
         'method_codes': 'GM-UPB', 'citations': OOTES_DOI,
         'timescale_eon': 'Proterozoic', 'timescale_era': 'Paleoproterozoic',
         'timescale_period': 'Orosirian',
         'description': ('Labine Group rhyolite tuff, weighted-mean 207Pb/206Pb '
                         'U-Pb 1876.4 +/- 2.4 Ma (95% conf.; Ootes et al., 2015)')},
        {'location': 'Labine Group', 'age': '1872.6', 'age_sigma': '0.45',
         'age_low': '', 'age_high': '', 'age_unit': 'Ma',
         'method_codes': 'GM-UPB', 'citations': OOTES_DOI,
         'timescale_eon': 'Proterozoic', 'timescale_era': 'Paleoproterozoic',
         'timescale_period': 'Orosirian',
         'description': ('Labine Group andesite, U-Pb 1872.6 +/- 0.9 Ma '
                         '(95% conf.; Ootes et al., 2015)')},
        {'location': 'Hottah Gabbro', 'age': '779.6', 'age_sigma': '0.7',
         'age_low': '', 'age_high': '', 'age_unit': 'Ma',
         'method_codes': 'GM-UPB', 'citations': HARLAN_DOI,
         'timescale_eon': 'Proterozoic', 'timescale_era': 'Neoproterozoic',
         'timescale_period': 'Tonian',
         'description': ('Gunbarrel gabbro (Gunbarrel Inlet, SE Great Bear '
                         'Lake), U-Pb baddeleyite 779.6 +/- 1.4 Ma (95% conf.; '
                         'Harlan et al., 2003)')},
    ]
    return pd.DataFrame(rows, columns=AGE_COLS)


def write_table(df, path, table_name):
    with open(path, 'w') as f:
        f.write(f'tab delimited\t{table_name}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')
    print(f'-I- wrote {path} ({len(df)} rows)')


def main():
    csv_path = os.path.join(HERE, 'Irving2004_Cleaver_Table1.csv')
    df = pd.read_csv(csv_path, comment='#', dtype=str)
    df['note'] = df['note'].fillna('')

    sites_df = build_sites(df)
    write_table(sites_df, os.path.join(OUT_DIR, 'sites.txt'), 'sites')

    loc_df = build_locations(df, sites_df)
    write_table(loc_df, os.path.join(OUT_DIR, 'locations.txt'), 'locations')

    ages_df = build_ages()
    write_table(ages_df, os.path.join(OUT_DIR, 'ages.txt'), 'ages')

    out_upload = write_combined_upload(loc_df, sites_df, ages_df)

    # validate the combined upload file before submission: PmagPy data-model
    # checks + the age-completeness rule MagIC enforces server-side. locations,
    # sites, and ages are all authored here (built from the published tables).
    _validate(out_upload, tables=['locations', 'sites', 'ages'])


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


def write_combined_upload(loc_df, sites_df, ages_df):
    """Write the combined MagIC upload file (locations + sites + ages) for the
    three Cleaver locations, using canonical ``tab delimited`` markers and the
    ``>>>>>>>>>>`` table separator. Named ``<Location-names>_<DD.Mon.YYYY>.txt``.
    """
    names = '_'.join(loc.replace(' ', '-') for loc in sorted(loc_df['location']))
    stamp = date.today().strftime('%d.%b.%Y')
    out = os.path.join(OUT_DIR, f'{names}_{stamp}.txt')
    with open(out, 'w') as f:
        for table, df in [('locations', loc_df), ('sites', sites_df),
                          ('ages', ages_df)]:
            if table != 'locations':
                f.write('>>>>>>>>>>\n')
            f.write(f'tab delimited\t{table}\n')
            df.to_csv(f, sep='\t', index=False)
    print(f'-I- wrote {out} (locations + sites + ages, canonical markers)')
    return out


if __name__ == '__main__':
    main()
