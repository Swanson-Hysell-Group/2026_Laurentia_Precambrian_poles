"""Build the MagIC contribution for the Freda Formation data and upper Freda pole of Fuentes et al. (2025).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026).

Fuentes, A. J., Fairchild, L. M., Hodgin, E. B., Alemu, T., & Swanson-Hysell, N. L. (2025).
Termination of Laurentia's rapid plate motion at the start of the Grenvillian Orogeny.
Journal of Geophysical Research: Solid Earth 130, e2025JB031794. doi:10.1029/2025JB031794

The specimen interpretations and measurements are those archived at earthref.org/MagIC/20259
(read from ``../previous_MagIC/magic_contribution_20259.txt``). The measurements table is
copied through unchanged. The study code is the repository of Fuentes et al. (2025)
(github.com/Swanson-Hysell-Group/Upper_Freda_Pole, doi:10.5281/zenodo.14728751).

Structure and interpretation
----------------------------
- Each sample was collected at its own stratigraphic level and is treated as an independent
  reading of the geomagnetic field (Fuentes et al., 2025, following Sapienza et al., 2023).
  Each sample is therefore a site (site name = sample name) and each sampled stratigraphic
  section is a location. Site rows give the specimen direction of each component (HT,
  high-temperature specularite detrital remanence; MT, mid-temperature pigmentary hematite
  chemical remanence) in geographic and tilt-corrected coordinates, with the virtual
  geomagnetic pole of the tilt-corrected direction and the sample's stratigraphic height in
  the section (m).
- Site coordinates are the section coordinates of Fuentes et al. (2025, Table S1).
- Specimen WRF.76a is linked to sample WRF.76a (it was linked to WRF.75a). Sites whose
  specimen gave no interpretable component, or whose fit is flagged 'b', are in no mean.
- Location rows: for every section, a row listing all of its sites and the tilt-corrected
  HT and MT section means (sites flagged 'b' excluded). The Wolf Point section also carries
  three poles, each the Fisher mean of site VGPs or derived from them:
  (1) HT as measured, the compilation pole (uncorrected pole of Fuentes et al., 2025,
      Table 1 note);
  (2) the SVEI inclination-shallowing-corrected Kent mean pole of Fuentes et al. (2025,
      Table 1; f 0.59-0.87 propagated following Pierce et al., 2022), stored as published;
      pole_alpha95 is the radius of the circle with the ellipse's area, sqrt(zeta95 * eta95);
  (3) MT chemical remanence pole (Table 1 note).
  The VGPs are calculated at the Table S1 coordinates; the poles equal the published values
  (the study code used the rounded locality 47.0 N, 89.0 W).
- Field tests: the bootstrap fold test on the combined Freda Formation data of this study
  and Henry et al. (1977) is positive for both HT (95-102% unfolding) and MT (97-104%)
  components, recorded as fold_test F+ on the Wolf Point poles.
- Ages are depositional ages from the Bayesian thermal-subsidence (SubsidenceChron.jl) age
  model of the Oronto Group of Fuentes et al. (2025): model mean and 95% credible interval
  at each section's stratigraphic height (Bad River 2200 m; White River Flowage 3720 m,
  spanning its 2620-4820 m height range; Lake Road 5800 m; Wolf Point 5000-6000 m, the
  pole age of Table 1). The ages table holds the CA-ID-TIMS detrital zircon maximum
  depositional age from the LRa section (Hodgin et al., 2024).

Outputs: ``../Fuentes2025_Freda_<DD.Mon.YYYY>.txt`` (combined upload file) and the
``sites``, ``samples``, ``specimens``, ``locations`` and ``ages`` tables in ``..``.
Run from any directory: ``python build_upper_freda.py``.
"""
import io
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
POLE_DIR = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(POLE_DIR))
SOURCE = os.path.join(POLE_DIR, 'previous_MagIC', 'magic_contribution_20259.txt')
NORDIC = os.path.join(REPO, 'data', 'nordic_summaries', '1045_Upper_Freda.csv')
SEP = '>>>>>>>>>>'
CIT = 'This study'
HODGIN2024 = '10.1029/2024TC008415'
DIR_COLS = ['dir_dec', 'dir_inc']

# Fuentes et al. (2025), Table S1: section coordinates (latitude N, longitude E)
SECTIONS = {
    'BRa': (46.3850, -90.6336), 'BRb': (46.3864, -90.6372), 'BRc': (46.3914, -90.6429),
    'BRd': (46.4135, -90.6393), 'WRF': (46.4990, -90.9078), 'LRa': (46.5381, -90.3761),
    'LRb': (46.5660, -90.4149), 'WP': (46.9995, -89.0279),
}
LOCATION_NAMES = {
    'BRa': 'Bad River BRa', 'BRb': 'Bad River BRb', 'BRc': 'Bad River BRc',
    'BRd': 'Bad River BRd', 'WRF': 'White River Flowage', 'LRa': 'Lake Road LRa',
    'LRb': 'Lake Road LRb', 'WP': 'Wolf Point',
}
# depositional ages (Ma): model mean, 95% credible bounds
AGES = {
    'BR': ('1075.7', '1064.0', '1085.0'),
    'WRF': ('1061.3', '1041.6', '1082.7'),
    'LR': ('1046.1', '1033.4', '1052.9'),
    'WP': ('1045.0', '1031.7', '1063.2'),
}
AGE_TEXT = {
    'BR': ('Depositional age at 2200 m in the Oronto Group age model of Fuentes et al. (2025): '
           'model mean and 95% credible interval.'),
    'WRF': ('Depositional age at 3720 m in the Oronto Group age model of Fuentes et al. (2025), '
            'with bounds spanning the 95% credible intervals over the section\'s possible '
            'stratigraphic heights (2620-4820 m).'),
    'LR': ('Depositional age at 5800 m in the Oronto Group age model of Fuentes et al. (2025): '
           'model mean and 95% credible interval.'),
    'WP': ('Depositional age of the upper Freda pole from the Oronto Group age model of '
           'Fuentes et al. (2025, Table 1), with bounds spanning the 95% credible intervals over '
           'the section\'s possible stratigraphic heights (5000-6000 m).'),
}
STATE = {'BR': 'Wisconsin', 'WRF': 'Wisconsin', 'LR': '', 'WP': 'Michigan'}
SECTION_TEXT = {
    'BR': ('Steeply dipping (about 80 degrees) lower Freda Formation of the Montreal River '
           'Monocline along the Bad River, Wisconsin, a few hundred meters above the contact '
           'with the Nonesuch Formation.'),
    'WRF': ('Moderately dipping (about 20 degrees) middle Freda Formation within the Ashland '
            'Syncline near the White River Flowage, Wisconsin.'),
    'LR': ('Steeply dipping uppermost Freda Formation of the Montreal River Monocline near the '
           'mouth of the Montreal River on the Michigan-Wisconsin border.'),
    'WP': ('Nearly flat-lying (about 4 degrees) upper Freda Formation at Wolf Point, Michigan, '
           'outside the Montreal River Monocline and Ashland Syncline.'),
}
GEOLOGY = dict(geologic_classes='Sedimentary', geologic_types='Sediment Layer',
               lithologies='Sandstone')
COMP_TEXT = {'HT': 'high-temperature (HT) specularite detrital remanent magnetization',
             'MT': 'mid-temperature (MT) pigmentary hematite chemical remanent magnetization'}

# Fuentes et al. (2025), Table 1: SVEI-corrected Kent mean pole (authors' calculation,
# full precision from the study code output)
KENT = dict(lat=3.6194515351092487, lon=168.7837450095129,
            zdec=262.50564916575206, zinc=45.74126931524531, zeta=2.8009635559872477,
            edec=75.2779329704423, einc=44.03028658650698, eta=2.3990116901247744)

# published values checked by the build: (plat, plon, A95, n) and directions (dec, inc, a95, k)
EXPECT_HT_POLE = (4.7, 169.8, 2.1, 90)
EXPECT_HT_DIR = (281.5, -7.9, 2.8, 28.4)
EXPECT_MT_POLE = (3.3, 170.0, 2.2, 55)


def read_blocks(path):
    """Return {table: body} of a combined MagIC file; body is the raw text after the marker."""
    with open(path, encoding='utf-8-sig', newline='') as f:
        text = f.read()
    blocks = {}
    for chunk in text.split(SEP):
        chunk = chunk.strip('\r\n')
        if chunk:
            marker, _, body = chunk.partition('\n')
            blocks[marker.strip().split('\t')[-1]] = body
    return blocks


def to_df(body):
    return pd.read_csv(io.StringIO(body), sep='\t', dtype=str, keep_default_na=False)


def group_of(section):
    return section[:2] if section[:2] in ('BR', 'LR') else section


def fmt(x, nd=1):
    return f'{x:.{nd}f}'


def fix_specimens(body):
    """Link specimen WRF.76a to its own sample (the archived table links it to WRF.75a)."""
    lines = body.split('\n')
    header = lines[0].split('\t')
    i_spec, i_samp = header.index('specimen'), header.index('sample')
    n_fixed = 0
    for k, line in enumerate(lines[1:], start=1):
        cells = line.split('\t')
        if len(cells) > i_samp and cells[i_spec] == 'WRF.76a' and cells[i_samp] == 'WRF.75a':
            cells[i_samp] = 'WRF.76a'
            lines[k] = '\t'.join(cells)
            n_fixed += 1
    assert n_fixed == 1, n_fixed
    return '\n'.join(lines)


def build_sites(samples, specimens):
    """One site per sample: specimen directions per component and coordinate system."""
    rows = []
    for _, smp in samples.iterrows():
        section = smp['sample'].split('.')[0]
        lat, lon = SECTIONS[section]
        age, age_low, age_high = AGES[group_of(section)]
        base = dict(site=smp['sample'], location=LOCATION_NAMES[section], samples=smp['sample'],
                    citations=CIT, **GEOLOGY, lat=fmt(lat, 4), lon=fmt(lon % 360, 4),
                    height=smp['height'], age=age, age_low=age_low, age_high=age_high,
                    age_unit='Ma')
        spec = specimens[specimens['sample'] == smp['sample']]
        fits = spec[spec['dir_comp'] != '']
        if fits.empty:
            row = dict(base, specimens=':'.join(spec['specimen']))
            if spec.empty:
                row['method_codes'] = smp['method_codes']
                row['description'] = 'No specimen of this sample was measured.'
            else:
                row['method_codes'] = ':'.join(sorted(set(spec['method_codes'])))
                row['description'] = ('No interpretable remanence component was isolated from '
                                      'this sample; it is in no mean.')
            rows.append(row)
            continue
        for _, fit in fits.sort_values(['dir_comp', 'dir_tilt_correction'],
                                       key=lambda s: s.astype(str)).iterrows():
            row = dict(base, specimens=fit['specimen'], result_quality=fit['result_quality'],
                       method_codes=fit['method_codes'],
                       dir_tilt_correction=fit['dir_tilt_correction'],
                       dir_dec=fit['dir_dec'], dir_inc=fit['dir_inc'],
                       dir_comp_name=fit['dir_comp'], dir_n_samples='1', dir_n_specimens='1')
            if fit['dir_tilt_correction'] == '100':
                plon, plat = pmag.dia_vgp(float(fit['dir_dec']), float(fit['dir_inc']), 0,
                                          lat, lon % 360)[:2]
                row.update(vgp_lat=fmt(plat, 2), vgp_lon=fmt(plon, 2))
            if fit['result_quality'] == 'b':
                row['description'] = ('This specimen fit is flagged as unreliable and is in no '
                                      'mean.')
            rows.append(row)
    return pd.DataFrame(rows).fillna('')


def tc_sites(sites, section, comp):
    """Tilt-corrected, 'g' site rows of one component in one section."""
    return sites[(sites['location'] == LOCATION_NAMES[section]) & (sites['dir_comp_name'] == comp)
                 & (sites['dir_tilt_correction'] == '100') & (sites['result_quality'] == 'g')]


def dir_mean(rows):
    m = ipmag.fisher_mean(dec=rows['dir_dec'].astype(float).tolist(),
                          inc=rows['dir_inc'].astype(float).tolist())
    return dict(dir_dec=fmt(m['dec']), dir_inc=fmt(m['inc']), dir_alpha95=fmt(m['alpha95']),
                dir_k=fmt(m['k']), dir_n_sites=str(int(m['n'])))


def pole_mean(rows):
    m = ipmag.fisher_mean(rows['vgp_lon'].astype(float).tolist(),
                          rows['vgp_lat'].astype(float).tolist())
    return m


def location_base(section, site_names):
    lat, lon = SECTIONS[section]
    group = group_of(section)
    age, age_low, age_high = AGES[group]
    base = dict(location=LOCATION_NAMES[section], location_type='Stratigraphic Section',
                citations=CIT, geologic_classes='Sedimentary', lithologies='Sandstone',
                lat_s=fmt(lat, 4), lat_n=fmt(lat, 4),
                lon_w=fmt(lon % 360, 4), lon_e=fmt(lon % 360, 4),
                age=age, age_low=age_low, age_high=age_high, age_unit='Ma',
                continent_ocean='North America', country='United States of America',
                state_province=STATE[group], sites=':'.join(site_names))
    return base


def build_locations(sites):
    rows = []
    for section in SECTIONS:
        names = sorted(sites.loc[sites['location'] == LOCATION_NAMES[section], 'site'].unique(),
                       key=lambda s: (len(s), s))
        group = group_of(section)
        base = location_base(section, names)
        no_pole = (' The site rows of this section are in no pole.' if section != 'WP' else
                   ' The HT and MT site rows are in the Wolf Point poles; sites with no '
                   'interpretable component are in no pole.')
        rows.append(dict(base, description=(
            f'Freda Formation section {section} of Fuentes et al. (2025): '
            + SECTION_TEXT[group] + ' Each sample is from its own stratigraphic level and is '
            'treated as a site.' + no_pole + ' ' + AGE_TEXT[group])))
        if section == 'WP':
            continue
        for comp in ('HT', 'MT'):
            sel = tc_sites(sites, section, comp)
            row = dict(base, sites=':'.join(sel['site']), result_type='a',
                       result_name=f'{section} {comp} tilt-corrected section mean',
                       dir_tilt_correction='100', **dir_mean(sel),
                       method_codes='LP-DIR-T:DE-FM:DA-DIR-TILT')
            row['description'] = (
                f'Fisher mean of the tilt-corrected {COMP_TEXT[comp]} directions of section '
                f'{section} (one specimen per site; sites flagged b excluded). No pole is '
                'reported because the sections of the Montreal River Monocline and Ashland '
                'Syncline record vertical-axis rotation (Fuentes et al., 2025).')
            rows.append(row)
    return rows


def wolf_point_poles(sites):
    ht = tc_sites(sites, 'WP', 'HT')
    mt = tc_sites(sites, 'WP', 'MT')
    ht_dir, mt_dir = dir_mean(ht), dir_mean(mt)
    ht_pole, mt_pole = pole_mean(ht), pole_mean(mt)

    def round_pole(p):
        return (round(p['inc'], 1), round(p['dec'], 1), round(p['alpha95'], 1), int(p['n']))
    assert round_pole(ht_pole) == EXPECT_HT_POLE, round_pole(ht_pole)
    assert round_pole(mt_pole) == EXPECT_MT_POLE, round_pole(mt_pole)
    got_dir = tuple(float(ht_dir[k]) for k in ('dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k'))
    assert got_dir == EXPECT_HT_DIR, got_dir
    nordic = pd.read_csv(NORDIC, encoding='utf-8-sig').iloc[0]
    assert (nordic['PLAT'], nordic['PLONG'], nordic['A95'], nordic['N']) == EXPECT_HT_POLE, \
        nordic[['PLAT', 'PLONG', 'A95', 'N']]
    print(f'-I- poles reproduced: HT {round_pole(ht_pole)}, MT {round_pole(mt_pole)}')

    base = location_base('WP', [])
    base.update(age=AGES['WP'][0], result_type='a', result_quality='g', fold_test='F+',
                dir_tilt_correction='100')
    age_note = ' ' + AGE_TEXT['WP']
    fold_note = (' The HT and MT components pass a bootstrap fold test on the combined Freda '
                 'Formation data of this study and Henry et al. (1977) (tightest grouping at '
                 '95-102% and 97-104% unfolding).')

    measured = dict(base, sites=':'.join(ht['site']), **ht_dir,
                    result_name='Upper Freda Formation (Wolf Point) HT specularite DRM pole',
                    method_codes='LP-DIR-T:DE-FM:DA-DIR-TILT:DE-VGP:ST-F', pole_comp_name='HT',
                    pole_lat=fmt(ht_pole['inc']), pole_lon=fmt(ht_pole['dec']),
                    pole_alpha95=fmt(ht_pole['alpha95']), pole_k=fmt(ht_pole['k']),
                    pole_n_sites=str(int(ht_pole['n'])), pole_reversed_perc='0')
    measured['description'] = (
        'Paleomagnetic pole of the upper Freda Formation: Fisher mean of the virtual geomagnetic '
        'poles of the tilt-corrected HT specularite detrital remanent magnetization of the Wolf '
        'Point section, one specimen from each of 90 stratigraphic levels, with no '
        'inclination-shallowing correction (Fuentes et al., 2025, Table 1 note). The section is '
        'nearly flat-lying and outside the structures that record vertical-axis rotation. The HT '
        'specularite component of lower Freda siltstone intraclasts passes an intraformational '
        'conglomerate test (Swanson-Hysell, Fairchild and Slotznick, 2019). All directions are '
        'of one polarity.' + fold_note + age_note)

    kent = dict(base, sites=':'.join(ht['site']),
                result_name='Upper Freda Formation (Wolf Point) SVEI-corrected Kent mean pole',
                method_codes='LP-DIR-T:DE-FM:DA-DIR-TILT:DE-VGP:DA-DIR-FLAT:ST-F',
                pole_comp_name='HT', pole_lat=fmt(KENT['lat']), pole_lon=fmt(KENT['lon']),
                pole_alpha95=fmt((KENT['zeta'] * KENT['eta']) ** 0.5),
                pole_n_sites=str(int(ht_pole['n'])), pole_reversed_perc='0',
                pole_conf=':'.join(fmt(KENT[k]) for k in
                                   ('zdec', 'zinc', 'zeta', 'edec', 'einc', 'eta')))
    kent['description'] = (
        'Inclination-shallowing-corrected pole of the upper Freda Formation, the preferred pole '
        'of Fuentes et al. (2025, Table 1): the Wolf Point HT specimen directions corrected with '
        'flattening factors drawn from the SVEI 95% range of 0.59-0.87 (midpoint 0.73; Tauxe et '
        'al., 2024), with the flattening and directional uncertainty propagated following '
        'Pierce et al. (2022) and summarized as a Kent mean pole. pole_conf gives the Kent 95% '
        'ellipse (major axis declination, inclination and semi-angle, then minor axis); '
        'pole_alpha95 is the radius of the circle with the same area.' + fold_note + age_note)

    crm = dict(base, sites=':'.join(mt['site']), **mt_dir,
               result_name='Upper Freda Formation (Wolf Point) MT pigmentary hematite CRM pole',
               method_codes='LP-DIR-T:DE-FM:DA-DIR-TILT:DE-VGP:ST-F', pole_comp_name='MT',
               pole_lat=fmt(mt_pole['inc']), pole_lon=fmt(mt_pole['dec']),
               pole_alpha95=fmt(mt_pole['alpha95']), pole_k=fmt(mt_pole['k']),
               pole_n_sites=str(int(mt_pole['n'])), pole_reversed_perc='0')
    crm['description'] = (
        'Pole of the MT pigmentary hematite chemical remanent magnetization of the Wolf Point '
        'section: Fisher mean of the tilt-corrected site virtual geomagnetic poles (Fuentes et '
        'al., 2025, Table 1 note). The chemical remanence formed after deposition and before '
        'tilting during structural inversion (Hodgin et al., 2024); its direction is '
        'indistinguishable from the Wolf Point HT direction. The age is that of the rock.'
        + fold_note + age_note)
    return [measured, kent, crm]


def build_ages():
    return pd.DataFrame([dict(
        location=LOCATION_NAMES['LRa'], age='1051.58', age_sigma='0.555', age_unit='Ma',
        method_codes='GM-UPB', citations=HODGIN2024, timescale_eon='Proterozoic',
        timescale_era='Mesoproterozoic', timescale_period='Stenian',
        description=('Maximum depositional age of the upper Freda Formation: youngest CA-ID-TIMS '
                     'detrital zircon 206Pb/238U date from samples LRa-67.01 and LRa-38.39 of '
                     'the LRa section, 1051.58 +/- 1.11 Ma (2 sigma; Hodgin et al., 2024).'))])


def table_text(name, df):
    return f'tab delimited\t{name}\n' + df.to_csv(sep='\t', index=False)


def main():
    blocks = read_blocks(SOURCE)
    specimens_body = fix_specimens(blocks['specimens'])
    specimens = to_df(specimens_body)
    samples = to_df(blocks['samples'])
    samples['site'] = samples['sample']

    sites = build_sites(samples, specimens)
    loc_rows = build_locations(sites)
    wp_rows = wolf_point_poles(sites)
    # Wolf Point pole rows follow its row listing every site
    i_wp = next(i for i, r in enumerate(loc_rows) if r['location'] == LOCATION_NAMES['WP'])
    loc_rows[i_wp + 1:i_wp + 1] = wp_rows
    locations = pd.DataFrame(loc_rows).fillna('')
    ages = build_ages()

    site_cols = ['site', 'location', 'samples', 'specimens', 'result_quality', 'method_codes',
                 'citations', 'geologic_classes', 'geologic_types', 'lithologies', 'lat', 'lon',
                 'height', 'age', 'age_low', 'age_high', 'age_unit', 'dir_comp_name',
                 'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_n_samples', 'dir_n_specimens',
                 'vgp_lat', 'vgp_lon', 'description']
    sites = sites[site_cols]

    tables = [('locations', locations), ('sites', sites), ('samples', samples),
              ('specimens', specimens), ('ages', ages)]
    for name, df in tables:
        with open(os.path.join(POLE_DIR, f'{name}.txt'), 'w', encoding='utf-8', newline='') as f:
            f.write(table_text(name, df))

    stamp = date.today().strftime('%d.%b.%Y')
    upload = os.path.join(POLE_DIR, f'Fuentes2025_Freda_{stamp}.txt')
    parts = [table_text(name, df) for name, df in tables]
    measurements = blocks['measurements']
    parts.append('tab delimited\tmeasurements\n' + measurements
                 + ('' if measurements.endswith('\n') else '\n'))
    with open(upload, 'w', encoding='utf-8', newline='') as f:
        f.write((SEP + '\n').join(parts))
    print(f'-I- wrote {upload}: {len(locations)} location rows, {sites["site"].nunique()} sites '
          f'({len(sites)} rows), {len(samples)} samples, {len(specimens)} specimen rows, '
          f'{len(ages)} ages; measurements copied from 20259')

    sys.path.insert(0, os.path.join(REPO, 'scripts'))
    from validate_magic_contribution import validate_upload_file
    if not validate_upload_file(upload, server=True):
        sys.exit('-E- validation failed')


if __name__ == '__main__':
    main()
