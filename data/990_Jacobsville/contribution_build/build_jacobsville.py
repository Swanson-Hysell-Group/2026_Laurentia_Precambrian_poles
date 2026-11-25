"""Build the updated MagIC contribution for the Jacobsville Formation data of Zhang et al. (2024).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Zhang, Y., Hodgin, E. B., Alemu, T., Pierce, J., Fuentes, A., & Swanson-Hysell, N. L. (2024).
Tracking Rodinia into the Neoproterozoic: new paleomagnetic constraints from the Jacobsville
Formation. Tectonics 43, e2023TC007866. doi:10.1029/2023TC007866

Ages: Hodgin, E. B., Swanson-Hysell, N. L., et al. (2022). Final inversion of the
Midcontinent Rift during the Rigolet phase of the Grenvillian orogeny. Geology 50, 547-551.
doi:10.1130/G49439.1

The published contribution earthref.org/MagIC/19780 holds the sites (each a single specimen
component), samples, specimens and measurements of the study and four location-level rows.
This build reads the archived copy ``../previous_MagIC/magic_contribution_19780.txt``, copies
the samples, specimens and measurements tables through unchanged, and revises the sites and
locations tables and adds an ages table.

Sites
-----
- Every site of 19780 carries one coordinate (47.3 N, 271.8 E, the Snake Creek tributary
  section rounded). Each site is given the coordinate of its section, as in the study's
  section-location table (``SECTIONS`` below; the GPS positions of the section tops and bases
  are printed on Fig. 2 of the paper). The VGP of each site row is recalculated from its
  direction at that coordinate; vgp_dp and vgp_dm are left blank (single specimens, no
  alpha95).
- Location 'Jacobsville Formation' (the sections span 46.48-47.34 N; Agate Falls lies south of
  the Keweenaw Peninsula). Site ages are the depositional bracket 980-1000 Ma.
- Section codes: AF Agate Falls (AFC intraclast conglomerate), HF Hungarian side falls on Dover
  Creek (HFC intraclast conglomerate), JK Hammel Creek, NW Natural Wall, SC Snake Creek
  tributary.

Locations (one location, four rows)
-----------------------------------
1. A row without a pole listing every site.
2. Detrital remanence (hdt) pole as measured: the 307 tilt-corrected hdt specimen directions
   of Zhang et al. (2024) from AF1, AF2, HF1-HF6 (HF3-2 not used), JK2, JK3, NW2 and SC1.
   The mean direction unifies polarity with the authors' grouping (AF1, AF2, HF1 and HF5
   split by angle to each section's principal direction; HF3, NW2 and SC1 normal; HF2, HF4,
   HF6, JK2 and JK3 reversed); the pole is the Fisher mean of the specimen VGPs with those of
   specimens with declination below 180 taken to the antipode, as in the study's code. This
   reproduces the published pole (185.5 E, 14.0 S, A95 2.7) and is the pole of the
   Swanson-Hysell et al. (2026) compilation. conglomerate_test G+ and fold_test F+.
3. E/I-corrected pole: the grouped directions with inclinations unflattened at f = 0.64
   (95% range 0.51-0.83) give the corrected mean direction and a Fisher mean pole of
   183.4 E, 16.9 S (A95 3.1); the row carries the Kent mean pole and ellipse of Table 1
   (uncertainty propagated following Pierce et al., 2022) in pole_conf, with pole_alpha95
   the radius of the circle of equal area, sqrt(zeta95 * eta95).
4. Pigmentary hematite (hct) chemical remanence pole: the 186 geographic-coordinate hct
   specimen directions of AF1, AF2, AFC, HF1-HF4 (HF3-2 not used), HF6, HFC, JK2, JK3 and SC1,
   polarity assigned by angle to the principal direction, reproducing the published pole
   (179.6 E, 10.2 S, A95 3.6).

Ages table: the youngest CA-ID-TIMS detrital zircon dates at Agate Falls and Sandstone Creek
(maximum depositional ages) and the U-Pb date of late-kinematic calcite in the Keweenaw fault
(Hodgin et al., 2022).

Outputs: ``../sites.txt``, ``../locations.txt``, ``../ages.txt`` and
``../Zhang2024_Jacobsville_<DD.Mon.YYYY>.txt`` (combined upload file).
Run from any directory: ``python build_jacobsville.py``.
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
SOURCE = os.path.join(POLE_DIR, 'previous_MagIC', 'magic_contribution_19780.txt')
SEP = '>>>>>>>>>>'
PASSTHROUGH = ['samples', 'specimens', 'measurements']
CIT = 'This study'
HODGIN2022 = '10.1130/G49439.1'
LOCATION = 'Jacobsville Formation'
AGE = dict(age='990', age_low='980', age_high='1000', age_unit='Ma')
GEOG = dict(continent_ocean='North America', country='United States of America',
            state_province='Michigan')

# section coordinates (lat N, lon E) of the study's section-location table
SECTIONS = {
    'AF1': (46.48098, -89.09036), 'AF2': (46.48116, -89.09022), 'AFC': (46.48116, -89.09022),
    'SC1': (47.34248, -88.17128), 'NW1': (47.24099, -88.39637), 'NW2': (47.2402, -88.39146),
    'JK2': (47.20312, -88.41193), 'JK3': (47.20463, -88.41516), 'HF1': (47.17071, -88.4465),
    'HF2': (47.17034, -88.44574), 'HF3': (47.17047, -88.44443), 'HF4': (47.17019, -88.44321),
    'HF5': (47.17019, -88.44321), 'HF6': (47.16961, -88.44215), 'HF7': (47.16919, -88.44072),
    'HFC': (47.17071, -88.4465)}
SECTION_NAMES = {'AF': 'Agate Falls', 'HF': 'Hungarian side falls, Dover Creek',
                 'JK': 'Hammel Creek', 'NW': 'Natural Wall', 'SC': 'Snake Creek tributary'}

# detrital remanence: sections and polarity grouping of the study's code
DRM_SPLIT = ['AF1', 'AF2', 'HF1', 'HF5']           # by angle to the section's principal direction
DRM_NORMAL = ['HF3', 'NW2', 'SC1']
DRM_REVERSED = ['HF2', 'HF4', 'HF6', 'JK2', 'JK3']
CRM_SECTIONS = ['AF1', 'AF2', 'AFC', 'HF1', 'HF2', 'HF3', 'HF4', 'HF6', 'HFC', 'JK2', 'JK3', 'SC1']
EXCLUDED_SPECIMEN = 'HF3-2'
F_EI = 0.64

# Zhang et al. (2024), Table 1: E/I-corrected Kent mean pole of the detrital remanence
KENT = dict(lat=-16.9, lon=183.4, zeta=4.1, eta=3.1)

# published values reproduced (Table 1 note; 19780)
EXPECT = {'drm_pole': (185.5, -14.0, 2.7, 307), 'drm_dir': (256.9, -13.4, 3.3),
          'ei_pole': (183.4, -16.9, 3.1), 'ei_dir': (256.7, -18.6, 3.9),
          'crm_pole': (179.6, -10.2, 3.6, 186), 'crm_dir': (263.8, -15.1, 4.2)}


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


def section_of(site):
    return site.split('-')[0]


def vgp(dec, inc, section):
    lat, lon = SECTIONS[section]
    plon, plat, _, _ = pmag.dia_vgp(dec, inc, 0, lat, lon % 360)
    return plon, plat


def polarity_by_principal(block):
    """'N' for directions more than 90 deg from the principal direction, else 'R'."""
    pp = pmag.doprinc(np.asarray(block))
    return ['N' if pmag.angle(list(d), [pp['dec'], pp['inc']])[0] > 90 else 'R' for d in block]


def specimen_directions(specimens, samples, comp, tilt):
    sp = specimens.merge(samples[['sample', 'site']], on='sample', how='left')
    tc = pd.to_numeric(sp['dir_tilt_correction'], errors='coerce')
    out = sp[(sp['dir_comp'] == comp) & (tc == tilt) & (sp['specimen'] != EXCLUDED_SPECIMEN)].copy()
    out['section'] = out['site'].map(section_of)
    out['dec'] = pd.to_numeric(out['dir_dec'])
    out['inc'] = pd.to_numeric(out['dir_inc'])
    return out


def rounded(m, keys, nd=1):
    return tuple(round(float(m[k]), nd) for k in keys)


def drm_results(specimens, samples):
    h = specimen_directions(specimens, samples, 'hdt', 100)
    h = h[h['section'].isin(DRM_SPLIT + DRM_NORMAL + DRM_REVERSED)].copy()
    pol = pd.Series('N', index=h.index)
    for sec in DRM_SPLIT:
        sub = h[h['section'] == sec]
        pol[sub.index] = polarity_by_principal(sub[['dec', 'inc']].to_numpy())
    pol[h['section'].isin(DRM_REVERSED)] = 'R'
    rev = pol == 'R'
    D = np.where(rev, (h['dec'] + 180) % 360, h['dec'])
    I = np.where(rev, -h['inc'], h['inc'])
    direction = ipmag.fisher_mean(D.tolist(), I.tolist())
    # pole: VGPs of specimens with declination below 180 taken to the antipode
    vgps = []
    for (_, r) in h.iterrows():
        plon, plat = vgp(r['dec'], r['inc'], r['section'])
        vgps.append([(plon + 180) % 360, -plat] if r['dec'] < 180 else [plon, plat])
    pole = ipmag.fisher_mean(di_block=vgps)
    # E/I correction of the grouped directions
    I_cor = np.array(ipmag.unsquish(I.tolist(), F_EI))
    dir_cor = ipmag.fisher_mean(D.tolist(), I_cor.tolist())
    pole_cor = ipmag.fisher_mean(di_block=[list(vgp(d, i, s)) for d, i, s in
                                           zip(D, I_cor, h['section'])])
    return h, direction, pole, dir_cor, pole_cor


def crm_results(specimens, samples):
    h = specimen_directions(specimens, samples, 'hct', 0)
    h = h[h['section'].isin(CRM_SECTIONS)].copy()
    block = h[['dec', 'inc']].to_numpy()
    pol = polarity_by_principal(block)
    vgps = []
    for (_, r), p in zip(h.iterrows(), pol):
        plon, plat = vgp(r['dec'], r['inc'], r['section'])
        vgps.append([plon, plat] if p == 'N' else [(plon + 180) % 360, -plat])
    pole = ipmag.fisher_mean(di_block=vgps)
    direction = ipmag.fisher_mean(di_block=ipmag.do_flip(di_block=pmag.flip(block.tolist(),
                                                                           combine=True)))
    return h, direction, pole


def revise_sites(sites):
    s = sites.copy()
    sec = s['site'].map(section_of)
    assert set(sec) == set(SECTIONS), sorted(set(sec) ^ set(SECTIONS))
    s['location'] = LOCATION
    s['lat'] = [f'{SECTIONS[x][0]:.5f}' for x in sec]
    s['lon'] = [f'{SECTIONS[x][1] % 360:.5f}' for x in sec]
    s['age_low'], s['age_high'], s['age_unit'] = AGE['age_low'], AGE['age_high'], 'Ma'
    vlat, vlon = [], []
    for (_, r), x in zip(s.iterrows(), sec):
        if r['dir_dec'] == '':
            vlat.append('')
            vlon.append('')
            continue
        plon, plat = vgp(float(r['dir_dec']), float(r['dir_inc']), x)
        vlat.append(f'{plat:.1f}')
        vlon.append(f'{plon:.1f}')
    s['vgp_lat'], s['vgp_lon'] = vlat, vlon
    s['vgp_dp'], s['vgp_dm'] = '', ''
    return s


def check(name, got, expected):
    assert got == expected, (name, got, expected)


def main():
    blocks = read_blocks(SOURCE)
    specimens = to_df(blocks['specimens'])
    samples = to_df(blocks['samples'])
    sites = to_df(blocks['sites'])
    published = to_df(blocks['locations'])
    assert len(published) == 4

    h, d, p, d_cor, p_cor = drm_results(specimens, samples)
    check('drm_pole', rounded(p, ('dec', 'inc', 'alpha95')) + (int(p['n']),), EXPECT['drm_pole'])
    check('drm_dir', rounded(d, ('dec', 'inc', 'alpha95')), EXPECT['drm_dir'])
    check('ei_pole', rounded(p_cor, ('dec', 'inc', 'alpha95')), EXPECT['ei_pole'])
    check('ei_dir', rounded(d_cor, ('dec', 'inc', 'alpha95')), EXPECT['ei_dir'])
    assert set(h['site']) == set(published.iloc[0]['sites'].split(':'))
    c, cd, cp = crm_results(specimens, samples)
    check('crm_pole', rounded(cp, ('dec', 'inc', 'alpha95')) + (int(cp['n']),), EXPECT['crm_pole'])
    check('crm_dir', rounded(cd, ('dec', 'inc', 'alpha95')), EXPECT['crm_dir'])
    print(f'-I- published poles reproduced: DRM {EXPECT["drm_pole"]}, E/I {EXPECT["ei_pole"]}, '
          f'CRM {EXPECT["crm_pole"]}')

    new_sites = revise_sites(sites)
    all_sites = list(dict.fromkeys(new_sites['site']))
    lats = pd.to_numeric(new_sites['lat'])
    lons = pd.to_numeric(new_sites['lon'])
    bbox = dict(lat_s=f'{lats.min():.5f}', lat_n=f'{lats.max():.5f}',
                lon_w=f'{lons.min():.5f}', lon_e=f'{lons.max():.5f}')
    common = dict(location=LOCATION, location_type='Region', citations=f'{CIT}:{HODGIN2022}',
                  geologic_classes='Sedimentary', lithologies='Siltstone:Sandstone',
                  result_type='a', result_quality='g', **AGE, **GEOG)
    age_text = ('Age: deposition after the youngest detrital zircon dates of 1003.21 +/- 2.23 Ma '
                '(Agate Falls) and 992.51 +/- 0.64 Ma (Sandstone Creek) and before late-kinematic '
                'calcite in the Keweenaw fault (985.5 +/- 35.8 Ma) and the ca. 980 Ma end of the '
                'Rigolet phase of the Grenvillian orogeny (Hodgin et al., 2022; Zhang et al., '
                '2024).')

    def fmt(x, nd=1):
        return f'{float(x):.{nd}f}'

    def lists(df):
        names = list(dict.fromkeys(df['site']))
        return dict(sites=':'.join(names), samples=':'.join(dict.fromkeys(df['sample'])),
                    specimens=':'.join(dict.fromkeys(df['specimen'])))

    no_pole = dict(
        common, **bbox, sites=':'.join(all_sites), method_codes='LP-DIR-T:DE-BFL:DE-FM',
        description=(
            'Jacobsville Formation, northern Michigan: all sites of Zhang et al. (2024), each a '
            'single specimen component, from the Agate Falls (AF), Hungarian side falls on '
            'Dover Creek (HF), Hammel Creek (JK), Natural Wall (NW) and Snake Creek tributary '
            '(SC) sections. Sites in no pole: the siltstone intraclasts of the conglomerate '
            'tests (AFC, HFC) except in the chemical remanence pole, sections HF7 and NW1, the '
            'low-temperature components, specimen HF3-2 and specimens without a fit. Site '
            'coordinates are those of the sections. ' + age_text))

    drm_bbox_sites = new_sites[new_sites['site'].isin(h['site'])]

    def bbox_of(df):
        la, lo = pd.to_numeric(df['lat']), pd.to_numeric(df['lon'])
        return dict(lat_s=f'{la.min():.5f}', lat_n=f'{la.max():.5f}',
                    lon_w=f'{lo.min():.5f}', lon_e=f'{lo.max():.5f}')

    drm = dict(
        common, **bbox_of(drm_bbox_sites), **lists(h),
        result_name='Jacobsville Formation detrital remanence (hdt) pole',
        method_codes='LP-DIR-T:DE-BFL:DE-FM:DE-VGP:ST-G:ST-F', pole_comp_name='hdt',
        dir_tilt_correction='100', dir_dec=fmt(d['dec']), dir_inc=fmt(d['inc']),
        dir_alpha95=fmt(d['alpha95']), dir_k=fmt(d['k']), dir_n_sites=str(int(d['n'])),
        pole_lat=fmt(p['inc']), pole_lon=fmt(p['dec']), pole_alpha95=fmt(p['alpha95']),
        pole_k=fmt(p['k']), pole_n_sites=str(int(p['n'])),
        conglomerate_test='G+', fold_test='F+',
        description=(
            'Paleomagnetic pole of the Jacobsville Formation detrital remanent magnetization, '
            'without inclination-shallowing correction (Zhang et al., 2024, Table 1 note): the '
            'Fisher mean of the virtual geomagnetic poles of 307 tilt-corrected high-temperature '
            'detrital hematite (hdt) specimen directions from sections AF1, AF2, HF1-HF6, JK2, '
            'JK3, NW2 and SC1, each specimen treated as an independent reading of the field, '
            'calculated at the section coordinates and reported as the southern pole. The mean '
            'direction is given in reversed polarity. The detrital remanence passes '
            'intraformational conglomerate tests on siltstone intraclasts at Agate Falls and '
            'Dover Creek and a bootstrap fold test at the Snake Creek tributary; directions are '
            'of both polarities, and the reversal test is negative (McFadden and McElhinny, '
            '1990). ' + age_text))

    area_a95 = (KENT['zeta'] * KENT['eta']) ** 0.5
    kent_conf = published.iloc[3]['pole_conf']
    assert kent_conf.startswith('255.53:45.17:4.14')
    ei = dict(
        common, **bbox_of(drm_bbox_sites), **lists(h),
        result_name='Jacobsville Formation detrital remanence (hdt) E/I-corrected pole',
        method_codes='LP-DIR-T:DE-BFL:DE-FM:DE-VGP:DA-DIR-FLAT-EI:ST-G:ST-F',
        pole_comp_name='hdt', dir_tilt_correction='100', dir_dec=fmt(d_cor['dec']),
        dir_inc=fmt(d_cor['inc']), dir_alpha95=fmt(d_cor['alpha95']), dir_k=fmt(d_cor['k']),
        dir_n_sites=str(int(d_cor['n'])), pole_lat=fmt(KENT['lat']), pole_lon=fmt(KENT['lon']),
        pole_alpha95=f'{area_a95:.1f}', pole_conf=kent_conf, pole_n_sites=str(int(p_cor['n'])),
        conglomerate_test='G+', fold_test='F+',
        description=(
            'Inclination-shallowing-corrected paleomagnetic pole of the Jacobsville Formation '
            'detrital remanence (Zhang et al., 2024, Table 1): the 307 tilt-corrected detrital '
            'hematite directions corrected with the elongation/inclination method (Tauxe and '
            'Kent, 2004), flattening factor f = 0.64 with a 95% range of 0.51-0.83, with the '
            'flattening uncertainty propagated to the pole following Pierce et al. (2022) and '
            'summarized as a Kent mean pole. The direction is the corrected mean direction at '
            'f = 0.64 (reversed polarity); the Fisher mean of the corrected VGPs at that f is '
            f'{abs(p_cor["inc"]):.1f}S {p_cor["dec"]:.1f}E with A95 {p_cor["alpha95"]:.1f}. '
            'pole_conf gives the Kent 95% ellipse (major axis declination, inclination and '
            'semi-angle, then minor axis); pole_alpha95 is the radius of the circle with the '
            'same area as that ellipse. ' + age_text))

    crm = dict(
        common, **bbox_of(new_sites[new_sites['site'].isin(c['site'])]), **lists(c),
        result_name='Jacobsville Formation pigmentary hematite (hct) chemical remanence pole',
        method_codes='LP-DIR-T:DE-BFL:DE-FM:DE-VGP', pole_comp_name='hct',
        dir_tilt_correction='0', dir_dec=fmt(cd['dec']), dir_inc=fmt(cd['inc']),
        dir_alpha95=fmt(cd['alpha95']), dir_k=fmt(cd['k']), dir_n_sites=str(int(cd['n'])),
        pole_lat=fmt(cp['inc']), pole_lon=fmt(cp['dec']), pole_alpha95=fmt(cp['alpha95']),
        pole_k=fmt(cp['k']), pole_n_sites=str(int(cp['n'])),
        description=(
            'Pole of the chemical remanent magnetization carried by pigmentary hematite (hct) '
            '(Zhang et al., 2024, Table 1 note): the Fisher mean of the VGPs of 186 specimen '
            'directions in geographic coordinates from sections AF1, AF2, HF1-HF4, HF6, JK2, JK3 '
            'and SC1 and the intraclasts AFC and HFC, calculated at the section coordinates and '
            'reported as the southern pole. The component fails the intraformational '
            'conglomerate and fold tests and postdates tilting; its similarity to the detrital '
            'remanence indicates growth of the pigmentary hematite soon after deposition '
            '(Zhang et al., 2024). The age given is that of deposition. ' + age_text))

    locs = pd.DataFrame([no_pole, drm, ei, crm]).fillna('')

    ages = pd.DataFrame([
        dict(location=LOCATION, age='1003.21', age_sigma='1.115', age_unit='Ma',
             method_codes='GM-UPB', citations=HODGIN2022, timescale_eon='Proterozoic',
             timescale_era='Mesoproterozoic', timescale_period='Stenian',
             description=('Maximum depositional age: youngest CA-ID-TIMS detrital zircon date of '
                          'sample AF1-29.3, Agate Falls section (46.48046N, 89.09026W), '
                          '1003.21 +/- 2.23 Ma (2 sigma; Hodgin et al., 2022).')),
        dict(location=LOCATION, age='992.51', age_sigma='0.32', age_unit='Ma',
             method_codes='GM-UPB', citations=HODGIN2022, timescale_eon='Proterozoic',
             timescale_era='Neoproterozoic', timescale_period='Tonian',
             description=('Maximum depositional age: youngest CA-ID-TIMS detrital zircon date of '
                          'sample SC18-1, Sandstone Creek (46.69132N, 89.16658W), 992.51 +/- '
                          '0.64 Ma (2 sigma; Hodgin et al., 2022).')),
        dict(location=LOCATION, age='985.5', age_sigma='17.9', age_unit='Ma',
             method_codes='GM-UPB', citations=HODGIN2022, timescale_eon='Proterozoic',
             timescale_era='Neoproterozoic', timescale_period='Tonian',
             description=('Minimum age constraint: LA-ICP-MS U-Pb lower-intercept date of '
                          'late-kinematic vein calcite in the brecciated Keweenaw fault (sample '
                          'LLB2, 47.38929N, 87.97621W), 985.5 +/- 35.8 Ma (2 sigma; Hodgin et '
                          'al., 2022); the Jacobsville Formation is folded in the footwall of '
                          'the fault.')),
    ])

    def table(kind, df):
        return f'tab delimited\t{kind}\n' + df.to_csv(sep='\t', index=False)

    for kind, df in (('sites', new_sites), ('locations', locs), ('ages', ages)):
        with open(os.path.join(POLE_DIR, f'{kind}.txt'), 'w', encoding='utf-8', newline='') as f:
            f.write(table(kind, df))

    stamp = date.today().strftime('%d.%b.%Y')
    upload = os.path.join(POLE_DIR, f'Zhang2024_Jacobsville_{stamp}.txt')
    parts = [table('locations', locs), table('sites', new_sites)]
    for kind in PASSTHROUGH:
        body = blocks[kind]
        parts.append(f'tab delimited\t{kind}\n' + body + ('' if body.endswith('\n') else '\n'))
    parts.append(table('ages', ages))
    with open(upload, 'w', encoding='utf-8', newline='') as f:
        f.write((SEP + '\n').join(parts))
    print(f'-I- wrote {upload} ({len(locs)} locations rows, {len(new_sites)} sites rows; '
          f'{", ".join(PASSTHROUGH)} copied from 19780)')

    sys.path.insert(0, os.path.join(REPO, 'scripts'))
    from validate_magic_contribution import validate_upload_file
    if not validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True):
        sys.exit('-E- validation failed')


if __name__ == '__main__':
    main()
