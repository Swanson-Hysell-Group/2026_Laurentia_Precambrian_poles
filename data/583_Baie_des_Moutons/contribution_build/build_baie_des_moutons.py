"""Build the MagIC contribution for the Baie des Moutons syenite complex (McCausland et al., 2011).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

McCausland, P. J. A., Hankard, F., Van der Voo, R., & Hall, C. M. (2011). Ediacaran
paleogeography of Laurentia: Paleomagnetism and 40Ar-39Ar geochronology of the 583 Ma Baie
des Moutons syenite, Quebec. Precambrian Research 187, 58-78.
doi:10.1016/j.precamres.2011.02.004

Only site means are published (the online appendix holds the 40Ar-39Ar step-heating data),
so this is a site-level contribution built from Table 1 (p. 61; site means, counts and
coordinates), Table 2 (p. 67; 40Ar-39Ar plateau ages) and sections 2-5 of the paper.

Sources and checks
------------------
- Table 1 was transcribed twice, independently, and checked against the rendered page;
  the two transcriptions agree for all 44 sites (rock type, unit, susceptibility, NRM,
  counts, directions, k, alpha95, ChRM type, coordinates and footnotes).
- Counts are specimen counts: N (specimens demagnetized) -> dir_n_total_specimens; e+g
  (endpoint and great-circle specimens used) -> dir_n_specimens, dir_n_specimens_lines,
  dir_n_specimens_planes. The number of samples per site mean is not reported, so
  dir_n_samples is left blank. The N column sums to 515; the text reports 530 specimens
  prepared.
- MB15 is a great-circle-only site mean (0+5) for which Table 1 gives k but no alpha95;
  its alpha95 and VGP dp/dm are left blank. MB27 carries the direction of its single
  ChRM-bearing specimen (footnote b), not a site mean.
- Site coordinates are NAD-27 decimal minutes added to 50 N and 58 W (Table 1), not
  converted; site VGPs are calculated from each site direction at its coordinates.
- The ChRM A mean of Table 1 (alpha95 6.6) is reproduced; the abstract and section 3.4
  give alpha95 6.5. The six-site ChRM B mean is reproduced exactly.
- Ages: the individual 40Ar/39Ar plateau ages are 1 sigma (section 4.1 and the Table 2
  heading). The 583.4 +/- 2.0 Ma weighted mean of the four hornblende ages is shown on
  Fig. 1 as 583 +/- 2 Ma, whose caption quotes age errors at 2 sigma (the 1-sigma
  standard error of the four 1-sigma ages is 0.85 Ma), so its age_sigma is 1.0.

Structure and interpretation
----------------------------
- One location, 'Baie des Moutons syenite complex': an umbrella row holding all 44 sites
  of Table 1, and two pole rows, one per characteristic component (pole_comp_name 'A' and
  'B', following the 'Type' column of Table 1). The two components are never observed
  together in one site, so the A and B site sets are disjoint. Sites without a site mean
  are entered with result_quality 'b' and the reason from Table 1.
- Poles: Fisher means of the in-situ site VGPs, polarity unified (MB20, of the polarity
  opposite to the other ChRM B sites, is inverted), reported in the hemispheres of the
  paper (A north, B south). Both are poles of the Swanson-Hysell et al. (2026)
  compilation ('Baie des Moutons complex A' and 'B'). The paper's poles BMa and BMb are
  poles of the mean directions (dp/dm). The four-dyke ChRM B pole and the five-site
  ChRM B mean of the paper are subsets of the BMb sites and are described, not entered.
- No tilt correction: the steep dykes of varied strike indicate that the complex has not
  been significantly tilted since emplacement (section 2). No structured field-test
  column: the six baked-contact tests were inconclusive and there is no reversal test.
- Ages are ages of the rock: every site and pole row carries the complex age, 583.4 Ma
  with the 2-sigma range 581.4-585.4 Ma. The ages table holds the plateau ages (Table 2)
  and the weighted mean.
"""
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))

CIT = '10.1016/j.precamres.2011.02.004'
LOCATION = 'Baie des Moutons syenite complex'
AGE = dict(age=583.4, age_low=581.4, age_high=585.4, age_unit='Ma')  # 583.4 +/- 2.0 (2 sigma)
GEOG = dict(continent_ocean='North America', country='Canada', state_province='Quebec',
            region='Quebec Lower North Shore, Gulf of St. Lawrence',
            terranes='Laurentia:Grenville Province', tectonic_settings='Continental Rift')

# Table 1 (p. 61). site, rock type (as printed), intrusive unit, susceptibility (1e-6 SI),
# NRM (mA/m), N (specimens demagnetized), e, g (endpoint / great-circle specimens used),
# Dec, Inc, alpha95, k, ChRM type, lat minutes (+50 N), lon minutes (+58 W), footnotes.
T1 = [
    ('MB01', 'Fsp-por dyke', 'd', 2819, 1000, 21, None, None, None, None, None, None, None, 46.68, 59.99, 'a'),
    ('MB02', 'Syenite xeno', '3', 1250, 151, 9, None, None, None, None, None, None, None, 52.13, 59.91, ''),
    ('MB03', 'Fol grey syenite', '2', 1862, 196, 13, 10, 1, 109.2, 80.8, 7.2, 41.5, 'A', 48.44, 58.62, ''),
    ('MB04', 'Fol grey syenite', '2', 1951, 211, 13, 7, 2, 82.9, 79.4, 7.1, 54.7, 'A', 48.65, 58.61, ''),
    ('MB05', 'Pink fsp syenite', '1', 1917, 3240, 8, 5, 0, 141.5, 82.8, 9.3, 68.8, 'A', 48.68, 58.63, ''),
    ('MB06', 'Aplite dyke', '1', 31, 3, 9, 7, 2, 183.8, -5.9, 5.9, 80.1, 'B', 48.42, 59.49, ''),
    ('MB07', 'Fsp-por dyke', 'd', 829, 80, 10, None, None, None, None, None, None, None, 48.48, 59.45, ''),
    ('MB08', 'Green syenite', '1', 1577, 251, 14, None, None, None, None, None, None, None, 48.37, 59.88, 'a'),
    ('MB09', 'Grey syenite', '2', 2177, 169, 9, None, None, None, None, None, None, None, 47.39, 59.97, ''),
    ('MB10', 'Aplite dyke', '2', 1661, 335, 12, None, None, None, None, None, None, None, 47.21, 59.98, ''),
    ('MB11', 'Red fsp dyke', 'd', 43, 5, 9, 4, 1, 162.8, -3.7, 4.4, 325.0, 'B', 47.83, 59.35, ''),
    ('MB12', 'Red fsp dyke', 'd', 20, 4, 9, 3, 2, 163.2, 23.7, 9.1, 92.1, 'B', 47.83, 59.35, ''),
    ('MB13', 'Syenite dyke', '2', 3634, 413, 14, None, None, None, None, None, None, None, 47.83, 59.35, 'c'),
    ('MB14', 'Syenite dyke', '2', 2806, 265, 18, None, None, None, None, None, None, None, 47.83, 59.35, 'c'),
    ('MB15', 'Grey syenite', '2', 193, 19, 19, 0, 5, 94.5, 64.4, None, 10.8, 'A', 47.75, 59.62, 'c'),
    ('MB16', 'Red fsp dyke', 'd', 102, 2, 11, None, None, None, None, None, None, None, 48.09, 59.35, 'ac'),
    ('MB17', 'Aplite dyke', '1', 5515, 1630, 10, None, None, None, None, None, None, None, 49.29, 58.30, ''),
    ('MB18', 'Red syenite', '3', 218, 417, 10, None, None, None, None, None, None, None, 49.38, 58.05, 'a'),
    ('MB19', 'Aplite dyke', '1', 2881, 404, 8, 5, 0, 46.4, 73.6, 11.4, 46.1, 'A', 51.64, 59.67, ''),
    ('MB20', 'Red syenite', '3', 703, 121, 10, 9, 0, 321.1, 22.0, 6.7, 59.6, '-B', 50.20, 58.56, ''),
    ('MB21', 'Aplite dyke', '2', 5495, 403, 10, None, None, None, None, None, None, None, 54.96, 57.01, ''),
    ('MB22', 'Aplite dyke', '2', 2652, 1430, 7, None, None, None, None, None, None, None, 54.96, 57.01, 'c'),
    ('MB23', 'Finer granite', '2', 4224, 562, 9, None, None, None, None, None, None, None, 55.14, 56.34, ''),
    ('MB24', 'Finer granite', '2', 3600, 504, 9, None, None, None, None, None, None, None, 55.14, 56.34, ''),
    ('MB25', 'Aplite + host', '2', 514, 91, 8, 5, 0, 78.7, 79.7, 9.7, 62.9, 'A', 54.61, 57.53, ''),
    ('MB26', 'Mafic dyke', '3', 3360, 258, 5, 5, 0, 136.0, 80.7, 8.5, 81.5, 'A', 54.61, 57.53, 'c'),
    ('MB27', 'Green syenite', '3', 625, 294, 8, 1, 0, 80.1, 57.4, None, None, None, 52.10, 58.88, 'b'),
    ('MB28', 'Xeno red syenite', '3', 3343, 390, 12, None, None, None, None, None, None, None, 50.12, 58.44, ''),
    ('MB29', 'Fol red syenite', '3', 1388, 303, 12, None, None, None, None, None, None, None, 48.97, 57.45, ''),
    ('MB30', 'Xeno in syenite', '3', 1940, 519, 15, None, None, None, None, None, None, None, 48.99, 57.50, ''),
    ('MB31', 'Fol grey syenite', '2', 2918, 1130, 12, None, None, None, None, None, None, None, 48.65, 51.78, ''),
    ('MB32', 'Fol grey syenite', '2', 2215, 396, 12, None, None, None, None, None, None, None, 48.62, 51.88, 'a'),
    ('MB33', 'Pink fs syenite', '1', 2267, 291, 12, None, None, None, None, None, None, None, 48.80, 51.05, ''),
    ('MB34', 'Aplite dyke', '1', 97, 192, 17, None, None, None, None, None, None, None, 48.79, 51.06, ''),
    ('MB35', 'Aplite dyke', '1', 43, 58, 16, None, None, None, None, None, None, None, 50.05, 56.01, ''),
    ('MB36', 'Red fsp dyke', 'd', 77, 4, 10, 9, 0, 173.1, 9.7, 10.3, 25.8, 'B', 48.09, 59.17, 'c'),
    ('MB37', 'Red fsp dyke', 'd', 209, 18, 13, None, None, None, None, None, None, None, 48.00, 58.97, ''),
    ('MB38', 'Fol green syenite', '1', 160, 28, 14, None, None, None, None, None, None, None, 50.00, 57.02, ''),
    ('MB39', 'Fol green syenite', '1', 26, 2, 12, None, None, None, None, None, None, None, 49.98, 57.00, ''),
    ('MB40', 'Green syenite', '1', 44, 3, 18, None, None, None, None, None, None, None, 49.52, 57.30, ''),
    ('MB41', 'Green syenite', '1', 208, 24, 13, None, None, None, None, None, None, None, 49.53, 57.29, ''),
    ('MB42', 'Fol green syenite', '1', 894, 75, 13, None, None, None, None, None, None, None, 50.45, 57.96, 'b'),
    ('MB43', 'Carbonatite', 'd', 201, 36, 9, 7, 0, 151.4, 33.4, 5.8, 108.0, 'B', 50.48, 57.96, ''),
    ('MB44', 'Green syenite', '1', 193, 23, 13, 7, 3, 124.0, 69.2, 7.3, 46.0, 'A', 49.69, 57.75, ''),
]
assert len(T1) == 44 and sum(r[5] for r in T1) == 515

# Table 1 rock type -> (geologic_types, lithologies, description stem)
ROCK = {
    'Fsp-por dyke': ('Volcanic Dike', 'Feldspar Porphyry', 'Feldspar porphyry dyke'),
    'Red fsp dyke': ('Volcanic Dike', 'Feldspar Porphyry',
                     'Red, NNE-trending, co-magmatic feldspar porphyry dyke'),
    'Syenite xeno': ('Pluton', 'Syenite', 'Syenite xenolith'),
    'Fol grey syenite': ('Pluton', 'Syenite', 'Foliated grey syenite'),
    'Grey syenite': ('Pluton', 'Syenite', 'Grey syenite'),
    'Pink fsp syenite': ('Pluton', 'Syenite', 'Pink feldspar syenite'),
    'Pink fs syenite': ('Pluton', 'Syenite', 'Pink feldspar syenite'),
    'Green syenite': ('Pluton', 'Syenite', 'Green syenite'),
    'Fol green syenite': ('Pluton', 'Syenite', 'Foliated green syenite'),
    'Red syenite': ('Pluton', 'Syenite', 'Red syenite'),
    'Fol red syenite': ('Pluton', 'Syenite', 'Foliated red syenite'),
    'Xeno red syenite': ('Pluton', 'Syenite', 'Xenolith-bearing red syenite'),
    'Xeno in syenite': ('Pluton', 'Syenite', 'Xenolith in syenite'),
    'Syenite dyke': ('Volcanic Dike', 'Syenite', 'Syenite dyke'),
    'Aplite dyke': ('Volcanic Dike', 'Aplite', 'Aplite dyke'),
    'Aplite + host': ('Volcanic Dike:Pluton', 'Aplite:Syenite', 'Aplite dyke and its syenite host'),
    'Mafic dyke': ('Volcanic Dike', 'Mafic Dike', 'Mafic dyke'),
    'Finer granite': ('Pluton', 'Granite', 'Finer-grained granite'),
    'Carbonatite': ('Volcanic Dike', 'Carbonatite', 'Ribbon carbonatite dyke'),
}
UNIT = {'1': 'intrusive unit 1 of the Baie des Moutons complex (Table 1)',
        '2': 'intrusive unit 2 of the Baie des Moutons complex (Table 1)',
        '3': 'intrusive unit 3 of the Baie des Moutons complex (Table 1)',
        'd': ('one of the late dykes that are the last magmatic event of the Baie des '
              'Moutons complex (Table 1)')}
FOOT = {'a': 'No site mean: two or more ChRM-bearing specimens have inconsistent directions.',
        'b': 'No site mean: too few ChRM-bearing specimens to calculate a site mean direction.',
        'c': 'Site of an attempted but inconclusive baked-contact test (no stable, '
             'non-viscous remanence in the host syenite).'}
NO_MEAN = ('No site mean: no detectable or stable ancient remanence; most such sites are '
           'dominated by a viscous remanence of recent origin carried by multidomain magnetite.')
AR_DATED = {'MB44': 'Hornblende, biotite and K-feldspar 40Ar/39Ar dated (Table 2).',
            'MB20': 'Hornblende, biotite and K-feldspar 40Ar/39Ar dated (Table 2).',
            'MB36': 'Hornblende and K-feldspar 40Ar/39Ar dated (Table 2).'}
SITE_NOTE = {
    'MB06': ('One of the four sites (MB06, MB11, MB12, MB36) of the four-dyke ChRM B mean, '
             'which the text describes as NNE-trending red feldspar porphyry dykes.'),
    'MB15': 'Site mean from great circles only; Table 1 gives k but no alpha95.',
    'MB20': ('Carries ChRM B in the polarity opposite to the other ChRM B sites (northwesterly '
             'and down); inverted by the authors for the six-site ChRM B mean.'),
    'MB27': ('Direction of the single ChRM-bearing endpoint specimen, not a site mean; not '
             'assigned to ChRM A or B.'),
    'MB43': ('Included by the authors in the ChRM B mean because its direction is similar to '
             'that of the dykes and it is late in the magmatic sequence.'),
}


def latlon(lat_min, lon_min):
    """Table 1 gives NAD-27 decimal minutes added to 50 N and 58 W."""
    return round(50 + lat_min / 60, 4), round(360 - (58 + lon_min / 60), 4)


rows = []
for (site, rock, unit, chi, nrm, n_tot, e, g, dec, inc, a95, k, typ, lat_m, lon_m,
     fn) in T1:
    gtype, lith, stem = ROCK[rock]
    lat, lon = latlon(lat_m, lon_m)
    desc = [f'{stem}, {UNIT[unit]}. Mean susceptibility {chi} x 1e-6 SI, mean NRM {nrm} '
            f'mA/m.']
    accepted = typ in ('A', 'B', '-B')
    r = {'site': site, 'location': LOCATION, 'result_type': 'i',
         'result_quality': 'g' if accepted else 'b',
         'citations': CIT, 'geologic_classes': 'Intrusive', 'geologic_types': gtype,
         'lithologies': lith, 'lat': lat, 'lon': lon, 'geographic_precision': 0.0002,
         **AGE, 'dir_n_total_specimens': n_tot}
    if dec is not None:
        codes = ['LP-DIR-AF', 'LP-DIR-T', 'DE-BFL'] if e else ['LP-DIR-AF', 'LP-DIR-T']
        if g:
            codes += ['DE-BFP', 'DE-FM-LP']
        elif e + g > 1:
            codes += ['DE-FM']
        r.update({'method_codes': ':'.join(codes), 'dir_tilt_correction': 0,
                  'dir_dec': dec, 'dir_inc': inc, 'dir_alpha95': a95, 'dir_k': k,
                  'dir_n_specimens': e + g, 'dir_n_specimens_lines': e,
                  'dir_n_specimens_planes': g})
        if accepted:
            plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95 if a95 else 0.0, lat, lon)
            r.update({'dir_comp_name': typ.lstrip('-'),
                      'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
                      'vgp_dp': round(dp, 1) if a95 else None,
                      'vgp_dm': round(dm, 1) if a95 else None})
            desc.append(f"Site mean of ChRM {typ.lstrip('-')}.")
    else:
        r['method_codes'] = 'LP-DIR-AF:LP-DIR-T'
        if not set(fn) & {'a', 'b'}:
            desc.append(NO_MEAN)
    if site in SITE_NOTE:
        desc.append(SITE_NOTE[site])
    for f in fn:
        desc.append(FOOT[f])
    if site in AR_DATED:
        desc.append(AR_DATED[site])
    r['description'] = ' '.join(desc)
    rows.append(r)

COLS = ['site', 'location', 'result_type', 'result_quality', 'method_codes', 'citations',
        'geologic_classes', 'geologic_types', 'lithologies', 'lat', 'lon',
        'geographic_precision', 'age', 'age_low', 'age_high', 'age_unit', 'dir_comp_name',
        'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k',
        'dir_n_specimens', 'dir_n_specimens_lines', 'dir_n_specimens_planes',
        'dir_n_total_specimens', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description']
sites = pd.DataFrame(rows)[COLS]
INT_COLS = ['dir_n_specimens', 'dir_n_specimens_lines', 'dir_n_specimens_planes',
            'dir_n_total_specimens', 'dir_tilt_correction', 'dir_n_sites', 'pole_n_sites',
            'pole_vv_q']

good = sites[sites['result_quality'] == 'g']
A = good[good['dir_comp_name'] == 'A']
B = good[good['dir_comp_name'] == 'B']
assert len(A) == 8 and len(B) == 6 and (sites['result_quality'] == 'b').sum() == 30


def minority_mask(sub):
    """True for sites whose direction points away from the majority polarity."""
    d = ipmag.fisher_mean(di_block=pmag.flip(
        ipmag.make_di_block(sub['dir_dec'].tolist(), sub['dir_inc'].tolist()), combine=True))
    ref = pmag.dir2cart([d['dec'], d['inc'], 1.])
    away = np.array([np.dot(pmag.dir2cart([r.dir_dec, r.dir_inc, 1.]), ref) < 0
                     for r in sub.itertuples()])
    return away if away.sum() <= len(away) / 2 else ~away


def majority_dir(sub):
    m = minority_mask(sub)
    decs = np.where(m, (sub['dir_dec'] + 180) % 360, sub['dir_dec'])
    incs = np.where(m, -sub['dir_inc'], sub['dir_inc'])
    return ipmag.fisher_mean(di_block=ipmag.make_di_block(list(decs), list(incs))), m


def unified_pole(sub, south):
    blk = pmag.flip(ipmag.make_di_block(sub['vgp_lon'].tolist(), sub['vgp_lat'].tolist()),
                    combine=True)
    p = ipmag.fisher_mean(di_block=blk)
    if (p['inc'] < 0) != south:
        p = ipmag.fisher_mean(di_block=ipmag.do_flip(di_block=blk))
    return p


def box(sub):
    return dict(lat_s=sub['lat'].min(), lat_n=sub['lat'].max(),
                lon_w=sub['lon'].min(), lon_e=sub['lon'].max())


def liths(sub):
    return ':'.join(dict.fromkeys(l for s in sub['lithologies'] for l in s.split(':')))


UMBRELLA_DESC = (
    "The Baie des Moutons (Mutton Bay) syenite complex, a ca. 20 km diameter alkalic complex of "
    "nested cones emplaced in Grenville gneiss on the Quebec north shore of the Gulf of St. "
    "Lawrence during late Neoproterozoic Iapetan rifting (McCausland et al., 2011). Its three "
    "main intrusive units (Lalonde and Martin, 1983) range from dark green feldspar "
    "porphyry-hornblende syenite (unit 1) through grey hornblende syenite (unit 2) to "
    "trachytic quartz syenite (unit 3); all are cut by aplite and microsyenite dykes and, "
    "last, by co-magmatic red NNE-trending feldspar porphyry dykes and minor carbonatite "
    "dykes. All 44 sites of Table 1 (530 specimens, typically six drill-core samples per "
    "site). Pilot specimens were demagnetized by stepwise AF to 120 mT or thermally to 580 C; "
    "characteristic directions were isolated by AF cleaning to 15-20 mT or low-temperature "
    "demagnetization, followed by thermal demagnetization to at least 450 C, and fitted by "
    "least-squares lines and, where a viscous overprint overlapped, great circles. Site means "
    "required five or more specimens from at least three samples and alpha95 < 15 deg. "
    "Low-temperature demagnetization was found to create a spurious present-field-like "
    "remanence in some specimens, which were excluded at sites with a large viscous "
    "component. Eight sites carry the steep component ChRM A and six the shallow component "
    "ChRM B; the two are nowhere observed together. Site MB27 has a single ChRM-bearing "
    "specimen, and the other 29 sites gave no site mean, mostly because a recent viscous "
    "remanence carried by multidomain magnetite dominates. Directions are in geographic "
    "coordinates: the steep, variably striking dykes indicate that the complex has not been "
    "significantly tilted since emplacement. Baked-contact tests at sites MB13-14, MB15, "
    "MB16, MB22, MB26 and MB36 were inconclusive because the host syenite lacks a stable "
    "non-viscous remanence. Age of the rock: the weighted mean of the four hornblende "
    "40Ar/39Ar plateau ages from sites MB44 (unit 1) and MB20 (unit 3), which the authors "
    "take as the intrusive age of the complex; the age range is its uncertainty, quoted at 2 "
    "sigma in Fig. 1 of the paper. Overlapping biotite and K-feldspar ages indicate rapid "
    "cooling after emplacement. Site coordinates are NAD-27 (Table 1).")
POLE_COMMON = (" The pole is the Fisher mean of the site virtual geomagnetic poles in "
               "geographic coordinates. The authors' pole, calculated from the mean direction "
               "at 50 deg 50.4' N, 58 deg 57.6' W, is {}. The age is that of the complex "
               "(weighted mean hornblende 40Ar/39Ar age; 2-sigma range).")
LOCS = [
    dict(comp='A', sub=A, south=False, gpmdb=9364,
         rname='Baie des Moutons syenite complex ChRM A pole',
         desc=("Paleomagnetic pole from the steep, easterly, down-directed characteristic "
               "component ChRM A, carried by pseudo-single-domain magnetite with unblocking "
               "temperatures above 450 C and coercivities above 20 mT, in syenite, aplite dyke "
               "and mafic dyke sites of intrusive units 1 to 3. The component is of a single "
               "polarity and lies close to the present-field direction, but was isolated with "
               "demagnetization designed to screen out present-field overprints and does not "
               "resemble Phanerozoic directions for the locality; the authors interpret it as "
               "most simply a primary thermoremanence acquired before ChRM B."
               + POLE_COMMON.format("BMa, 42.6 N 332.7 E (dp/dm 11.7/12.4)"))),
    dict(comp='B', sub=B, south=True, gpmdb=9365,
         rname='Baie des Moutons syenite complex ChRM B pole',
         desc=("Paleomagnetic pole from the shallow, south-southeasterly or northwesterly "
               "characteristic component ChRM B, removed between 260 C and 450-520 C in the "
               "feldspar porphyry dykes, found in late co-magmatic dykes (three red feldspar "
               "porphyry dykes, an aplite dyke and a ribbon carbonatite dyke) and in the late "
               "unit-3 red syenite MB20, whose direction is of the opposite polarity and is "
               "inverted for the mean. The authors interpret ChRM B as acquired during the "
               "waning magmatic stages of the intrusion; it resembles the ca. 565 Ma Sept-Iles "
               "A direction and Early Carboniferous directions for Laurentia. The authors also "
               "give a pole from the four NNE-trending dyke sites alone (MB06, MB11, MB12, "
               "MB36; -35.5 N 312.3 E, dp/dm 8.5/17.0) and a five-site mean adding MB43."
               + POLE_COMMON.format("BMb, -34.2 N 321.5 E (dp/dm 10.9/21.8)"))),
]

loc_rows = [{'location': LOCATION, 'location_type': 'Region', **GEOG,
             'sites': ':'.join(sites['site']), **box(sites), **AGE,
             'geologic_classes': 'Intrusive', 'lithologies': liths(sites),
             'method_codes': 'LP-DIR-AF:LP-DIR-T:LP-LT', 'citations': CIT,
             'description': UMBRELLA_DESC}]
poles = {}
for L in LOCS:
    sub = L['sub']
    p = unified_pole(sub, L['south'])
    d, minority = majority_dir(sub)
    poles[L['comp']] = (p, d, sub)
    codes = ['LP-DIR-AF', 'LP-DIR-T', 'LP-LT', 'DE-BFL', 'DE-BFP', 'DE-FM-LP', 'DE-VGP']
    loc_rows.append({
        'location': LOCATION, 'location_type': 'Region', **GEOG,
        'result_name': L['rname'], 'result_type': 'a', 'result_quality': 'g',
        'sites': ':'.join(sub['site']), **box(sub), **AGE,
        'geologic_classes': 'Intrusive', 'lithologies': liths(sub),
        'method_codes': ':'.join(codes), 'citations': CIT,
        'dir_tilt_correction': 0,
        'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
        'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1),
        'dir_n_sites': int(d['n']),
        'pole_comp_name': L['comp'],
        'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
        'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
        'pole_n_sites': int(p['n']),
        'pole_reversed_perc': round(100 * minority.sum() / len(sub), 1),
        'pole_vv_q': 5, 'external_database_ids': f"GPMDB[{L['gpmdb']}]",
        'description': L['desc']})
locs = pd.DataFrame(loc_rows)

# ---- ages table: 40Ar/39Ar plateau ages (Table 2; 1 sigma) and the weighted mean ----------
TS = dict(timescale_eon='Proterozoic', timescale_era='Neoproterozoic',
          timescale_period='Ediacaran')
AGE_ROWS = [
    ('MB44', 584.1, 1.6, 'hornblende grain (Fig. 6a)', 'unit 1 green syenite'),
    ('MB44', 582.5, 1.6, 'replicate hornblende grain', 'unit 1 green syenite'),
    ('MB44', 591.6, 1.7, 'biotite grain C81a (Fig. 6b)', 'unit 1 green syenite'),
    ('MB44', 586.2, 1.7, 'biotite grain C81b', 'unit 1 green syenite'),
    ('MB20', 582.0, 1.8, 'hornblende grain B66a (Fig. 6d)', 'unit 3 red syenite'),
    ('MB20', 585.0, 1.8, 'hornblende grain B66b', 'unit 3 red syenite'),
    ('MB20', 583.1, 1.6, 'biotite grain B68a (Fig. 6e)', 'unit 3 red syenite'),
    ('MB20', 584.4, 1.6, 'biotite grain B68b', 'unit 3 red syenite'),
    ('MB36', 583.5, 4.2, 'hornblende grain C83c (Fig. 6h)', 'late red feldspar porphyry dyke'),
]
ages = [{'location': LOCATION, 'site': s, 'age': a, 'age_sigma': sig, 'age_unit': 'Ma',
         'method_codes': 'GM-ARAR-AP', 'citations': CIT, **TS,
         'description': (f'40Ar/39Ar single-grain laser step-heating plateau age, {what}, '
                         f'site {s} ({unit}); 1 sigma, excluding decay-constant and '
                         f'standard uncertainties (MMhb-1 monitor at 520.4 Ma).')}
        for s, a, sig, what, unit in AGE_ROWS]
ages.append({'location': LOCATION, 'site': None, 'age': 583.4, 'age_sigma': 1.0,
             'age_unit': 'Ma', 'method_codes': 'GM-ARAR-AP', 'citations': CIT, **TS,
             'description': ('Weighted mean of the four hornblende 40Ar/39Ar plateau ages '
                             'from sites MB44 (unit 1) and MB20 (unit 3), 583.4 +/- 2.0 Ma, '
                             'interpreted by the authors as the best estimate of the intrusive '
                             'age of the complex. Fig. 1 of the paper shows it as 583 +/- 2 Ma '
                             'with age errors quoted at 2 sigma; age_sigma is 1 sigma. Excludes '
                             'decay-constant and standard uncertainties.')})
ages = pd.DataFrame(ages)


def fmt(df):
    df = df.dropna(axis=1, how='all').copy()
    for c in df.columns:
        if c in INT_COLS:
            df[c] = df[c].map(lambda x: '' if pd.isna(x) else str(int(x)))
        elif df[c].dtype.kind == 'f':
            df[c] = df[c].map(lambda x: '' if pd.isna(x) else f'{x:.10g}')
    return df.fillna('')


def table_text(df, kind):
    return f'tab delimited\t{kind}\n' + fmt(df).to_csv(sep='\t', index=False)


TABLES = (('locations', locs), ('sites', sites), ('ages', ages))
for kind, df in TABLES:
    with open(os.path.join(OUT, f'{kind}.txt'), 'w', encoding='utf-8') as f:
        f.write(table_text(df, kind))
stamp = date.today().strftime('%d.%b.%Y')
upload = os.path.join(OUT, f'McCausland2011_Baie_des_Moutons_{stamp}.txt')
with open(upload, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {upload}')

for comp, (p, d, sub) in poles.items():
    print(f'  ChRM {comp}: B {int(p["n"])} N {int(sub["dir_n_specimens"].sum())} pole '
          f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} K {p["k"]:.1f} | dir '
          f'{d["dec"]:.1f}/{d["inc"]:.1f} a95 {d["alpha95"]:.1f} k {d["k"]:.1f}')

# poles of the Swanson-Hysell et al. (2026) compilation (lat, lon, A95, B, N specimens) and
# the published mean directions (dec, inc, alpha95, k)
EXPECT = {'A': ((42.9, 331.9, 11.5, 8, 55), (98.6, 78.0, 6.6, 71.8)),
          'B': ((-33.9, 321.9, 15.7, 6, 44), (163.1, 6.0, 21.7, 10.5))}
for comp, (p, d, sub) in poles.items():
    got = ((round(p['inc'], 1), round(p['dec'], 1), round(p['alpha95'], 1), int(p['n']),
            int(sub['dir_n_specimens'].sum())),
           (round(d['dec'], 1), round(d['inc'], 1), round(d['alpha95'], 1), round(d['k'], 1)))
    assert got == EXPECT[comp], (comp, got)
assert locs.loc[2, 'pole_reversed_perc'] == 16.7
print('-I- both poles match the compilation')

sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True)
