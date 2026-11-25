"""Build the MagIC contribution for the Adirondack Highlands paleomagnetic data of Brown and McEnroe (2012).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Brown, L. L., & McEnroe, S. A. (2012). Paleomagnetism and magnetic mineralogy of Grenville
metamorphic and igneous rocks, Adirondack Highlands, USA. Precambrian Research 212-213, 57-74.
doi:10.1016/j.precamres.2012.04.012

Age of the metamorphic anorthosite pole: Zhang, Y., Hodgin, E. B., Schmitz, M. D., Schwartz,
D. M., Mohr, M. T., Crowley, J. L., & Swanson-Hysell, N. L. (2026). Protracted high-grade
metamorphism and slow cooling in the Adirondacks recalibrate Neoproterozoic global
paleogeography. Tectonics 45, e2025TC009333. doi:10.1029/2025TC009333
Fayalite granite emplacement ages: McLelland, J., Hamilton, M., Selleck, B., McLelland, J.,
Walker, D., & Orrell, S. (2001). Precambrian Research 109, 39-72.
doi:10.1016/S0301-9268(01)00141-3 (uncertainties at 95% confidence; age_sigma is half).

Sources and checks
------------------
- Table 1 (p. 65; site means and site VGPs of 36 sites) was transcribed from the page image.
  The 13 ADH sites are from an earlier unpublished study of R. Hargraves (per. comm., 2002),
  reported in Table 1 with N only; the 23 AD sites were measured at the University of
  Massachusetts Amherst (section 4.1).
- Site VGPs are the authors' tabulated values. Those of the AD sites follow from the site
  directions at a common point near 44.0N, 76.0W, and those of the ADH sites from the site
  directions at the individual (estimated) site positions (asserted below, within 1.1 deg).
- Two ADH sites do not follow this pattern: Table 1 prints the same direction for ADH 12 and
  ADH 34 (D 240.7, I -73.3), and ADH 34's VGP follows from it, ADH 12's does not; the VGP of
  ADH 10 does not follow from its direction either. Both are stored as printed and described.
  As printed, the anorthosite (11 sites) and combined (14 sites) mean directions of Table 2
  are not reproduced from the Table 1 directions; the microcline gneiss, other metamorphic
  rocks and fayalite granite means are (asserted below).
- Site coordinates: the paper gives none. The 11 anorthosite sites carry the positions
  provided by L. Brown (personal communication) and plotted in Zhang et al. (2026, Fig. 2b).
  The other sites carry the position of their sampling area on Fig. 1 (letters R, O, S, W,
  georeferenced from the map's graticule ticks), the village of Au Sable Forks for ADH 20,
  the centroid of R, O and S for the gneiss sites the paper does not place in an area, and
  the middle of the sampled traverse for the three associated metamorphic rock sites;
  geographic_precision gives the uncertainty.

Structure and interpretation
----------------------------
- Three locations, one per unit. No tilt correction. Every site is in a pole.
- Adirondack metamorphic anorthosite: the 11 metamorphosed anorthosite sites and the three
  associated metamorphic rock sites ('Combined meta. rocks' of Table 2; 'metamorphic
  anorthosites and other rocks' of Table 3), with two pole rows: the Fisher mean of the
  tabulated site VGPs (polarity unified, reported in the reversed, majority polarity), which
  is the pole used in the Swanson-Hysell et al. (2026) compilation, and the authors' pole
  calculated from the mean direction (Table 2; GPMDB 9394). The positive inverse
  baked-contact test of Zhang et al. (2026) at a mafic dike (site GD1) cutting the Marcy
  anorthosite (site MA11) shows the host remanence predates the dike (contact_test IC+).
- Adirondack microcline gneiss and Adirondack fayalite granite: the authors' poles
  calculated from the mean directions (Table 2; GPMDB 9393 and 9395).
- Ages. The remanence of these rocks was acquired on cooling after granulite-facies
  metamorphism, so the age given for each site and pole is the time of cooling through the
  remanence acquisition temperatures, not the protolith or emplacement age:
  - magnetite in the anorthosite and associated rocks: 887 Ma (864-910 Ma), the remanence
    acquisition age of Zhang et al. (2026) from a thermal history inversion of U-Pb zircon,
    garnet, titanite, rutile and apatite data from the Marcy massif (cooling through about
    570-525 C);
  - ilmeno-hematite in the microcline gneiss: 850-960 Ma, a broad bracket because the
    timing of lamellar exsolution is uncertain. 960 Ma is the authors' maximum acquisition
    age (lamellae form below about 520 C; 500 C taken as the mean acquisition temperature on
    the published cooling curves for the western Highlands); exsolution continued to about
    420 C (sections 6.1 and 7). 850 Ma is the youngest of the rutile and apatite U-Pb
    cooling dates of Zhang et al. (2026) from the Marcy massif;
  - magnetite in the fayalite granites: 888-1047 Ma, from U-Pb zircon emplacement
    (McLelland et al., 2001) to the young 95% limit of the time at which the Marcy massif
    cooled through 570 C in the thermal history of Zhang et al. (2026; 896 +14/-8 Ma); the
    authors argue that the granites reached 570 C before the hotter metamorphic rocks of the
    central Highlands.
  The authors' cooling-curve estimates (granite ca. 990 Ma, anorthosite ca. 970 Ma, gneiss
  ca. 960 Ma) are given in the descriptions.
- The ages table holds measured dates: the U-Pb zircon emplacement ages of the two fayalite
  granite bodies and the ID-TIMS dates of Zhang et al. (2026) on the Marcy anorthosite (the
  three concordant igneous zircon dates, the youngest metamorphic zircon and the ranges of
  the rutile and apatite cooling dates).
"""
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(OUTDIR))

CIT = '10.1016/j.precamres.2012.04.012'
ZHANG2026 = '10.1029/2025TC009333'
MCLELLAND2001 = '10.1016/S0301-9268(01)00141-3'
MEZGER1991 = '10.1086/629503'
GEOG = dict(continent_ocean='North America', country='United States of America',
            state_province='New York', terranes='Laurentia',
            geological_province_sections='Grenville Province')
AD_MC = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM'  # section 4.1
ADH_MC = 'DE-FM'  # Hargraves sites: laboratory procedures not given

# Table 1 (p. 65): site, unit, N, No (AD sites), susceptibility (SI), NRM (A/m), INC, DEC, k,
# alpha95, VGP lat, VGP lon. Units: GMS microcline gneiss, ANT metamorphosed anorthosites,
# OMR associated metamorphic rocks, FG post-metamorphic fayalite granites.
TABLE1 = [
    ('AD1', 'GMS', 6, 7, '0.00059', '5.16', -62.3, 232.1, 100, 6.7, -53.1, 175.8),
    ('AD2', 'GMS', 4, 6, '0.00090', '0.86', -62.9, 257.8, 97, 9.4, -36.4, 164.3),
    ('AD5', 'GMS', 10, 10, '0.00707', '1.49', -51.5, 292.2, 22, 10.6, -8.0, 157.4),
    ('AD6', 'GMS', 10, 10, '0.04192', '4.92', -67.8, 271.2, 260, 3.0, -31.9, 152.1),
    ('AD11', 'GMS', 7, 12, '0.02508', '2.21', 44.8, 98.1, 120, 5.5, 12.6, 349.3),
    ('AD15', 'GMS', 3, 3, '0.00457', '0.15', 68.9, 138.7, 351, 6.6, 12.7, 308.5),
    ('AD18', 'GMS', 5, 5, '0.00224', '0.09', 55.3, 111.6, 165, 6.0, 11.0, 334.2),
    ('AD21', 'GMS', 7, 7, '0.12257', '7.02', -64.1, 318.4, 655, 2.4, -7.1, 131.8),
    ('AD22', 'GMS', 5, 5, '0.08142', '6.39', -63.3, 282.3, 113, 7.2, -22.4, 152.6),
    ('AD26', 'GMS', 5, 5, '0.01485', '1.69', 76.9, 162.9, 760, 2.8, 19.9, 291.6),
    ('AD34', 'GMS', 7, 8, '0.00163', '7.61', -61.1, 262.7, 62, 7.7, -32.3, 164.4),
    ('AD35', 'GMS', 3, 3, '0.02948', '1.32', -58.6, 309.9, 238, 8.0, -4.8, 140.6),
    ('AD36', 'GMS', 5, 5, '0.00011', '0.15', -51.2, 316.9, 55, 10.4, 4.6, 139.6),
    ('AD40', 'GMS', 3, 3, '0.00103', '0.02', -55.6, 294.2, 130, 10.9, -9.8, 152.4),
    ('AD12', 'ANT', 5, 5, '0.08083', '1.88', -61.5, 304.4, 44, 11.7, -10.0, 142.0),
    ('AD28', 'ANT', 5, 5, '0.00462', '0.20', 51.2, 101.6, 67, 9.4, 14.1, 343.1),
    ('AD32', 'ANT', 3, 3, '0.00158', '1.24', -87.5, 182.2, 926, 4.1, -49.1, 104.3),
    ('AD33', 'ANT', 5, 5, '0.00055', '0.14', -65.6, 290.4, 153, 6.2, -20.3, 146.2),
    ('ADH9', 'ANT', 5, None, None, None, -75.6, 312.8, 210, 4.3, -23.5, 128.0),
    ('ADH10', 'ANT', 4, None, None, None, -76.0, 279.5, 177, 5.3, -19.5, 105.6),
    ('ADH12', 'ANT', 5, None, None, None, -73.3, 240.7, 798, 2.2, -41.6, 130.1),
    ('ADH32', 'ANT', 6, None, None, None, 43.9, 100.3, 107, 5.5, 10.7, 350.7),
    ('ADH33', 'ANT', 6, None, None, None, 46.6, 82.5, 16, 14.3, 24.0, 359.7),
    ('ADH34', 'ANT', 6, None, None, None, -73.3, 240.7, 49, 8.2, -51.0, 151.4),
    ('ADH35', 'ANT', 6, None, None, None, 53.9, 108.0, 25, 11.5, 11.9, 339.1),
    ('ADH8', 'OMR', 4, None, None, None, -75.6, 283.1, 117, 5.3, -33.4, 138.0),
    ('ADH18', 'OMR', 4, None, None, None, 72.8, 146.3, 33, 12.2, 16.5, 301.9),
    ('ADH19', 'OMR', 4, None, None, None, -68.9, 264.9, 72, 8.2, -36.5, 154.7),
    ('AD9', 'FG', 7, 7, '0.00616', '0.69', -75.6, 271.8, 289, 3.6, -37.4, 139.1),
    ('AD10', 'FG', 4, 4, '0.01989', '1.04', -78.3, 269.7, 2409, 1.9, -40.0, 133.9),
    ('AD20', 'FG', 4, 4, '0.00553', '0.69', -76.4, 317.3, 2157, 2.0, -23.3, 122.7),
    ('AD37', 'FG', 5, 5, '0.00552', '0.49', -74.5, 291.0, 99, 7.7, -28.9, 135.1),
    ('AD38', 'FG', 7, 7, '0.00201', '0.51', -69.9, 320.7, 244, 3.9, -13.4, 126.6),
    ('ADH4', 'FG', 4, None, None, None, -75.7, 286.3, 847, 2.4, -31.9, 136.0),
    ('ADH5', 'FG', 4, None, None, None, -78.7, 309.8, 1406, 1.9, -28.4, 124.1),
    ('ADH20', 'FG', 5, None, None, None, 71.7, 118.7, 177, 13.7, 23.2, 318.1),
]
# Table 2 (p. 67) means reproduced from the Table 1 directions: N, DEC, INC, alpha95
TABLE2 = {'GMS': (14, 289.2, -62.8, 7.6), 'OMR': (3, 289.5, -74.1, 15.2),
          'FG': (8, 297.0, -75.8, 3.9)}

# Site positions (lat N, lon E, geographic_precision deg, source)
PERS = 'provided by L. Brown (personal communication) and plotted in Zhang et al. (2026, Fig. 2b)'
FIG1 = 'on Fig. 1, georeferenced from the map graticule'
AREA = {
    'R': (44.31, 284.91, 0.15, f'the Russell area (letter R {FIG1}); section 3.3.1 places AD 1-6 '
                               f'and 34 in the Russell area'),
    'O': (44.17, 284.95, 0.1, f'the Oswegatchie area (letter O {FIG1}); the captions of Figs. 7c-d '
                              f'and 11e place AD6 in the Oswegatchie area'),
    'S': (44.13, 285.32, 0.1, f'the Sabattis area (letter S {FIG1}); section 3.3.1 places AD 15 and '
                              f'18 in the Sabattis area'),
    'GMS': (44.21, 285.06, 0.3, f'the centroid of the gneiss sampling areas (letters R, O and S {FIG1}); '
                                f'the paper does not say in which area this site lies'),
    'W': (44.17, 285.09, 0.1, f'the Wanakena granite (letter W {FIG1}); all Wanakena samples were '
                              f'collected west of the town of Wanakena along State Highway 3'),
    'A': (44.44, 286.32, 0.1, 'the village of Au Sable Forks, near which the paper places the second '
                              'fayalite granite body'),
    'HL': (44.25, 285.65, 0.8, 'the middle of the sampled east-west Highlands traverse; the paper '
                               'gives no location for the associated metamorphic rocks'),
}
SITE_POS = {
    'AD12': (44.23, 286.12), 'AD28': (44.25, 285.75), 'AD32': (44.29, 286.47),
    'AD33': (44.30, 286.60), 'ADH9': (44.30, 286.52), 'ADH10': (44.27, 286.38),
    'ADH12': (44.00, 286.00), 'ADH32': (44.12, 286.25), 'ADH33': (43.61, 286.13),
    'ADH34': (44.02, 285.95), 'ADH35': (43.67, 285.95)}
SITE_AREA = {'AD1': 'R', 'AD2': 'R', 'AD5': 'R', 'AD34': 'R', 'AD6': 'O', 'AD15': 'S',
             'AD18': 'S', 'AD11': 'GMS', 'AD21': 'GMS', 'AD22': 'GMS', 'AD26': 'GMS',
             'AD35': 'GMS', 'AD36': 'GMS', 'AD40': 'GMS', 'AD9': 'W', 'AD10': 'W', 'AD20': 'W',
             'AD37': 'W', 'AD38': 'W', 'ADH4': 'W', 'ADH5': 'W', 'ADH20': 'A', 'ADH8': 'HL',
             'ADH18': 'HL', 'ADH19': 'HL'}
AD_VGP_POINT = (44.0, 284.0)  # common point at which the AD-site VGPs follow from the directions

LOC = {'GMS': 'Adirondack microcline gneiss', 'ANT': 'Adirondack metamorphic anorthosite',
       'OMR': 'Adirondack metamorphic anorthosite', 'FG': 'Adirondack fayalite granite'}
# Ages of cooling through the remanence acquisition temperatures (see the docstring)
AGE = {'GMS': dict(age=None, age_low=850, age_high=960),
       'ANT': dict(age=887, age_low=864, age_high=910),
       'OMR': dict(age=887, age_low=864, age_high=910),
       'FG': dict(age=None, age_low=888, age_high=1047)}
AGE_DESC = {
    'ANT': ('Age: remanence acquisition on cooling through about 570-525 C, 887 Ma (864-910 Ma), '
            'from the thermal history of the Marcy massif of Zhang et al. (2026).'),
    'GMS': ('Age: cooling through the ilmeno-hematite exsolution temperatures (below about '
            '520 C, continuing to about 420 C), 850-960 Ma; not a measured age. 960 Ma is the '
            'maximum acquisition age of Brown and McEnroe (2012) from published cooling '
            'curves; 850 Ma is the youngest rutile and apatite U-Pb cooling date of Zhang et '
            'al. (2026) from the Marcy massif.'),
    'FG': ('Age: cooling of magnetite through about 570 C, 888-1047 Ma; not a measured age. '
           'From U-Pb zircon emplacement (McLelland et al., 2001) to the young limit of the '
           'time the Marcy massif cooled through 570 C (Zhang et al., 2026); Brown and '
           'McEnroe (2012) estimated ca. 990 Ma from published cooling curves.'),
}
AGE_DESC['OMR'] = AGE_DESC['ANT']
ROCK = {
    'GMS': ('Metamorphic', 'Not Specified', 'Gneiss',
            'Microcline gneiss (GMS unit; sillimanite-microcline gneiss), western Adirondack '
            'Highlands; remanence carried by ilmeno-hematite.'),
    'ANT': ('Metamorphic', 'Pluton', 'Meta Anorthosite',
            'Metamorphosed anorthosite, Marcy massif or an outlying body; remanence carried by '
            'magnetite.'),
    'OMR': ('Metamorphic', 'Not Specified', 'Not Specified',
            'Associated metamorphic rock (Table 1; the other units sampled were metamorphosed '
            'syenites, pyroxene diorites and granitic gneisses, and Table 1 does not say which).'),
    'FG': ('Intrusive', 'Pluton', 'Granite',
           'Post-metamorphic fayalite-ferrohedenbergite granite; remanence carried by magnetite.'),
}
SPECIAL = {
    'ADH10': ('The tabulated VGP does not follow from the tabulated direction at the site '
              'position (they differ by 33 deg); both are stored as printed.'),
    'ADH12': ('Table 1 prints the same direction for ADH 12 and ADH 34; the tabulated VGP of '
              'ADH 12 does not follow from it at the site position (17 deg apart). Direction '
              'and VGP are stored as printed.'),
    'ADH34': 'Table 1 prints the same direction for ADH 12 and ADH 34.',
}

rows = []
for site, unit, N, No, susc, nrm, inc, dec, k, a95, vlat, vlon in TABLE1:
    cls, typ, lith, what = ROCK[unit]
    desc = [what]
    if site.startswith('ADH'):
        desc.append('Site of R. Hargraves (per. comm., 2002), reported in Table 1.')
    else:
        desc.append(f'Site mean susceptibility {susc} SI, NRM {nrm} A/m (Table 1).')
    if site in SITE_POS:
        (lat, lon), prec = SITE_POS[site], None
        desc.append(f'Coordinates {PERS}.')
    else:
        lat, lon, prec, where = AREA[SITE_AREA[site]]
        desc.append(f'The paper gives no site coordinates; the site position is that of {where}.')
    if site in SPECIAL:
        desc.append(SPECIAL[site])
    elif site.startswith('ADH'):
        desc.append('VGP as tabulated by the authors.')
    else:
        desc.append('VGP as tabulated by the authors; like those of the other AD sites, it '
                    'follows from the site direction at 44.0N, 76.0W.')
    desc.append(AGE_DESC[unit])
    r = {'site': site, 'location': LOC[unit], 'result_type': 'i', 'result_quality': 'g',
         'method_codes': ADH_MC if site.startswith('ADH') else AD_MC, 'citations': CIT,
         'geologic_classes': cls, 'geologic_types': typ, 'lithologies': lith,
         'lat': lat, 'lon': lon, 'geographic_precision': prec, **AGE[unit], 'age_unit': 'Ma',
         'dir_tilt_correction': 0, 'dir_dec': dec, 'dir_inc': inc, 'dir_k': float(k),
         'dir_alpha95': a95, 'dir_n_samples': N, 'dir_n_total_samples': No,
         'dir_polarity': 'n' if inc > 0 else 'r', 'vgp_lat': vlat, 'vgp_lon': vlon,
         'description': ' '.join(desc), '_unit': unit}
    rows.append(r)
sites = pd.DataFrame(rows)
sites['dir_n_total_samples'] = [None if pd.isna(x) else str(int(x))
                                 for x in sites['dir_n_total_samples']]


def misfit(r, lat, lon):
    plon, plat, _, _ = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], 1, lat, lon)
    return pmag.angle([r['vgp_lon'], r['vgp_lat']], [plon, plat])[0]


# tabulated VGPs follow from the directions (AD: common point; ADH: site positions)
for _, r in sites.iterrows():
    if r['site'] in SPECIAL and r['site'] != 'ADH34':
        assert misfit(r, r['lat'], r['lon']) > 15, r['site']
    elif r['site'].startswith('ADH') and r['site'] in SITE_POS:
        assert misfit(r, r['lat'], r['lon']) < 0.2, r['site']
    elif r['site'].startswith('ADH'):
        assert misfit(r, r['lat'], r['lon']) < 2.0, r['site']
    else:
        assert misfit(r, *AD_VGP_POINT) < 1.1, r['site']


def unified(sub, cols, reverse=True):
    """Site vectors (lon/dec, lat/inc) with normal-polarity sites taken to the reversed mode."""
    out = []
    for _, r in sub.iterrows():
        a, b = r[cols[0]], r[cols[1]]
        if r['dir_polarity'] == 'n' and reverse:
            a, b = (a + 180) % 360, -b
        out.append([a, b])
    return out


# Table 2 means reproduced (microcline gneiss, other metamorphic rocks, fayalite granite)
for unit, (n, D, I, a95) in TABLE2.items():
    d = pmag.fisher_mean(unified(sites[sites['_unit'] == unit], ['dir_dec', 'dir_inc']))
    assert int(d['n']) == n and abs(d['dec'] - D) < 0.06 and abs(d['inc'] - I) < 0.06, unit
    assert abs(d['alpha95'] - a95) < 0.15, unit

# Table 2 (p. 67) means and the authors' poles calculated from them: N, DEC, INC, k, alpha95,
# pole lat, pole lon, dp, dm ('Combined meta. rocks' for the anorthosite location)
TABLE2_POLES = {'ANT': (14, 283.9, -67.3, 28, 7.7, -25.1, 149.0, 10.6, 12.7),
                'GMS': (14, 289.2, -62.8, 29, 7.6, -18.4, 151.1, 9.3, 11.8),
                'FG': (8, 297.0, -75.8, 199, 3.9, -28.4, 132.7, 6.7, 7.2)}
# the published dp/dm follow from the mean inclinations and alpha95 (within rounding)
for unit, (n, D, I, k, a95, plat, plon, dp, dm) in TABLE2_POLES.items():
    lam = np.radians(pmag.plat(I))  # paleolatitude
    cdp = a95 * (1 + 3 * np.sin(lam) ** 2) / 2
    cdm = a95 * np.cos(lam) / np.cos(np.radians(I))
    assert abs(cdp - dp) < 0.15 and abs(cdm - dm) < 0.15, (unit, cdp, cdm)

AGE_REMARK = {
    'ANT': ('Zhang et al. (2026) report a positive inverse baked-contact test at a mafic dike '
            'cutting the Marcy anorthosite, showing that the anorthosite magnetite remanence '
            'predates the dike. From U-Pb petrochronology, thermochronology and a thermal '
            'history inversion of the Marcy massif, which crystallized ca. 1159 Ma, they date '
            'the remanence acquisition, on cooling through about 570-525 C, to 910-864 Ma '
            '(887 +/- 23 Ma); Brown and McEnroe (2012) had estimated ca. 970 Ma from published '
            'cooling curves.'),
    'GMS': ('The remanence was acquired as ilmenite lamellae exsolved from titanohematite '
            'below about 520 C, continuing to about 420 C. Brown and McEnroe (2012) estimated '
            'ca. 960 Ma (maximum, for 500 C) and ca. 940 Ma (for 420 C) from published '
            'cooling curves for the western Highlands. The age bracket given, 850-960 Ma, '
            'extends to the youngest rutile and apatite U-Pb cooling dates of Zhang et al. '
            '(2026), whose thermal history of the Marcy massif places magnetite remanence '
            'acquisition there ca. 85 Myr later than the cooling-curve estimate.'),
    'FG': ('The age bracket, 888-1047 Ma, runs from emplacement to the young limit of the time '
           'at which the Marcy massif cooled through 570 C in the thermal history of Zhang et '
           'al. (2026); Brown and McEnroe (2012) argue that the granites reached the magnetite '
           'blocking temperature before the hotter metamorphic rocks of the central Highlands, '
           'and estimated ca. 990 Ma from published cooling curves.'),
}
UNIT_TEXT = {
    'ANT': ('the metamorphosed anorthosites of the Marcy massif and outlying bodies together '
            'with three sites of associated metamorphic rocks (the combined metamorphic rocks '
            'of Brown and McEnroe, 2012), Adirondack Highlands, New York. Ten of the 14 sites '
            'are from an earlier unpublished study of R. Hargraves (per. comm., 2002). '
            'Granulite-facies metamorphism reset earlier remanence; after removal of a '
            'present-field overprint by AF or thermal demagnetization (400-525 C), the '
            'characteristic remanence is carried by near end-member multidomain and '
            'pseudo-single-domain magnetite, with subsidiary hemo-ilmenite, and was acquired on '
            'cooling. Both polarities are present, and the normal and reversed site means do '
            'not pass the reversal test at the 95% level (Brown and McEnroe, 2012).'),
    'GMS': ('the sillimanite-microcline gneisses (GMS unit) of the Russell, Oswegatchie and '
            'Sabattis areas, western Adirondack Highlands, New York. The characteristic '
            'remanence is very stable to AF demagnetization and unblocks between 600 and 680 C; '
            'it is a chemical remanence (lamellar magnetism) carried by ilmeno-hematite, with '
            'multidomain magnetite also present at Oswegatchie. Both polarities are present, '
            'and the normal and reversed site means narrowly fail the reversal test at the 95% '
            'level (Brown and McEnroe, 2012).'),
    'FG': ('the post-metamorphic fayalite-ferrohedenbergite granites of the Adirondack '
           'Highlands, New York: the Wanakena granite (seven reversed sites) and the fayalite '
           'granite near Au Sable Forks (one normal site, ADH 20), both undeformed and dated at '
           'ca. 1047 Ma by U-Pb zircon (McLelland et al., 2001). The remanence is carried by '
           'magnetite, with hemo-ilmenite, and unblocks mostly below 575 C. The single normal '
           'site and the reversed sites pass the reversal test at class C (Brown and McEnroe, '
           '2012).'),
}
RESULT = {'ANT': 'Adirondack metamorphic anorthosite and associated rocks pole',
          'GMS': 'Adirondack microcline gneiss pole', 'FG': 'Adirondack fayalite granite pole'}
CITS = {'ANT': f'{CIT}:{ZHANG2026}', 'GMS': f'{CIT}:{ZHANG2026}',
        'FG': f'{CIT}:{MCLELLAND2001}:{ZHANG2026}'}
GPMDB = {'ANT': 9394, 'GMS': 9393, 'FG': 9395}
EXTRA = {'ANT': dict(contact_test='IC+'), 'GMS': {}, 'FG': dict(reversal_test='Rc')}
TEST_MC = {'ANT': ':ST-C-I', 'GMS': '', 'FG': ':ST-R'}
LAB_MC = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM'

loc_rows, poles = [], {}
for unit in ('ANT', 'GMS', 'FG'):
    sub = sites[sites['location'] == LOC[unit]]
    vgps = unified(sub, ['vgp_lon', 'vgp_lat'])
    p = pmag.fisher_mean(vgps)
    # same pole as pmag.flip(combine=True) followed by a flip to the southern pole
    q = pmag.fisher_mean(pmag.flip(sub[['vgp_lon', 'vgp_lat']].values.tolist(), combine=True))
    assert pmag.angle([p['dec'], p['inc']], [(q['dec'] + 180) % 360, -q['inc']])[0] < 1e-6
    d = pmag.fisher_mean(unified(sub, ['dir_dec', 'dir_inc']))
    poles[unit] = (p, d, sub)
    rev = round(100 * (sub['dir_polarity'] == 'r').mean(), 1)
    common = {'location': LOC[unit], 'location_type': 'Region', 'result_type': 'a',
              'result_quality': 'g', 'sites': ':'.join(sub['site']), 'citations': CITS[unit],
              'geologic_classes': ':'.join(dict.fromkeys(sub['geologic_classes'])),
              'lithologies': ':'.join(dict.fromkeys(sub['lithologies'])),
              'lat_s': sub['lat'].min(), 'lat_n': sub['lat'].max(),
              'lon_w': sub['lon'].min(), 'lon_e': sub['lon'].max(), **AGE[unit],
              'age_unit': 'Ma', 'dir_tilt_correction': 0, 'pole_reversed_perc': rev,
              **EXTRA[unit], **GEOG}
    n, D, I, k, a95, plat, plon, dp, dm = TABLE2_POLES[unit]
    vgp_mean = (f'The Fisher mean of the tabulated site VGPs is {abs(p["inc"]):.1f}S '
                f'{p["dec"]:.1f}E (A95 {p["alpha95"]:.1f}).')
    if unit == 'ANT':
        # the compilation pole: Fisher mean of the tabulated site VGPs
        loc_rows.append({
            **common, 'result_name': RESULT[unit] + ' (site VGPs)',
            'method_codes': LAB_MC + ':DE-VGP' + TEST_MC[unit],
            'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
            'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1),
            'dir_n_sites': int(d['n']),
            'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
            'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
            'pole_n_sites': int(p['n']),
            'description': (f'Paleomagnetic pole for {UNIT_TEXT[unit]} {AGE_REMARK[unit]} The '
                            'pole is the Fisher mean of the tabulated site VGPs with polarity '
                            'unified, without tilt correction; the directions and VGPs of '
                            'sites ADH 10 and ADH 12 are inconsistent as tabulated (see the '
                            'site descriptions). The mean direction is that of the tabulated '
                            'site directions.'),
            '_unit': unit})
        authors = ('Pole of Brown and McEnroe (2012, Table 2) for the combined metamorphic rocks '
                   'of this location (11 metamorphosed anorthosite and 3 associated metamorphic '
                   'rock sites), calculated from their mean direction; mean direction, pole and '
                   'dp/dm as published. The age is the remanence acquisition age of Zhang et al. '
                   '(2026), as for the pole from the site VGPs in this location.')
    else:
        authors = (f'Paleomagnetic pole of Brown and McEnroe (2012, Table 2) for '
                   f'{UNIT_TEXT[unit]} {AGE_REMARK[unit]} The pole is calculated from the mean '
                   f'direction, without tilt correction; mean direction, pole and dp/dm as '
                   f'published. {vgp_mean}')
    loc_rows.append({
        **common, 'result_name': RESULT[unit] + (' (Brown and McEnroe, 2012)' if unit == 'ANT'
                                                 else ''),
        'method_codes': LAB_MC + TEST_MC[unit],
        'dir_dec': D, 'dir_inc': I, 'dir_alpha95': a95, 'dir_k': float(k), 'dir_n_sites': n,
        'pole_lat': plat, 'pole_lon': plon, 'pole_dp': dp, 'pole_dm': dm, 'pole_n_sites': n,
        'external_database_ids': f'GPMDB[{GPMDB[unit]}]', 'description': authors,
        '_unit': unit})
locs = pd.DataFrame(loc_rows)

AGE_ROWS = [
    ('FG', None, 1047, 5.0, None, None, 'GM-UPB', MCLELLAND2001, 'Stenian',
     'U-Pb zircon, SHRIMP weighted mean of the zircon rims, Wanakena fayalite granite (sample '
     'CLFG; the pluton of sites AD9, AD10, AD20, AD37, AD38, ADH4 and ADH5): 1047 +/- 10 Ma '
     '(95% confidence; McLelland et al., 2001).'),
    ('FG', 'ADH20', 1047.3, 1.1, None, None, 'GM-UPB', MCLELLAND2001, 'Stenian',
     'U-Pb zircon, single-grain TIMS upper intercept, fayalite granite near Au Sable Forks '
     '(sample AM-86-11): 1047.3 +/- 2.2 Ma (95% confidence; McLelland et al., 2001).'),
    ('ANT', None, 1158.90, 0.285, None, None, 'GM-UPB', ZHANG2026, 'Stenian',
     'U-Pb zircon CA-ID-TIMS, concordant igneous zircon MA1-z6 of Marcy massif anorthositic '
     'sample MA1: 1158.90 +/- 0.57 Ma (2 sigma; Zhang et al., 2026), dating crystallization of '
     'the massif ca. 1159 Ma.'),
    ('ANT', None, 1158.70, 0.46, None, None, 'GM-UPB', ZHANG2026, 'Stenian',
     'U-Pb zircon CA-ID-TIMS, concordant igneous zircon fragment MA1-z10a of Marcy massif '
     'anorthositic sample MA1: 1158.70 +/- 0.92 Ma (2 sigma; Zhang et al., 2026).'),
    ('ANT', None, 1159.73, 0.46, None, None, 'GM-UPB', ZHANG2026, 'Stenian',
     'U-Pb zircon CA-ID-TIMS, concordant igneous zircon fragment MA1-z10b of Marcy massif '
     'anorthositic sample MA1: 1159.73 +/- 0.92 Ma (2 sigma; Zhang et al., 2026).'),
    ('ANT', None, 995.81, 0.37, None, None, 'GM-UPB', ZHANG2026, 'Tonian',
     'U-Pb zircon CA-ID-TIMS, youngest metamorphic zircon of Marcy massif anorthositic sample '
     'MA5: 995.81 +/- 0.74 Ma (2 sigma; Zhang et al., 2026), dating protracted high-grade '
     'metamorphism.'),
    ('ANT', None, None, None, 850.11, 890.08, 'GM-UPB', ZHANG2026, 'Tonian',
     'U-Pb rutile ID-TIMS 206Pb/238U cooling dates of six grains, Marcy massif sample MA1: '
     '890.08 +/- 0.50 to 850.11 +/- 0.64 Ma (2 sigma; Zhang et al., 2026), correlated with '
     'grain size (volume diffusion).'),
    ('ANT', None, None, None, 851.45, 881.78, 'GM-UPB', ZHANG2026, 'Tonian',
     'U-Pb apatite ID-TIMS 206Pb/238U cooling dates of 22 grains, Marcy massif samples MA1 and '
     'MA5: 881.78 +/- 1.47 to 851.45 +/- 0.90 Ma (2 sigma; Zhang et al., 2026), correlated with '
     'grain size (volume diffusion).'),
]
ages = pd.DataFrame([{'site': s, 'location': None if s else LOC[u], 'age': a, 'age_sigma': sg,
                      'age_low': lo, 'age_high': hi, 'age_unit': 'Ma', 'method_codes': mc,
                      'citations': cit, 'timescale_eon': 'Proterozoic',
                      'timescale_era': 'Mesoproterozoic' if per == 'Stenian' else 'Neoproterozoic',
                      'timescale_period': per, 'description': desc, '_unit': u}
                     for u, s, a, sg, lo, hi, mc, cit, per, desc in AGE_ROWS])


def fmt(df):
    df = df.drop(columns='_unit').dropna(axis=1, how='all').copy()
    for c in df.columns:
        if df[c].dtype.kind == 'f':
            df[c] = df[c].map(lambda x: '' if pd.isna(x) else f'{x:.10g}')
    return df.fillna('')


def table_text(df, kind):
    return f'tab delimited\t{kind}\n' + fmt(df).to_csv(sep='\t', index=False)


TABLES = (('locations', locs), ('sites', sites), ('ages', ages))
for kind, df in TABLES:
    with open(os.path.join(OUTDIR, f'{kind}.txt'), 'w', encoding='utf-8') as f:
        f.write(table_text(df, kind))

stamp = date.today().strftime('%d.%b.%Y')
upload = os.path.join(OUTDIR, f'Brown2012_Adirondack_Highlands_{stamp}.txt')
with open(upload, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {upload}')

for unit, (p, d, sub) in poles.items():
    print(f'  {LOC[unit]:36s} B {int(p["n"]):2d} N {int(sub["dir_n_samples"].sum()):2d} '
          f'pole {p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} K {p["k"]:.1f} | dir '
          f'{d["dec"]:.1f}/{d["inc"]:.1f} k {d["k"]:.1f} a95 {d["alpha95"]:.1f}')

# the compilation pole (lat, lon, A95, B, N) and the Fisher means of the tabulated site VGPs
# quoted in the descriptions of the authors' gneiss and granite poles
EXPECT = {'ANT': (-27.5, 143.5, 12.5, 14, 68), 'GMS': (-19.0, 148.7, 11.2, 14, 80),
          'FG': (-28.5, 131.7, 7.1, 8, 40)}
for unit, (p, d, sub) in poles.items():
    got = (round(p['inc'], 1), round(p['dec'], 1), round(p['alpha95'], 1), int(p['n']),
           int(sub['dir_n_samples'].sum()))
    assert got == EXPECT[unit], (unit, got)
print('-I- site-VGP means match the expected values')

sys.path.insert(0, os.path.join(REPO, 'scripts'))
from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True)
