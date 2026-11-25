"""Build the MagIC contribution for the Franklin LIP poles of Denyszyn et al. (2009), including the Franklin event grand mean.

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Denyszyn, S. W., Halls, H. C., Davis, D. W., & Evans, D. A. D. (2009). Paleomagnetism and
U-Pb geochronology of Franklin dykes in High Arctic Canada and Greenland: a revised age and
paleomagnetic pole constraining block rotations in the Nares Strait region. Canadian Journal
of Earth Sciences 46, 689-705. doi:10.1139/E09-042

Previous studies whose site results make up the grand mean (each site row cites its own
study): Fahrig, Irving & Jackson 1971 (10.1139/e71-047); Fahrig & Schwarz 1973
(10.1139/e73-057); Park 1974 (10.1139/e74-098); Palmer & Hayatsu 1975 (10.1139/e75-130);
Park 1981 (10.1139/e81-151); Christie & Fahrig 1983 (10.1139/e83-025); Palmer, Baragar,
Fortier & Foster 1983 (10.1139/e83-131).
Also named in descriptions: Park 1973 (10.1139/e73-093; Daly Bay host rocks); Wanless,
Stevens, Lachance & Rimsaite 1966 (Geological Survey of Canada Paper 65-17, K-Ar isotopic
ages, report 6, 10.4095/100990; location of Fahrig et al. 1971 site 29); Roest & Srivastava 1989
(10.1130/0091-7613(1989)017<1000:SFSITL>2.3.CO;2; Greenland restoration); Evans et al.
2021 (10.1016/B978-0-12-818533-9.00007-2; compilation that uses D+09can).
Dates: Pu, J. P., et al. (2022). Emplacement of the Franklin large igneous province and
initiation of the Sturtian Snowball Earth. Science Advances 8, eadc9430.
doi:10.1126/sciadv.adc9430 (U-Pb zircon and baddeleyite); Dufour, F., et al. (2023). New
U-Pb CA-ID TIMS zircon ages implicate the Franklin LIP as the proximal trigger for the
Sturtian Snowball Earth event. Earth and Planetary Science Letters 618, 118259.
doi:10.1016/j.epsl.2023.118259 (U-Pb zircon); Denyszyn et al. (2009) Table 2 (U-Pb
baddeleyite). Uncertainties are 2 sigma; age_sigma is 1 sigma.

Sources and checks
------------------
- Denyszyn et al. (2009) Table 1 (p. 694; 27 sites) is transcribed below from the page
  image. Coordinates are printed as degrees and decimal minutes (75 43.1 = 75.7183 N).
  Reversed sites (SG2, NU1, HE) are listed with the antipole of the VGP; the VGP of the
  observed direction is stored. Daggered sites fail the authors' filter (N >= 4 samples
  and alpha95 <= 15 deg) and are entered with result_quality 'b'.
- The previous-study sites are read from Franklin_previous_studies_sites.csv, assembled by
  assemble_previous_studies.py from each paper's site table, transcribed from the page
  images (source_transcriptions/paper_*.csv): Fahrig et al. 1971 Table 2 (p. 459), Fahrig &
  Schwarz 1973 Table 1 (p. 578), Park 1974 Table 1 (p. 1014), Palmer & Hayatsu 1975
  Table 1 (p. 1442), Park 1981 Table 1 (p. 1640), Christie & Fahrig 1983 Table 1
  (pp. 280-281), Palmer et al. 1983 Tables 2-4 (pp. 1464-1467). Denyszyn et al. (2009) do
  not tabulate these sites; which sites enter the grand mean follows the site list of
  their compilation (D. A. D. Evans, compilation workbook; source_transcriptions/evans_*.csv;
  the unit VGPs are in ../Denyszyn2009.csv). The 78 units match one to one; the VGP
  differences are listed unit by unit in ../Denyszyn2009_updated.csv, written by
  make_denyszyn2009_updated.py. Differences in the workbook that do not enter any pole: Park
  (1981) counts are specimens with recalculated alpha95; Fahrig et al. (1971) site 1 is
  entered as 260/-1 (the printed direction is 080/-01); and the k values of nine of the
  Greenland sites differ from Table 1 of Denyszyn et al. (2009), which is followed here.
- Accepted units per study: Fahrig et al. 1971, 19 (5 tabulated Franklin sites fail the
  filter); Fahrig & Schwarz 1973, 8 (2); Park 1974, 2 (1, Ontario dyke 16); Palmer &
  Hayatsu 1975, 3 (13); Park 1981, 3 (2); Christie & Fahrig 1983, 7 sites plus 2 dyke means
  of 6 sites on two dykes; Palmer et al. 1983, 12 (2). Sites that fail the filter are
  entered with result_quality 'b'. Sites the original authors set aside (Fahrig et al.
  1971 Group B; Palmer & Hayatsu 1L, 9L; Palmer et al. 1983 lower-reverse and steep zones
  and dolerites; Christie & Fahrig Borden, Tremblay Sound and overprint groups; Park 1974
  Mackenzie dykes) are not entered. Robertson & Baragar (1972) and Jones & Fahrig (1978)
  contribute no site that passes the filter.
- Directions are stored as printed, in their observed polarity. Site VGPs are calculated
  from the directions with pmag.dia_vgp, including those of Fahrig et al. 1971 and Park
  1974, which print site coordinates and VGPs in whole degrees (the printed VGPs are quoted
  in the site descriptions and asserted to agree within 1.6 deg). Palmer & Hayatsu 1975 and
  Park 1981 give no site coordinates, so their tabulated VGPs are used, oriented to the
  observed direction. With the printed VGPs (and the printed location of Fahrig et al. 1971
  site 29), D+09can and D+09cg are the published poles (checked below).
- Coordinates where a paper gives none: Fahrig & Schwarz 1973, Christie & Fahrig 1983 (the
  values agree with its Fig. 1 to ca. 0.1 deg) and Palmer et al. 1983 use the coordinates
  of the Denyszyn et al. (2009) compilation; Palmer & Hayatsu 1975 and Park 1981 use area
  coordinates back-calculated from the tabulated VGPs.
- Points noted in the source tables (values stored as printed unless stated):
  Fahrig et al. 1971 site 1 prints 080/-01 (stored; reversed, so the inverted direction is
  260/+01). Site 29 is printed at 54.5 N, 117.0 W, and its tabulated VGP (13/153) follows
  from that latitude; the site is the NE-trending dike near Hardisty Lake (text, p. 457),
  and Table 1 ties it to K-Ar determination GSC 64-64, whose sample (FA-162-62) is at
  64 deg 27' N, 117 deg 00' W (Wanless et al., 1966, Geological Survey of Canada Paper
  65-17, p. 46). That latitude is stored, with the VGP calculated from the direction
  (12.7/150.3). Site 20 (301/-7,
  westerly, here normal) is labelled R on the Fahrig & Schwarz (1973) Fig. 1 map. Park
  1974 sites 31 and 33 are means of dyke and contact-rock samples; site 16, a dyke in southern Ontario, passes the quality filter but
  lies outside the region of Franklin magmatism (Denyszyn et al., 2009, p. 699). Park 1981
  directions are corrected for bedding (tilt correction 100). Its alpha95 values follow
  from the specimen count and R, except site 1 (printed 9, whereas N = 4, R = 3.91 give
  16.1); sample and specimen counts and the printed alpha95 are stored, and the filter is
  applied, as in the compilation, to specimen counts and recalculated alpha95. Fahrig &
  Schwarz 1973 and Palmer et al. 1983 give no alpha95 (it is calculated from k and N).
  Palmer et al. 1983 directions are the anchored-line means (their method 2); the method-1
  and averaged directions are in the site descriptions. Christie & Fahrig 1983 site 46
  prints alpha95 8 with k 81 and N 6.

Structure and interpretation
----------------------------
- Five locations, each a Fisher mean of polarity-unified site VGPs reported at the Pacific
  pole:
    1. Devon and Ellesmere islands, 12 dykes (5.8/189.7, A95 10.7, as published), with the
       positive baked-contact test at GR (contact_test C+). GPMDB 9502.
    2. northwest Greenland, 10 dykes and sills (8.8/178.5, A95 7.2). GPMDB 9503.
    3. both together, 22 sites with no internal rotation (7.2/184.5, A95 6.7).
    4. the Franklin event grand mean, previous Canadian studies (pole D+09can of Table 3):
       56 units of seven studies, none of them Denyszyn et al.'s own sites. This is the pole
       used in the Swanson-Hysell et al. (2026) compilation (6.6/162.1, A95 3.0; the
       published 6.7/162.1 with the printed VGPs and the printed location of Fahrig et al.
       1971 site 29).
       contact_test C+: Denyszyn et al. (2009, pp. 698-699 and Table 3) regard the combined
       Park (1973) host-rock and Park (1974) dyke data of the Daly Bay area (sites 31 and
       33) as a positive baked-contact test. D+09can is preferred to D+09cg (as in Evans
       et al., 2021) because it rests only on autochthonous sites; the position of D+09cg
       depends on the restoration rotations, and the microplate rotation is inferred in
       part from the Devon and Ellesmere data themselves.
    5. the grand mean with Greenland and the proposed Ellesmere microplate restored to
       North America (pole D+09cg, 78 sites; the pole of the abstract). The site VGPs of
       Denyszyn et al. are rotated about 72 N, 274 E by -20 deg (Devon and Ellesmere; the
       microplate rotation of the paper) and about 67.5 N, 241.5 E by -13.8 deg (Greenland;
       Roest & Srivastava, 1989, the parameters of the compilation workbook), which
       reproduces the published 8.4/163.8, A95 2.8. The restored VGPs in ../Denyszyn2009.csv
       use -35 deg for the microplate and give 9.0/161.7.
  The Devon-only and Ellesmere-only subset poles of the text and the individual poles of
  the earlier studies (Table 3) are not entered as poles.
- The grand-mean rows have no mean direction because their sites span ca. 2500 km, and so
  no dir_tilt_correction.
- Two groups of Christie & Fahrig (1983) sites lie on the same dyke; the grand mean uses
  their sample-weighted dyke means, which are entered as result_type 'a' site rows, with
  the constituent sites as 'i' rows.
- Site and location ages are the Franklin emplacement age of the compilation, 719 Ma
  (718-720 Ma), from the zircon dates of Pu et al. (2022). The ages table holds the
  measured dates of Pu et al. (2022), Dufour et al. (2023) and Denyszyn et al. (2009).
"""
import os, sys
from datetime import date
import numpy as np
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

CIT_D09 = '10.1139/E09-042'
CIT_PU22 = '10.1126/sciadv.adc9430'
CIT_E21 = '10.1016/B978-0-12-818533-9.00007-2'   # Evans et al. (2021), compilation using D+09can
AGE, AGE_LOW, AGE_HIGH = 719, 718, 720
GEOG = dict(continent_ocean='North America', terranes='Laurentia')

# ---- Denyszyn et al. (2009) Table 1 (p. 694) ------------------------------------
# site, flags, lat deg, lat min, lon(W) deg, lon min, strike, width (m), NS, D, I,
# Plat, Plon (antipole for reversed sites), k, a95, area
# flags: '*' magnetic compass only; 'R' reversed; 'x' excluded by the quality filter
# (dagger in Table 1); 'S' sill.
T1 = [
    ('BB', '', 75, 43.1, 82, 59.2, 285, 20, 9, 286.0, -18.0, -5.1, 169.3, 91, 5.4, 'Devon'),
    ('PH', 'x', 75, 0.0, 79, 36.0, 255, 30, 5, 282.1, 12.6, 9.3, 180.3, 8, 28.3, 'Devon'),
    ('GR', '', 76, 25.9, 83, 1.2, 315, 40, 8, 272.9, 23.3, 12.5, 187.1, 190, 4.0, 'Ellesmere'),
    ('BG', '', 75, 35.7, 80, 11.7, 270, 75, 7, 271.3, 34.5, 18.7, 193.4, 24, 12.6, 'Devon'),
    ('BP', '', 75, 46.1, 81, 19.4, 265, 34, 9, 279.2, -7.5, -1.4, 178.8, 57, 6.9, 'Devon'),
    ('BR', '', 75, 28.2, 81, 33.1, 275, 30, 7, 282.6, -7.5, -0.5, 175.3, 226, 4.0, 'Devon'),
    ('EG', '', 75, 49.4, 82, 4.2, 255, 21, 7, 280.1, 24.8, 15.1, 181.4, 134, 5.2, 'Devon'),
    ('LG', '', 78, 43.6, 75, 41.0, 260, 50, 7, 243.1, 13.7, 1.8, 222.0, 91, 6.4, 'Ellesmere'),
    ('CG', '', 78, 17.1, 77, 6.6, 280, 25, 7, 251.9, 1.8, -2.7, 210.8, 46, 9.0, 'Ellesmere'),
    ('CM', '*', 74, 38.0, 80, 28.0, 220, 25, 7, 282.5, 1.7, 4.2, 177.7, 55, 9.1, 'Devon'),
    ('EA', '*', 75, 49.4, 82, 4.7, 260, 40, 11, 273.4, 31.5, 17.3, 188.9, 79, 5.5, 'Devon'),
    ('HM', '*', 74, 38.3, 80, 12.6, 260, 55, 4, 244.2, -4.7, -8.9, 214.2, 240, 5.9, 'Devon'),
    ('PI', '*x', 75, 0.4, 79, 46.2, 225, 50, 5, 293.3, 35.8, 25.1, 172.5, 22, 16.7, 'Devon'),
    ('CW', '*x', 74, 28.6, 81, 51.9, 305, 30, 5, 302.1, -0.2, 8.2, 157.0, 6, 35.4, 'Devon'),
    ('SG2', '*R', 75, 42.0, 84, 0.0, 305, 20, 6, 101.2, -26.4, 16.3, 178.6, 103, 6.6, 'Devon'),
    ('KL', '', 78, 41.0, 70, 40.0, 270, 54, 8, 283.0, 10.9, 7.9, 187.6, 74, 6.5, 'Greenland'),
    ('TB', '', 76, 27.5, 69, 14.4, 295, 30, 9, 296.5, 26.4, 19.6, 178.0, 37, 8.6, 'Greenland'),
    ('NU2', '', 77, 22.7, 71, 29.0, 295, 40, 8, 289.3, 16.3, 12.3, 181.4, 123, 5.0, 'Greenland'),
    ('PK', '', 77, 56.3, 72, 12.5, 280, 65, 7, 288.7, -18.5, -5.5, 177.6, 132, 5.3, 'Greenland'),
    ('KC', '', 77, 33.9, 70, 15.4, 300, 16, 4, 290.6, 15.7, 12.2, 181.2, 327, 5.1, 'Greenland'),
    ('PW', '', 77, 11.9, 70, 50.1, 290, 10, 7, 293.6, 15.9, 13.0, 177.8, 144, 5.1, 'Greenland'),
    ('QA', '', 77, 29.1, 68, 51.9, 305, 33, 6, 298.8, 26.1, 19.5, 175.7, 106, 6.5, 'Greenland'),
    ('KA', '*x', 77, 12.3, 70, 46.6, 300, 55, 2, 307.8, -1.6, 7.0, 162.0, 435, 12.0, 'Greenland'),
    ('CA', 'S', 76, 44.3, 73, 13.4, None, 20, 5, 295.4, 9.3, 10.2, 173.0, 234, 5.0, 'Greenland'),
    ('GF', 'S', 76, 51.4, 69, 55.5, None, 10, 6, 295.6, 17.3, 14.4, 177.0, 68, 8.2, 'Greenland'),
    ('NU1', 'R', 77, 24.5, 71, 33.7, 290, 40, 8, 108.8, 37.7, -16.6, 175.6, 59, 7.3, 'Greenland'),
    ('HE', 'Rx', 77, 25.3, 70, 8.5, 280, 21, 5, 118.3, 11.9, 0.0, 171.0, 14, 21.5, 'Greenland'),
]


# ---- helpers ---------------------------------------------------------------------
def r1(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 1)


def dm2dd(d, m):
    return round(d + m / 60.0, 4)


def angle(a, b):
    """Angle (deg) between two (lon, lat) poles."""
    return pmag.angle([a[0], a[1]], [b[0], b[1]])[0]


def observed_vgp(dec, inc, a95, lat, lon, paper=None, tol=1.6, label=''):
    """VGP of the observed direction. If the paper tabulates a VGP (in whatever
    polarity convention), use it -- oriented to the observed-direction VGP -- after
    checking that it reproduces within tol degrees."""
    plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95 if a95 and not np.isnan(a95) else 0., lat, lon)
    if paper is not None and not any(np.isnan(v) for v in paper):
        plat_p, plon_p = paper
        if angle((plon, plat), (plon_p, plat_p)) > 90:           # tabulated as antipole
            plat_p, plon_p = -plat_p, (plon_p + 180) % 360
        miss = angle((plon, plat), (plon_p, plat_p))
        assert miss < tol, (label, miss, (plat, plon), paper)
        plat, plon = plat_p, plon_p
    ok = a95 is not None and not np.isnan(a95)
    return r1(plat), r1(plon % 360), (r1(dp) if ok else None), (r1(dm) if ok else None)


site_rows, loc_rows = [], []

# ---- Denyszyn et al. (2009) sites -------------------------------------------------
D09_LOC = {'Devon': 'Franklin dykes, Devon and Ellesmere islands',
           'Ellesmere': 'Franklin dykes, Devon and Ellesmere islands',
           'Greenland': 'Franklin dykes and sills, northwest Greenland'}
for (site, fl, la, lam, lo, lom, strike, width, ns, dec, inc, plat, plon, k, a95, area) in T1:
    lat, lon = dm2dd(la, lam), round(360 - (lo + lom / 60.0), 4)
    vlat, vlon, dp, dm = observed_vgp(dec, inc, a95, lat, lon, (plat, plon), label=site)
    sill = 'S' in fl
    codes = ['SO-MAG'] if '*' in fl else ['SO-SUN', 'SO-MAG']
    codes += ['LP-DIR-AF'] + (['LP-DIR-T'] if site == 'GR' else []) + ['DE-BFL', 'DE-FM']
    island = {'Devon': 'on Devon Island', 'Ellesmere': 'on Ellesmere Island',
              'Greenland': 'in the Thule region, northwest Greenland'}[area]
    desc = (f"Franklin {'sill' if sill else 'dyke'} {island}"
            + ('' if sill else f", strike {strike} deg") + f", width {width} m.")
    if 'R' in fl:
        desc += ' Reversed (easterly) polarity; the paper lists the antipole of the VGP.'
    if '*' in fl:
        desc += ' Samples oriented by magnetic compass only.'
    if site == 'GR':
        desc += (' Site of the positive baked-contact test into ca. 1.97 Ga anorthosite '
                 '(Denyszyn et al., 2006, 2009).')
    if 'x' in fl:
        desc += (" Excluded from the means by the authors' quality filter (N >= 4 and "
                 'alpha95 <= 15 deg required).')
    site_rows.append({
        'site': site, 'location': D09_LOC[area], 'result_type': 'i',
        'result_quality': 'b' if 'x' in fl else 'g', 'method_codes': ':'.join(codes),
        'citations': CIT_D09, 'geologic_classes': 'Intrusive',
        'geologic_types': 'Sill' if sill else 'Volcanic Dike', 'lithologies': 'Diabase',
        'lat': lat, 'lon': lon, 'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH,
        'age_unit': 'Ma', 'dir_tilt_correction': 0, 'dir_dec': dec, 'dir_inc': inc,
        'dir_alpha95': a95, 'dir_k': k, 'dir_n_samples': ns,
        'vgp_lat': vlat, 'vgp_lon': vlon, 'vgp_dp': dp, 'vgp_dm': dm,
        'description': desc})

# ---- previous Canadian studies (D+09can inputs) ------------------------------------
PRIOR = pd.read_csv(os.path.join(HERE, 'Franklin_previous_studies_sites.csv'),
                    dtype={'site_printed': str})
GM = 'Franklin event grand mean, previous Canadian studies'
STUDY = {
    'Fahrig1971': 'Fahrig et al. (1971)', 'FahrigSchwarz1973': 'Fahrig & Schwarz (1973)',
    'Park1974': 'Park (1974)', 'PalmerHayatsu1975': 'Palmer & Hayatsu (1975)',
    'Park1981': 'Park (1981)', 'ChristieFahrig1983': 'Christie & Fahrig (1983)',
    'Palmer1983': 'Palmer et al. (1983)'}
MCODES = {'AF blanket': 'LP-DIR-AF:DE-BLANKET:DE-FM',
          'AF+TH PCA': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM',
          'anchored PCA': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-BFL-A:DE-FM'}
SO = {'Fahrig1971': 'SO-SUN:', 'FahrigSchwarz1973': 'SO-SUN:', 'ChristieFahrig1983': 'SO-SUN:',
      'Palmer1983': 'SO-SUN:SO-MAG:SO-SIGHT:', 'PalmerHayatsu1975': 'SO-SUN:SO-SIGHT:'}
# studies that print site coordinates: VGPs are calculated from the directions (the printed
# VGPs are in whole degrees and are checked against the calculated ones)
CALC_VGP = ('Fahrig1971', 'Park1974')
for r in PRIOR.itertuples():
    pol_paper = (r.vgp_lat_paper, r.vgp_lon_paper)
    printed = not any(np.isnan(v) for v in pol_paper)
    if r.study in CALC_VGP:
        if printed:
            observed_vgp(r.dec, r.inc, r.a95, r.lat, r.lon, pol_paper, label=r.site)
        vlat, vlon, dp, dm = observed_vgp(r.dec, r.inc, r.a95, r.lat, r.lon)
    else:
        vlat, vlon, dp, dm = observed_vgp(r.dec, r.inc, r.a95, r.lat, r.lon, pol_paper,
                                          label=r.site)
    lith_class = 'Extrusive' if r.geologic_types == 'Lava Flow' else 'Intrusive'
    desc = f"{r.group}; site {r.site_printed.replace(' Mean', '')} of {STUDY[r.study]}; {r.treatment}."
    if isinstance(r.note, str) and r.note:
        desc += ' ' + r.note
    if r.study in CALC_VGP and printed:
        desc += (' The VGP is calculated from the direction at the tabulated coordinates; '
                 f"{STUDY[r.study]} {'tabulate' if 'et al.' in STUDY[r.study] else 'tabulates'} "
                 f'a VGP of {r.vgp_lat_paper:g}, '
                 f'{r.vgp_lon_paper:g} (latitude, longitude; whole degrees).')
    desc += ' Site coordinates: ' + r.coord_source + '.'
    if r.d09 == 'rejected':
        if isinstance(r.reject_reason, str) and r.reject_reason:
            desc += ' ' + r.reject_reason
        else:
            # the authors' filter: N >= 4 and alpha95 <= 15 deg (alpha95 from k and N
            # where the paper gives none)
            n = r.n_samples
            R = n - (n - 1) / r.k
            a95 = r.a95 if not np.isnan(r.a95) else np.degrees(np.arccos(
                1 - (n - R) / R * (20 ** (1 / (n - 1)) - 1)))
            assert r.n_samples < 4 or a95 > 15, (r.site, r.n_samples, a95)
            desc += (' Not used in the grand mean of Denyszyn et al. (2009), which requires '
                     'N >= 4 and alpha95 <= 15 deg.')
    if r.d09 == 'constituent':
        desc += (f' Averaged with the other sites on the same dyke ({r.unit}) for the '
                 'grand mean; see the dyke-mean row.')
    site_rows.append({
        'site': r.site, 'location': GM, 'result_type': 'i',
        'result_quality': 'b' if r.d09 == 'rejected' else 'g',
        'method_codes': SO.get(r.study, '') + MCODES[r.method], 'citations': r.citation,
        'geologic_classes': lith_class, 'geologic_types': r.geologic_types,
        'lithologies': r.lithologies, 'lat': r.lat, 'lon': r.lon,
        'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': 0 if np.isnan(r.tilt_correction) else int(r.tilt_correction),
        'dir_dec': r.dec, 'dir_inc': r.inc, 'dir_alpha95': r1(r.a95), 'dir_k': r1(r.k),
        'dir_n_samples': r.n_samples,
        'dir_n_total_samples': None if np.isnan(r.n_total_samples) else int(r.n_total_samples),
        'dir_n_specimens': None if np.isnan(r.n_specimens) else int(r.n_specimens),
        'vgp_lat': vlat, 'vgp_lon': vlon, 'vgp_dp': dp, 'vgp_dm': dm,
        'description': desc})

# Christie & Fahrig (1983) dyke units combined by Denyszyn et al. (2009):
# sample-weighted Fisher mean of the site means on one dyke
for unit, sub in PRIOR[PRIOR.d09 == 'constituent'].groupby('unit'):
    blk = [[d, i, 1.] for d, i, n in zip(sub.dec, sub.inc, sub.n_samples) for _ in range(int(n))]
    f = pmag.fisher_mean(blk)
    lat, lon = sub.lat.iloc[0], sub.lon.iloc[0]
    vlat, vlon, _, _ = observed_vgp(f['dec'], f['inc'], None, lat, lon)
    name = 'CF83-' + unit.replace('/', '-')
    site_rows.append({
        'site': name, 'location': GM, 'result_type': 'a', 'result_quality': 'g',
        'method_codes': 'SO-SUN:LP-DIR-AF:DE-BLANKET:DE-FM', 'citations': '10.1139/e83-025',
        'geologic_classes': 'Intrusive', 'geologic_types': 'Volcanic Dike',
        'lithologies': 'Diabase', 'lat': lat, 'lon': lon,
        'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': 0, 'dir_dec': r1(f['dec']), 'dir_inc': r1(f['inc']),
        'dir_n_samples': int(sub.n_samples.sum()), 'vgp_lat': vlat, 'vgp_lon': vlon,
        'description': (
            f"Dyke mean of Christie & Fahrig (1983) sites {unit.replace('/', ', ')}, which lie "
            "on the same Franklin dyke, as combined for the grand mean of Denyszyn et al. "
            "(2009): Fisher mean of the site means weighted by the number of samples. "
            "Not a published result of Christie & Fahrig (1983).")})

sites = pd.DataFrame(site_rows)
for c in ['dir_n_samples', 'dir_n_total_samples', 'dir_n_specimens']:
    sites[c] = pd.to_numeric(sites[c]).astype('Int64')


# ---- poles --------------------------------------------------------------------------
def minority_mask(sub):
    d = ipmag.fisher_mean(di_block=pmag.flip(
        ipmag.make_di_block(sub['dir_dec'].tolist(), sub['dir_inc'].tolist()), combine=True))
    ref = pmag.dir2cart([d['dec'], d['inc'], 1.])
    away = np.array([np.dot(pmag.dir2cart([a, b, 1.]), ref) < 0
                     for a, b in zip(sub['dir_dec'], sub['dir_inc'])])
    return away if away.sum() <= len(away) / 2 else ~away


def poles(sub):
    blk = pmag.flip(ipmag.make_di_block(sub['vgp_lon'].tolist(), sub['vgp_lat'].tolist()),
                    combine=True)
    p = ipmag.fisher_mean(di_block=blk)
    if p['inc'] < 0 or p['dec'] < 90:          # report the Pacific (northern) pole
        p = ipmag.fisher_mean(di_block=ipmag.do_flip(di_block=blk))
    m = minority_mask(sub)
    decs = np.where(m, (sub['dir_dec'] + 180) % 360, sub['dir_dec'])
    incs = np.where(m, -sub['dir_inc'], sub['dir_inc'])
    d = ipmag.fisher_mean(di_block=ipmag.make_di_block(list(decs), list(incs)))
    return p, d, round(100 * m.sum() / len(sub), 1)


good = sites[sites['result_quality'] == 'g']
d09 = good[good['citations'] == CIT_D09]
canada = d09[d09['location'] == D09_LOC['Devon']]
green = d09[d09['location'] == D09_LOC['Greenland']]
gm = good[(good['location'] == GM) & ~good['site'].isin(
    PRIOR.loc[PRIOR.d09 == 'constituent', 'site'])]
assert len(canada) == 12 and len(green) == 10 and len(gm) == 56


def bbox(sub):
    return dict(lat_s=round(sub['lat'].min(), 4), lat_n=round(sub['lat'].max(), 4),
                lon_w=round(sub['lon'].min(), 4), lon_e=round(sub['lon'].max(), 4))


def union_codes(sub, extra):
    out = []
    for m in sub['method_codes']:
        for c in m.split(':'):
            if c not in out:
                out.append(c)
    return ':'.join(out + [e for e in extra if e not in out])


def loc_row(name, rname, sub, cit, desc, with_dir=True, extra_codes=('DE-VGP',), **kw):
    p, d, rev = poles(sub)
    row = {'location': name, 'location_type': 'Region', 'result_name': rname,
           'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(sub['site']),
           'method_codes': union_codes(sub, list(extra_codes)), 'citations': cit,
           'geologic_classes': ':'.join(sorted(set(':'.join(sub['geologic_classes']).split(':')))),
           'lithologies': ':'.join(sorted(set(':'.join(sub['lithologies']).split(':')))),
           **bbox(sub), 'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH,
           'age_unit': 'Ma'}
    if with_dir:
        row.update({'dir_tilt_correction': 0, 'dir_dec': r1(d['dec']), 'dir_inc': r1(d['inc']),
                    'dir_alpha95': r1(d['alpha95']), 'dir_k': r1(d['k']),
                    'dir_n_sites': int(d['n'])})
    row.update({'pole_lat': r1(p['inc']), 'pole_lon': r1(p['dec']),
                'pole_alpha95': r1(p['alpha95']), 'pole_k': r1(p['k']),
                'pole_n_sites': int(p['n']), 'pole_reversed_perc': rev, **GEOG, **kw,
                'description': desc})
    return row


D09_METHODS = (
    "Samples were drilled or taken as oriented blocks, oriented by sun and magnetic compass "
    "(magnetic compass only at seven sites), demagnetized by detailed alternating field and, "
    "for selected specimens, thermal methods at the University of Toronto, and analysed by "
    "principal components (MAD <= 10 deg). A stable, high-coercivity, high-unblocking-"
    "temperature component carried by near-stoichiometric magnetite is interpreted as "
    "primary. Sites with fewer than four samples or alpha95 above 15 deg are excluded. The "
    "authors list the VGPs of reversed sites as antipoles; they are stored here as the VGPs of "
    "the observed directions.")
loc_rows.append(loc_row(
    D09_LOC['Devon'], 'Franklin dykes of Devon and Ellesmere islands ca. 719 Ma pole', canada,
    f'{CIT_D09}:{CIT_PU22}',
    ("Paleomagnetic pole for east-west-trending Franklin dykes of Devon Island (9 sites) and "
     "southern Ellesmere Island (3 sites), Arctic Canada (Denyszyn et al., 2009): Fisher mean "
     "of the in-situ site VGPs, in present-day coordinates. The authors note that the "
     "Ellesmere sites alone give a pole far from the Devon pole (statistically indistinct "
     "because of the small number of sites) and interpret the Ellesmere microplate as rotated "
     "by ca. 20 deg in the early Cenozoic. A baked-contact test at site GR (Ellesmere Island) "
     "into Paleoproterozoic anorthosite is positive: the dyke direction is found in the baked "
     "host and the distinct host direction in unbaked host rock, without hybrid "
     "magnetizations. " + D09_METHODS),
    extra_codes=('DE-VGP', 'ST-C'), contact_test='C+', country='Canada',
    state_province='Nunavut', external_database_ids='GPMDB[9502]'))
loc_rows.append(loc_row(
    D09_LOC['Greenland'], 'Franklin dykes and sills of northwest Greenland ca. 719 Ma pole',
    green, f'{CIT_D09}:{CIT_PU22}',
    ("Paleomagnetic pole for Franklin dykes and two Thule sills of the Thule region, "
     "northwest Greenland (Denyszyn et al., 2009): Fisher mean of the in-situ site VGPs, in "
     "present-day Greenland coordinates (not restored to North America). " + D09_METHODS),
    continent_ocean='Greenland', country='Greenland', external_database_ids='GPMDB[9503]'))
comb = pd.concat([canada, green])
loc_rows.append(loc_row(
    'Franklin dykes, Arctic Canada and northwest Greenland',
    'Franklin dykes of Arctic Canada and Greenland ca. 719 Ma pole (no internal rotations)',
    comb, f'{CIT_D09}:{CIT_PU22}',
    ("Paleomagnetic pole for all quality-filtered Franklin dyke and sill sites of Denyszyn et "
     "al. (2009) from Devon and Ellesmere islands and northwest Greenland, with no internal "
     "rotations applied (Fisher mean of in-situ site VGPs in present-day coordinates). The "
     "authors regard this and their subgroup poles as distinct from earlier Franklin results "
     "and attribute the difference to early Cenozoic rotation of Greenland and the Ellesmere "
     "microplate relative to North America. " + D09_METHODS),
    extra_codes=('DE-VGP', 'ST-C'), contact_test='C+'))
loc_rows.append(loc_row(
    GM, 'Franklin event grand mean ca. 719 Ma pole (previous Canadian studies)', gm,
    ':'.join([CIT_D09] + sorted(set(gm['citations'])) + [CIT_PU22, CIT_E21]),
    ("Grand-mean paleomagnetic pole of the Franklin igneous event from the quality-filtered "
     "site results of previous studies in autochthonous Canada (Victoria Island lavas and "
     "sills, Brock Inlier and Coronation sills, and dykes of the mainland, Melville "
     "Peninsula, Boothia and Baffin Island), as compiled by Denyszyn et al. (2009, their pole "
     "D+09can): Fahrig et al. (1971), Fahrig & Schwarz (1973), Park (1974, 1981), Palmer & "
     "Hayatsu (1975), Christie & Fahrig (1983) and Palmer et al. (1983). Sites with fewer than "
     "four samples or alpha95 above 15 deg are excluded, and two groups of Christie & "
     "Fahrig (1983) sites on the same dyke are averaged into dyke means. The pole is the "
     "Fisher mean of the site VGPs after unifying polarity (half of the sites are of each "
     "polarity), reported at the Pacific pole. Site VGPs are calculated from the site "
     "directions at the tabulated site coordinates (Fahrig et al., 1971; Park, 1974) or at "
     "the coordinates used by Denyszyn et al. (2009) (Fahrig & Schwarz, 1973; Christie & "
     "Fahrig, 1983; Palmer et al., 1983); for Palmer & Hayatsu (1975) and Park (1981), which "
     "give no site coordinates, the tabulated VGPs are used. The site of Fahrig et al. (1971) "
     "near Hardisty Lake (site 29) is placed at the location of its K-Ar sample, ca. 10 deg "
     "north of the latitude printed in their Table 2. With the printed location and the "
     "VGPs printed by the original authors, the pole is the published 6.7 N, 162.1 E. "
     "Palmer et al. (1983) "
     "directions are their anchored-line means. Most of the earlier studies used "
     "alternating-field blanket cleaning. No mean direction is given because the sites span "
     "ca. 2500 km. Denyszyn et al. (2009) regard the combined data of Park (1973) on "
     "Paleoproterozoic host rocks and Park (1974) on the Franklin dykes of the Daly Bay area "
     "(Park 1974 sites 31 and 33) as a positive, albeit weak, baked-contact test; the primary origin of "
     "the Franklin direction is further supported by the positive baked-contact test at "
     "site GR of Denyszyn et al. (2009) and by the dual polarity. The age is from U-Pb "
     "zircon dates on Franklin sills (Pu et al., 2022). This pole is used in the "
     "compilations of Evans et al. (2021) and Swanson-Hysell et al. (2026) in preference "
     "to the grand mean with Greenland and the Ellesmere microplate restored, because it "
     "rests only on sites in autochthonous Laurentia and not on those restorations."),
    with_dir=False, extra_codes=('DE-VGP', 'ST-C'), contact_test='C+', country='Canada'))
# D+09cg: Devon and Ellesmere (Ellesmere microplate) rotated by -20 deg about 72 N, 274 E
# (Denyszyn et al., 2009); Greenland by -13.8 deg about 67.5 N, 241.5 E (Roest &
# Srivastava, 1989, as in the compilation workbook); [Euler lat, Euler lon, angle] for
# pmag.pt_rot (positive angle counterclockwise)
EULER_MICROPLATE, EULER_GREENLAND = [72.0, 274.0, -20.0], [67.5, 241.5, -13.8]
restored = []
for sub, euler in [(canada, EULER_MICROPLATE), (green, EULER_GREENLAND)]:
    sub = sub.copy()
    lats, lons = pmag.pt_rot(euler, list(sub['vgp_lat']), list(sub['vgp_lon']))
    sub['vgp_lat'], sub['vgp_lon'] = lats, lons
    restored.append(sub)
allcg = pd.concat([gm] + restored)
loc_rows.append(loc_row(
    'Franklin event grand mean, Canada and Greenland restored',
    'Franklin event grand mean ca. 719 Ma pole (Canada and Greenland, restored)', allcg,
    ':'.join([CIT_D09] + sorted(set(gm['citations'])) + [CIT_PU22]),
    ("Grand-mean Franklin pole of Denyszyn et al. (2009, their pole D+09cg) combining the "
     "56 units of the previous Canadian studies with their 22 Devon, Ellesmere and Greenland "
     "sites after restoring Greenland and the proposed Ellesmere microplate (Ellesmere, "
     "Devon, Cornwallis and Somerset islands) to North America. The pole is the Fisher mean "
     "of the site VGPs after unifying polarity, reported at the Pacific pole, with the VGPs "
     "of the Devon and Ellesmere sites rotated by -20 deg about 72 N, 274 E (the "
     "microplate rotation of Denyszyn et al., 2009) and those of the Greenland sites by "
     "-13.8 deg about 67.5 N, 241.5 E (the Greenland reconstruction of Roest & Srivastava, "
     "1989, used by Denyszyn et al., 2009). The site VGPs in the sites table are in "
     "present-day coordinates. The site of Fahrig et al. (1971) near Hardisty Lake is "
     "placed at the location of its K-Ar sample (see the previous-studies grand mean), "
     "which does not change this pole at the precision published. The position of the "
     "pole depends on the restoration rotations: the microplate rotation of ca. 20 deg is "
     "inferred in part from these Devon and Ellesmere data (Denyszyn et al., 2009), and "
     "with a rotation of -35 deg about the same Euler pole the mean is 9.0 N, 161.7 E. "
     "Field tests as for the previous-studies grand mean, with the positive baked-contact "
     "test at site GR of Denyszyn et al. (2009) among the sites of this pole."),
    with_dir=False, extra_codes=('DE-VGP', 'ST-C'), contact_test='C+'))
locs = pd.DataFrame(loc_rows)
for c in ['dir_tilt_correction', 'dir_n_sites', 'pole_n_sites']:
    locs[c] = pd.to_numeric(locs[c]).astype('Int64')

# ---- ages table: measured dates ---------------------------------------------------
TONIAN = dict(timescale_eon='Proterozoic', timescale_era='Neoproterozoic',
              timescale_period='Tonian')
PU = [  # sample, area, 206Pb/238U weighted mean, X, Y, Z (2 sigma), n, note
    ('17RAT-R35B1', 'Coppermine area', 719.86, 0.21, 0.30, 0.79, '6 of 15',
     'older grains of the main cluster; the full data set spreads along concordia and no single '
     'weighted mean could be determined'),
    ('93JP-71JB', 'Minto Inlier, Victoria Island (sill)', 719.04, 0.19, 0.28, 0.79, '9 of 13', ''),
    ('93JP-93K', 'Duke of York Inlier, Victoria Island (sill)', 718.77, 0.30, 0.36, 0.82,
     '4 of 10', 'four oldest, most concordant analyses'),
    ('93JP-93L', 'Duke of York Inlier, Victoria Island (sill)', 718.96, 0.21, 0.29, 0.79,
     '11 of 13', ''),
    ('14RAT-513A', 'Brock Inlier (sill)', 718.61, 0.30, 0.36, 0.82, '9 of 10',
     'metamict zircon; nine oldest analyses'),
    ('F1966', 'Great Slave Lake (Douglas Peninsula sill area)', 719.08, 0.22, 0.30, 0.79,
     '6 of 18', 'six older overlapping analyses; the full data set spreads along concordia'),
]
age_rows = [{
    'location': GM, 'age': a, 'age_sigma': round(x / 2, 3), 'age_unit': 'Ma',
    'method_codes': 'GM-UPB', 'citations': CIT_PU22, **TONIAN,
    'description': (f'CA-ID-TIMS zircon weighted-mean 206Pb/238U date, Franklin LIP sample {s}, '
                    f'{area} (n = {n}); {a:.2f} +/- {x:.2f}/{y:.2f}/{z:.2f} Ma (2 sigma; internal / with tracer '
                    f'calibration / with decay-constant uncertainties; Pu et al., 2022)'
                    + (f'; {note}' if note else '') + '.')}
    for s, area, a, x, y, z, n, note in PU]
age_rows.append({
    'location': GM, 'age': 718.94, 'age_sigma': 0.8, 'age_unit': 'Ma',
    'method_codes': 'GM-UPB-CC-T0', 'citations': CIT_PU22, **TONIAN,
    'description': ('U-Pb baddeleyite upper-intercept date, Franklin sample FA700408, Baffin '
                    'Island, 718.94 +/- 1.60 Ma (95% confidence interval, n = 3 of 4, '
                    'constrained to a present-day Pb-loss line; Pu et al., 2022).')})
CIT_DU23 = '10.1016/j.epsl.2023.118259'
DUFOUR = [  # sample, location, area and coordinates, date, X, Y, Z (2 sigma), n, Bayesian age, note
    ('FH', GM, 'Fury and Hecla dyke, Baffin Island (70.0600 N, 85.8851 W)', 717.73, 0.72, 0.77,
     1.1, '3 of 7', '717.72 +0.60/-0.69', 'three oldest overlapping concordant grains'),
    ('DS', GM, 'Dybbol sill, Baffin Island (70.0611 N, 85.8838 W)', 718.43, 0.21, 0.36, 0.85,
     '6 of 11', '718.43 +/- 0.2', 'six oldest concordant grains'),
    ('TD', D09_LOC['Greenland'], 'Thule dyke, northwest Greenland (76.8388 N, 70.2579 W)',
     718.12, 0.94, 0.98, 1.2, '1 of 6', '718.08 +0.68/-0.89', 'oldest grain'),
    ('SLSC', D09_LOC['Greenland'], 'Steensby Land sill complex, northwest Greenland '
     '(76.5469 N, 68.5358 W)', 718.82, 0.26, 0.40, 0.87, '5 of 10', '718.82 +/- 0.25',
     'main cluster of concordant grains'),
]
for s, loc, area, a, x, y, z, n, bayes, note in DUFOUR:
    age_rows.append({
        'location': loc, 'age': a, 'age_sigma': round(x / 2, 3), 'age_unit': 'Ma',
        'method_codes': 'GM-UPB', 'citations': CIT_DU23, **TONIAN,
        'description': (f"CA-ID-TIMS zircon {'single-grain' if n.startswith('1 ') else 'weighted-mean'} "
                        f'206Pb/238U date, Franklin sample {s}, '
                        f'{area}; {a:.2f} +/- {x:.2f}/{y:.2f}/{z:.2f} Ma (2 sigma; internal / '
                        f'with tracer calibration / with decay-constant uncertainties; n = {n}, '
                        f'{note}; Dufour et al., 2023). The dykes cut the sills, and a Bayesian '
                        f'model with that order gives {bayes} Ma.')})
D09AGES = [  # site, age, 2 sigma, method, text
    ('CG', 721, 2, 'GM-UPB', 'average 206Pb/238U age of the three oldest baddeleyite fractions '
     '(1, 5 and 6); other fractions range down to 714 +/- 2 Ma'),
    ('QA', 721, 4, 'GM-UPB', 'average 206Pb/238U age of baddeleyite fractions'),
    ('BG', 726, 24, 'GM-UPB', 'average 207Pb/206Pb age of baddeleyite fractions'),
    ('GF', 712, 2, 'GM-UPB', 'average 206Pb/238U age of three baddeleyite fractions; the sill is '
     'paleomagnetically indistinguishable from the dykes but distinctly younger'),
]
for s, a, e2, mc, txt in D09AGES:
    loc = D09_LOC['Greenland'] if s in ('QA', 'GF') else D09_LOC['Devon']
    age_rows.append({
        'location': loc, 'site': s, 'age': a, 'age_sigma': e2 / 2, 'age_unit': 'Ma',
        'method_codes': mc, 'citations': CIT_D09, **TONIAN,
        'description': f'U-Pb baddeleyite, site {s}: {txt}; {a} +/- {e2} Ma (2 sigma; '
                       'Denyszyn et al., 2009).'})
ages = pd.DataFrame(age_rows)
# ages.sample would have to exist in a samples table, so sample IDs go in the description
ages = ages[['location', 'site', 'age', 'age_sigma', 'age_unit', 'method_codes',
             'citations', 'timescale_eon', 'timescale_era', 'timescale_period', 'description']]


def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False)


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
write_magic(ages, 'ages', os.path.join(OUT, 'ages.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Denyszyn2009_Franklin_LIP_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    f.write('tab delimited\tlocations\n'); locs.to_csv(f, sep='\t', index=False)
    f.write('>>>>>>>>>>\n')
    f.write('tab delimited\tsites\n'); sites.to_csv(f, sep='\t', index=False)
    f.write('>>>>>>>>>>\n')
    f.write('tab delimited\tages\n'); ages.to_csv(f, sep='\t', index=False)

for r in loc_rows:
    print(f"  {r['location'][:55]:55s} pole {r['pole_lat']}/{r['pole_lon']} A95 "
          f"{r['pole_alpha95']} K {r.get('pole_k')} N {r['pole_n_sites']} rev% "
          f"{r['pole_reversed_perc']}; dir {r.get('dir_dec')}/{r.get('dir_inc')}")
print(f'  sites: {len(sites)} rows; ages: {len(ages)} rows; wrote {combined}')

# poles as published by Denyszyn et al. (2009) and, for the grand mean, as used in the
# Swanson-Hysell et al. (2026) compilation (lat, lon, A95, n_sites)
EXPECT = {
    D09_LOC['Devon']: (5.8, 189.7, 10.7, 12),
    D09_LOC['Greenland']: (8.8, 178.5, 7.2, 10),
    'Franklin dykes, Arctic Canada and northwest Greenland': (7.2, 184.5, 6.7, 22),
    GM: (6.6, 162.1, 3.0, 56),       # published 6.7/162.1 (VGPs calculated; site 29 relocated)
    'Franklin event grand mean, Canada and Greenland restored': (8.4, 163.8, 2.8, 78),
}
for r in loc_rows:
    got = (r['pole_lat'], r['pole_lon'], r['pole_alpha95'], r['pole_n_sites'])
    assert got == EXPECT[r['location']], (r['location'], got)
assert len(gm) == 56 and int(gm['dir_n_samples'].sum()) == 324
# with the printed location and VGP of Fahrig et al. (1971) site 29, D+09can and D+09cg are
# the published poles
PRINTED_VGP = {r.site: observed_vgp(r.dec, r.inc, r.a95, r.lat, r.lon,
                                    (r.vgp_lat_paper, r.vgp_lon_paper))[:2]
               for r in PRIOR.itertuples()
               if r.study in CALC_VGP and not np.isnan(r.vgp_lat_paper)}
PRINTED_VGP['F71-29'] = (13.0, 153.0)     # printed VGP of the printed location (54.5 N)
for sub, published in [(gm, (6.7, 162.1, 3.0)), (allcg, (8.4, 163.8, 2.8))]:
    sub = sub.copy()
    for site, (la, lo) in PRINTED_VGP.items():
        sub.loc[sub['site'] == site, ['vgp_lat', 'vgp_lon']] = [la, lo]
    p, _, _ = poles(sub)
    assert (r1(p['inc']), r1(p['dec']), r1(p['alpha95'])) == published, p
print('-I- all five poles match the expected values')

from validate_magic_contribution import validate_upload_file
validate_upload_file(combined, tables=['locations', 'sites', 'ages'], server=True)
