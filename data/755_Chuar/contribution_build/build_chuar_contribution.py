"""Build the MagIC contribution for the Chuar Group (Grand Canyon, Arizona) of Eyster et al. (2020).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Eyster, A., Weiss, B. P., Karlstrom, K., & Macdonald, F. A. (2020). Paleomagnetism of the Chuar
Group and evaluation of the late Tonian Laurentian apparent polar wander path with implications
for the makeup and breakup of Rodinia. GSA Bulletin 132, 710-738. doi:10.1130/B32012.1
Supplemental material: GSA Data Repository item 2019239 (doi:10.1130/2019239); a text dump of it
is `Eyster2020_SI_DR2.txt` in this folder.

Weil, A. B., Geissman, J. W., & Van der Voo, R. (2004). Paleomagnetism of the Neoproterozoic
Chuar Group, Grand Canyon Supergroup, Arizona: implications for Laurentia's Neoproterozoic APWP
and Rodinia break-up. Precambrian Research 129, 71-92. doi:10.1016/j.precamres.2003.09.016

Dates: Rooney, A. D., et al. (2018). Coupled Re-Os and U-Pb geochronology of the Tonian Chuar
Group, Grand Canyon. GSA Bulletin 130, 1085-1098. doi:10.1130/B31768.1. Dehler, C., et al.
(2023). Precise U-Pb age models refine Neoproterozoic western Laurentian rift initiation,
correlation, and Earth system changes. Precambrian Research 396, 107156.
doi:10.1016/j.precamres.2023.107156 (Table 3, Bayesian age model). Uncertainties of dates are
2 sigma; age_sigma is 1 sigma.

Sources and checks
------------------
- Eyster et al. (2020) Table 1 (poles), supplement Table S1 (locality coordinates and bedding),
  Table S2 (sample line and great-circle fits; transcribed from the page images, including the
  exclusion marks) and Table S5 (the Weil et al., 2004 sites as used by Eyster et al., 2020).
  Weil et al. (2004) Table 1 (site means, bedding, N/N0), transcribed from the page images.
- The 17 Eyster samples of the HT Carbon Butte-Awatubi result are the linear HT fits with
  MAD < 15 deg that are not marked as excluded in Table S2. Samples A1303-3A and A1303-3B each
  have two HT fits; the lower-temperature fit is marked as excluded, so the higher-temperature
  fit is used. Each sample is entered as a site with dir_n_samples = 1, as the authors treat
  every sample as a site (each horizon may be a spot reading of the field). The 17 reproduce
  Table 1: in situ 108.0/-31.2 (k 16.3), tilt-corrected 98.1/-19.3 (k 45.8, alpha95 5.3),
  pole 12.5 N, 161.6 E (A95 4.0) (asserted below).
- The six Weil et al. (2004) sites of the Carbon Butte-Awatubi pole are JG-53-2 (Awatubi Member)
  and AW-13-4, AW-12-10, AW-12-11, AW-12-14A and AW-12-14B (Carbon Butte Member); AW-12-13 is
  excluded (alpha95 > 20; Table S5). These six give 16.2 N, 166.2 E (A95 6.8), as in the
  Table 1 footnote, and with the 17 samples the Combined HT Carbon Butte-Awatubi pole
  (asserted below).
- Weil et al. (2004) Table 1 prints tilt-corrected directions of 303/29 for Walcott site
  JG-60-18 and 290/30 for JG-60-6 (as does Table S5). Restoring the in-situ directions about
  each site's bedding gives 290.0/30.0 for JG-60-18 and 302.7/29.4 for JG-60-6, so each site
  carries the tilt-corrected direction that its bedding reproduces. Neither site is in a pole.
- Weil et al. (2004) k and alpha95 are stored as printed; for AW-12-11 (k 258, alpha95 16) and
  AW-12-13 (k 20, alpha95 21) they are not mutually consistent. For sample A1307-2B, restoring
  the in-situ direction with the locality bedding gives a direction 5.7 deg from the tabulated
  tilt-corrected one; both are stored as published.
- Coordinates: Eyster localities from Table S1 (degrees W converted to E). Weil et al. (2004)
  give no site coordinates; the Weil sites carry approximate coordinates digitized from their
  sampling map (Fig. 1). Eyster et al. (2020) calculated VGPs at one point (36.30 N, 248.10 E);
  the VGPs here are calculated at each site, which changes the poles by less than 0.1 deg.
- The MT result of Table 1 (N = 20) is not reproduced by the qualifying MT fits of Table S2
  (19 fits, 2.4 deg from the published mean) and is not entered.

Structure and interpretation
----------------------------
- Reference: Eyster et al. (2020). The Weil et al. (2004) site rows cite Weil et al. (2004).
- One location, 'Chuar Group', holds every site, including those in no pole (the HTs and
  present-day-field overprints, Weil site AW-12-13, and the Walcott and Jupiter sites). Each
  site carries the age of its member. The location has a row without a pole that lists every
  site (age bracket from the base of the Jupiter Member to the top of the Walcott Member,
  Dehler et al., 2023, model 95% ranges), and three result rows, each listing its sites:
    1. the combined Chuar Group pole used in the Swanson-Hysell et al. (2026) compilation: the
       Fisher mean of the tilt-corrected VGPs of the 37 sites of the two member poles of
       Eyster et al. (2020, Table 1), Carbon Butte-Awatubi (23) and Carbon Canyon (14). Its age,
       755 Ma (748-761 Ma), spans the median ages of the Dehler et al. (2023) age model for the
       base of the Awatubi Member (748.35 Ma) and the base of the Carbon Canyon Member
       (760.68 Ma);
    2. and 3. the two member poles of Table 1, which the authors interpret as primary and ca.
       10 Myr apart in age; both are kept because they record apparent polar wander. Their ages
       span the Dehler et al. (2023) model 95% ranges of the member boundaries (model ages, not
       dates).
- The steep HTs component of 10 Eyster samples and the nine Weil et al. (2004) sites that carry
  only a present-day-field overprint are described rather than given pole rows. Like every
  site, they carry the depositional age of their member; the authors' interpretation of the
  overprint ages (recent, Laramide or Neogene) is given in the descriptions. The two Walcott and two Jupiter sites of Weil et al. (2004) are in no pole of Eyster
  et al. (2020).
- Poles: Fisher mean of the polarity-unified tilt-corrected site VGPs, in the ca. 163 E
  hemisphere of Eyster et al. (2020). Location mean directions: Fisher mean in the majority
  (observed) polarity. Directions and poles are not corrected for inclination shallowing; the
  authors' E-I flattening factors and corrected poles are given in the descriptions.
- Field tests: fold test F+ for both member poles (Tauxe and Watson, 1994, bootstrap; Eyster et
  al., 2020) and for the combined pole, all of whose sites are in those two tests. Reversal test
  Rb for the Carbon Butte-Awatubi pole: the McFadden and McElhinny (1990) test on its 23
  tilt-corrected directions is class B (means 7.5 deg apart, critical angle 7.7 deg); Eyster et
  al. (2020) report class C for their 17 samples and class B for the 23 after their f = 0.9
  flattening correction. No reversal test for the combined pole, whose polarity is largely
  confounded with member.
- Sample lithologies: Eyster et al. (2020) describe the Carbon Butte samples as spanning the
  basal sandstone and the white sandstone (entered as Sandstone:Siltstone, as the interval
  includes siltstone interbeds) and do not state the lithology of the Awatubi samples (Not
  Specified).
- The ages table holds the dates of Rooney et al. (2018) and the maximum depositional age of
  the Chuar Group of Dehler et al. (2023).
"""
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(OUT))

EYS = '10.1130/B32012.1'
WEIL = '10.1016/j.precamres.2003.09.016'
ROONEY = '10.1130/B31768.1'
DEHLER = '10.1016/j.precamres.2023.107156'

ALL = 'Chuar Group'
CBA = 'Chuar Group, Carbon Butte-Awatubi members'
CCM = 'Chuar Group, Carbon Canyon Member'
WAL = 'Chuar Group, Walcott Member'
JUP = 'Chuar Group, Jupiter Member'
CHRM = 'characteristic'
HTS = 'HTs overprint'
PDF = 'present-day field overprint'

# location: (age, age_low, age_high); Dehler et al. (2023) Table 3 medians and 2.5-97.5% quantiles
AGES = {ALL: (755, 748, 761),
        CBA: (749.0, 742.1, 754.9),     # base Carbon Butte 749.74 (743.90-754.87) to base Awatubi 748.35 (742.14-753.47)
        CCM: (756.1, 746.0, 766.0),     # base Carbon Canyon 760.68 (755.45-765.97) to base Duppa 751.44 (745.99-756.32)
        WAL: (734.4, 727.9, 746.1),     # base Walcott 739.92 (733.12-746.12) to top Walcott 728.93 (727.93-729.81)
        JUP: (762.7, 755.5, 768.4)}     # base Jupiter 764.70 (758.10-768.40) to base Carbon Canyon 760.68 (755.45-765.97)

# ---------------------------------------------------------------------------
# Eyster et al. (2020) Table S1: locality -> lat N, lon W, bedding strike, dip, member
EY_LOC = {
    'A1303': (36.268, 111.886, 133, 39, 'Carbon Butte'),
    'A1304': (36.268, 111.886, 140, 39, 'Awatubi'),
    'A1305': (36.274, 111.886, 157, 33, 'Carbon Butte'),
    'A1306': (36.276, 111.890, 70, 22, 'Carbon Butte'),
    'A1307': (36.274, 111.890, 65, 24, 'Awatubi'),
    'A1309': (36.274, 111.891, 44, 29, 'Carbon Butte'),
    'A1310': (36.167, 111.839, 304, 12, 'Carbon Butte'),
}
# Table S2: sample, geographic D, I, tilt-corrected D, I, unblocking range, MAD
EY_HT = [
    ('A1303-2A', 120.5, -29.5, 105.1, -14.7, '590-650 C', 5.7),
    ('A1303-2B', 109.9, -35.6, 94.0, -14.0, '530-680 C', 9.1),
    ('A1303-2C', 108.3, -36.9, 92.1, -14.4, '570-635 C', 8.8),
    ('A1303-3A', 104.3, -25.6, 95.8, -3.1, '575-635 C', 10.4),
    ('A1303-3B', 95.9, -31.8, 86.1, -4.5, '605-680 C', 12.9),
    ('A1304-1B', 123.6, -40.6, 101.8, -22.0, '400-570 C', 14.7),
    ('A1304-2A', 129.6, -50.3, 97.7, -32.0, '500-605 C', 14.2),
    ('A1304-2B', 122.2, -42.0, 100.0, -22.4, '500-590 C', 6.3),
    ('A1305-1A', 287.6, 51.2, 273.3, 23.2, '450-690 C', 14.6),
    ('A1305-4A', 298.7, 53.4, 279.1, 28.2, '670-690 C', 13.6),
    ('A1305-4B', 275.9, 54.8, 264.8, 24.2, '580-680 C', 5.2),
    ('A1306-1A', 104.9, 10.2, 106.1, -2.8, '450-660 C', 13.2),
    ('A1307-2B', 111.2, -11.9, 101.5, -33.4, '375-450 C', 5.7),
    ('A1310-2A', 102.4, -20.3, 107.0, -24.3, '300-665 C', 7.8),
    ('A1310-2B', 88.8, -8.8, 90.9, -15.6, '300-665 C', 8.3),
    ('A1310-3A', 98.6, -27.6, 105.0, -32.2, '550-660 C', 14.4),
    ('A1310-3B', 106.2, -10.7, 108.7, -14.1, '565-660 C', 5.7),
]
EY_HTS = [
    ('A1305-3B', 358.3, 56.9, 307.2, 54.2, '650-695 C', 6.2),
    ('A1306-1B', 333.1, 66.0, 284.3, 86.6, '450-660 C', 5.3),
    ('A1306-2A', 358.3, 53.8, 20.6, 73.4, '560-700 C', 11.9),
    ('A1306-2B', 351.2, 45.4, 0.3, 66.7, '660-700 C', 7.5),
    ('A1309-1A', 30.3, 58.7, 75.9, 54.3, 'NRM-605 C', 14.7),
    ('A1309-1B', 340.4, 19.0, 349.4, 44.4, '475-680 C', 14.1),
    ('A1309-2A', 10.4, 55.3, 59.1, 61.3, '100-680 C', 5.8),
    ('A1309-2B', 30.2, 42.2, 57.7, 42.9, '400-660 C', 6.9),
    ('A1309-3A', 341.2, 33.7, 358.6, 58.0, '620-660 C', 5.8),
    ('A1309-3B', 341.9, 45.4, 12.3, 68.0, '200-675 C', 5.2),
]

# ---------------------------------------------------------------------------
# Weil et al. (2004) Table 1: site, member, strike, dip, n, N0, in-situ D, I, alpha95, k,
# tilt-corrected D, I (None where only the in-situ direction is given), lat, lon (digitized
# from Fig. 1). JG-60-18 / JG-60-6: tilt-corrected directions as reproduced by their bedding.
W04 = [
    ('JG-60-18', 'Walcott', 200, 24, 7, 7, 290, 54, 18, 27, 290, 30, 36.203796, 248.149985),
    ('JG-60-6', 'Walcott', 223, 21, 7, 7, 299, 50, 13, 28, 303, 29, 36.205515, 248.152709),
    ('JG-53-2', 'Awatubi', 139, 30, 6, 7, 107, -31, 11, 76, 98, -13, 36.258118, 248.110988),
    ('T-52-1', 'Carbon Butte', 148, 22, 8, 8, 5, 64, 8, 56, None, None, 36.274493, 248.118863),
    ('T-65-1', 'Carbon Butte', 340, 18, 5, 8, 4, 54, 11, 50, None, None, 36.17234, 248.157581),
    ('T-65-2', 'Carbon Butte', 340, 18, 6, 8, 5, 53, 7, 107, None, None, 36.17559, 248.154966),
    ('JG-60-2', 'Carbon Butte', 195, 30, 6, 8, 0, 61, 6, 122, None, None, 36.207338, 248.170274),
    ('JG-60-3', 'Carbon Butte', 185, 46, 6, 8, 5, 55, 8, 151, None, None, 36.204852, 248.169978),
    ('JG-53-4', 'Carbon Butte', 155, 29, 6, 8, 13, 57, 10, 51, None, None, 36.274011, 248.121659),
    ('JG-53-5', 'Carbon Butte', 155, 29, 9, 9, 345, 63, 10, 30, None, None, 36.272798, 248.125496),
    ('AW-13-4', 'Carbon Butte', 5, 20, 8, 8, 272, 3, 9, 39, 272, 23, 36.26907, 248.103964),
    ('AW-13-5', 'Carbon Butte', 45, 22, 8, 8, 6, 51, 7.8, 61, None, None, 36.267261, 248.102825),
    ('AW-12-10', 'Carbon Butte', 313, 14, 6, 8, 88, -25, 12, 65, 94, -34, 36.247153, 248.133426),
    ('AW-12-11', 'Carbon Butte', 332, 21, 6, 8, 103, -13, 16, 258, 108, -28, 36.2495, 248.133917),
    ('AW-12-12', 'Carbon Butte', 334, 22, 8, 8, 11, 68, 9, 61, None, None, 36.24478, 248.12932),
    ('AW-12-13', 'Carbon Butte', 16, 9, 7, 8, 85, -28, 21, 20, 83, -36, 36.242285, 248.12589),
    ('AW-12-14A', 'Carbon Butte', 341, 21, 6, 8, 92, -20, 12, 60, 97, -39, 36.242383, 248.121909),
    ('AW-12-14B', 'Carbon Butte', 341, 21, 8, 8, 97, -19, 12, 39, 103, -37, 36.242021, 248.118195),
    ('AW-12-8', 'Carbon Canyon', 162, 43, 6, 8, 260, 60, 13, 27, 256, 17, 36.255328, 248.141498),
    ('AW-12-9', 'Carbon Canyon', 160, 60, 7, 8, 256, 48, 9, 90, 254, -11, 36.253493, 248.142241),
    ('AW-12-15', 'Carbon Canyon', 303, 12, 7, 8, 277, 31, 15, 24, 285, 35, 36.236527, 248.119935),
    ('AW-12-19', 'Carbon Canyon', 310, 39, 7, 7, 299, -8, 16, 16, 296, 1, 36.163766, 248.147031),
    ('AW-12-20', 'Carbon Canyon', 290, 16, 6, 7, 284, 5, 7, 53, 285, 7, 36.16525, 248.146738),
    ('AW-13-1', 'Carbon Canyon', 110, 25, 7, 7, 264, 6, 11, 31, 264, -5, 36.279811, 248.113842),
    ('AW-13-2', 'Carbon Canyon', 105, 20, 7, 7, 264, 18, 15, 16, 259, 10, 36.277885, 248.113658),
    ('AW-13-3', 'Carbon Canyon', 85, 20, 7, 8, 254, 26, 14, 21, 245, 22, 36.278323, 248.11634),
    ('AW-13-6', 'Carbon Canyon', 35, 30, 7, 8, 241, -20, 12, 31, 247, -5, 36.273822, 248.106226),
    ('AW-13-7', 'Carbon Canyon', 37, 29, 6, 7, 260, -19, 17, 22, 263, 2, 36.273822, 248.106226),
    ('AW-13-8', 'Carbon Canyon', 25, 13, 7, 7, 267, -2, 14, 24, 267, 10, 36.273822, 248.106226),
    ('AW-13-9', 'Carbon Canyon', 30, 30, 6, 7, 252, -29, 15, 27, 259, -7, 36.273822, 248.106226),
    ('AW-13-10', 'Carbon Canyon', 25, 23, 7, 7, 264, -6, 14, 25, 263, 13, 36.273822, 248.106226),
    ('AW-13-11', 'Carbon Canyon', 30, 30, 7, 8, 263, -30, 19, 13, 269, -5, 36.273822, 248.106226),
    ('AW-12-22', 'Jupiter', 287, 17, 8, 8, 268, -11, 16, 20, 266, -5, 36.149317, 248.156054),
    ('AW-12-23', 'Jupiter', 75, 1, 8, 8, 270, 3, 6, 112, 272, 5, 36.146884, 248.155693),
]
WEIL_CBA_USED = ['JG-53-2', 'AW-13-4', 'AW-12-10', 'AW-12-11', 'AW-12-14A', 'AW-12-14B']
WEIL_PDF = ['T-52-1', 'T-65-1', 'T-65-2', 'JG-60-2', 'JG-60-3', 'JG-53-4', 'JG-53-5',
            'AW-13-5', 'AW-12-12']
MEMBER_LOC = {'Walcott': WAL, 'Awatubi': CBA, 'Carbon Butte': CBA,
              'Carbon Canyon': CCM, 'Jupiter': JUP}
MEMBER_FM = {'Walcott': 'Kwagunt', 'Awatubi': 'Kwagunt', 'Carbon Butte': 'Kwagunt',
             'Carbon Canyon': 'Galeros', 'Jupiter': 'Galeros'}

MC_W = 'SO-MAG:LP-DIR-T:DE-BFL:DE-FM'
MC_E = 'SO-SM:LT-LT-Z:LP-DIR-AFT:DE-BFL'
MC_LOC_E = 'SO-SM:LT-LT-Z:LP-DIR-AFT:DE-BFL:LP-IRM:LP-IRM-DCD'
MC_LOC_W = 'SO-MAG:LP-DIR-T:LP-DIR-CHEM:LP-IRM-3D:DE-BFL:DE-FM'
LITH_W = 'Sandstone:Siltstone:Mudstone'
GEOG = dict(continent_ocean='North America', country='United States of America',
            state_province='Arizona', terranes='Laurentia')


def r1(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 1)


def vgp(dec, inc, a95, lat, lon):
    plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95 if a95 is not None else np.nan, lat, lon)
    return r1(plat), r1(plon), r1(dp), r1(dm)


def codes(*parts):
    """Join colon-delimited method-code strings without duplicates."""
    out = []
    for p in parts:
        for c in p.split(':'):
            if c and c not in out:
                out.append(c)
    return ':'.join(out)


rows = []


def add(**kw):
    group = kw['location']                 # member group, for the site age and pole membership
    age, lo, hi = AGES[group]
    plat, plon, dp, dm = vgp(kw['dir_dec'], kw['dir_inc'], kw.get('dir_alpha95'),
                             kw['lat'], kw['lon'])
    rows.append({
        'site': kw['site'], 'location': ALL, '_group': group, 'result_type': 'i',
        'result_quality': kw.get('quality', 'g'), 'method_codes': kw['mc'],
        'citations': kw['cit'], 'geologic_classes': 'Sedimentary',
        'geologic_types': 'Sediment Layer', 'lithologies': kw['lith'],
        'bed_dip': kw['bed_dip'], 'bed_dip_direction': kw['bed_dip_direction'],
        'lat': kw['lat'], 'lon': kw['lon'],
        'age': age, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
        'dir_comp_name': kw['dir_comp_name'],
        'dir_tilt_correction': kw['tc'], 'dir_dec': kw['dir_dec'], 'dir_inc': kw['dir_inc'],
        'dir_alpha95': kw.get('dir_alpha95'), 'dir_k': kw.get('dir_k'),
        'dir_n_samples': kw['n'], 'dir_n_total_samples': kw.get('n0'),
        'dir_n_specimens': kw.get('nspec'),
        'vgp_lat': plat, 'vgp_lon': plon, 'vgp_dp': dp, 'vgp_dm': dm,
        'description': kw['desc']})


# ---- Eyster et al. (2020) samples, each entered as a site ---------------------------
def eyster_rows(data, comp):
    for (s, gd, gi, td, ti, rng, mad) in data:
        la, lonw, strike, dip, member = EY_LOC[s[:5]]
        lon = round(360 - lonw, 3)
        # Eyster et al. (2020) give the Carbon Butte samples as spanning the basal and white
        # sandstones (with siltstone interbeds); the Awatubi lithology is not stated
        lith = 'Sandstone:Siltstone' if member == 'Carbon Butte' else 'Not Specified'
        what = ('shallow high-temperature (HT) component' if comp == CHRM else
                'steep high-temperature (HTs) component, interpreted by the authors as a '
                'chemical remanence overprint')
        area = 'Carbon Butte, Lava Chuar area' if s[:5] == 'A1310' else 'Nankoweap Butte'
        desc = (f"Individually oriented core from a block sample at locality {s[:5]}, "
                f"{member} Member, Kwagunt Formation, Chuar Group ({area}); treated as a site "
                f"by the authors. Principal-component line fit of the {what} (unblocking range "
                f"{rng}, MAD {mad}; Table S2). Coordinates are the locality values of Table S1.")
        if s in ('A1303-3A', 'A1303-3B') and comp == CHRM:
            desc += (' The sample has two HT fits in Table S2; the lower-temperature fit is '
                     'marked as excluded and the higher-temperature fit is entered.')
        if s == 'A1307-2B' and comp == CHRM:
            desc += (' Restoring the in-situ direction with the locality bedding gives a '
                     'direction ca. 6 degrees from the tabulated tilt-corrected direction; both '
                     'are stored as published.')
        for tc, d, i in ((0, gd, gi), (100, td, ti)):
            add(site=s, location=CBA, mc=MC_E, cit=EYS, lith=lith, bed_dip=dip,
                bed_dip_direction=(strike + 90) % 360, lat=la, lon=lon, dir_comp_name=comp,
                tc=tc, dir_dec=d, dir_inc=i, n=1, nspec=1, desc=desc)


eyster_rows(EY_HT, CHRM)
eyster_rows(EY_HTS, HTS)

# ---- Weil et al. (2004) sites ------------------------------------------------------
for (s, mem, strike, dip, n, n0, gd, gi, a95, k, td, ti, la, lo) in W04:
    desc = (f"{mem} Member, {MEMBER_FM[mem]} Formation, Chuar Group, eastern Grand Canyon; "
            f"site of multiple beds of hematite-cemented sandstone to mudstone (hand samples; "
            f"Weil et al., 2004, Table 1). Coordinates are approximate, digitized from the "
            f"sampling map (Weil et al., 2004, Fig. 1).")
    quality, comp = 'g', CHRM
    if s in WEIL_PDF:
        comp = PDF
        desc += (' All samples carry a stable north-directed, present-day-field overprint '
                 'unblocking up to 680 C and no recognizable ancient magnetization; only the '
                 'in-situ direction is tabulated.')
    if s == 'AW-12-13':
        quality = 'b'
        desc += (' Not used in the Carbon Butte-Awatubi pole of Eyster et al. (2020), which '
                 'excludes sites with alpha95 > 20 deg.')
    if mem == 'Walcott':
        desc += (' Not in the Carbon Butte-Awatubi or combined poles: Eyster et al. (2020) '
                 'suggest the Walcott directions may be biased by the steep overprint. Weil et '
                 'al. (2004, Table 1) print tilt-corrected directions of 303/29 for JG-60-18 '
                 'and 290/30 for JG-60-6; restoring each in-situ direction about its bedding '
                 'gives 290.0/30.0 for JG-60-18 and 302.7/29.4 for JG-60-6, and the '
                 'tilt-corrected row carries the direction that the bedding of this site '
                 'reproduces.')
    if mem == 'Jupiter':
        desc += (' Included in the Galeros Formation mean of Weil et al. (2004) but not in '
                 'the Carbon Canyon Member pole of Eyster et al. (2020).')
    if s == 'AW-13-4':
        desc += ' West-directed, downward (minority) polarity.'
    common = dict(site=s, location=MEMBER_LOC[mem], mc=MC_W, cit=WEIL, lith=LITH_W,
                  bed_dip=dip, bed_dip_direction=(strike + 90) % 360, lat=la, lon=lo,
                  dir_comp_name=comp, dir_alpha95=a95, dir_k=k, n=n, n0=n0, desc=desc,
                  quality=quality)
    add(tc=0, dir_dec=gd, dir_inc=gi, **common)
    if td is not None:
        add(tc=100, dir_dec=td, dir_inc=ti, **common)

sites = pd.DataFrame(rows)

# Walcott: the stored tilt-corrected directions follow from the in-situ direction and bedding
for s, want in (('JG-60-18', (290, 30)), ('JG-60-6', (303, 29))):
    r = sites[(sites.site == s) & (sites.dir_tilt_correction == 0)].iloc[0]
    d, i = pmag.dotilt(r.dir_dec, r.dir_inc, r.bed_dip_direction, r.bed_dip)
    assert pmag.angle([d, i], list(want))[0] < 1.5, s

# ---------------------------------------------------------------------------
tcg = sites[(sites.dir_tilt_correction == 100) & (sites.dir_comp_name == CHRM) &
            (sites.result_quality == 'g')]
cba = tcg[tcg._group == CBA]
ccm = tcg[tcg._group == CCM]
assert len(cba) == 23 and len(ccm) == 14
assert set(cba[cba.citations == WEIL].site) == set(WEIL_CBA_USED)
comb = pd.concat([cba, ccm])


def unify(pairs, ref):
    rc = pmag.dir2cart([ref[0], ref[1], 1.])
    return [p if np.dot(pmag.dir2cart([p[0], p[1], 1.]), rc) > 0 else
            [(p[0] + 180) % 360, -p[1]] for p in pairs]


def pole(sub, ref=(163, 5)):
    """Fisher mean of polarity-unified VGPs, in the hemisphere of ref (lon, lat)."""
    return pmag.fisher_mean(unify(sub[['vgp_lon', 'vgp_lat']].values.tolist(), ref))


def majority_dir(sub):
    """Fisher mean of site directions in the majority (observed) polarity, and the
    percentage of sites of the minority polarity."""
    pairs = sub[['dir_dec', 'dir_inc']].values.tolist()
    m0 = pmag.doprinc(pairs)
    rc = pmag.dir2cart([m0['dec'], m0['inc'], 1.])
    same = np.array([np.dot(pmag.dir2cart([d, i, 1.]), rc) > 0 for d, i in pairs])
    if same.sum() < len(same) / 2:
        same = ~same
        m0 = {'dec': (m0['dec'] + 180) % 360, 'inc': -m0['inc']}
    m = pmag.fisher_mean(unify(pairs, (m0['dec'], m0['inc'])))
    return m, round(100 * (~same).sum() / len(same), 1)


# Eyster et al. (2020) Table 1 checks: 17 samples alone, and the 6 Weil sites alone
ey17 = cba[cba.citations == EYS]
d17, _ = majority_dir(ey17)
p17 = pole(ey17)
assert (round(d17['dec'], 1), round(d17['inc'], 1), round(d17['k'], 1),
        round(d17['alpha95'], 1)) == (98.1, -19.3, 45.8, 5.3), d17
assert (round(p17['inc'], 1), round(p17['dec'], 1), round(p17['alpha95'], 1)) == \
    (12.5, 161.6, 4.0), p17
pw6 = pole(cba[cba.citations == WEIL])   # VGPs at the site coordinates: within 0.1 deg
assert pmag.angle([pw6['dec'], pw6['inc']], [166.2, 16.2])[0] < 0.1 and \
    round(pw6['alpha95'], 1) == 6.8, pw6


def locrow(group, rname, sub, p, d, rev, mc, cit, desc, liths=LITH_W, tests=None,
           ext=None):
    a, lo, hi = AGES[group]
    r = {'location': ALL, 'location_type': 'Region', 'result_name': rname,
         'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(sub.site),
         'method_codes': mc, 'citations': cit, 'geologic_classes': 'Sedimentary',
         'lithologies': liths,
         'lat_s': round(sub.lat.min(), 3), 'lat_n': round(sub.lat.max(), 3),
         'lon_w': round(sub.lon.min(), 3), 'lon_e': round(sub.lon.max(), 3),
         'age': a, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
         'dir_tilt_correction': 100, 'dir_dec': r1(d['dec']), 'dir_inc': r1(d['inc']),
         'dir_alpha95': r1(d['alpha95']), 'dir_k': r1(d['k']), 'dir_n_sites': int(d['n']),
         'pole_comp_name': CHRM, 'pole_lat': r1(p['inc']), 'pole_lon': r1(p['dec']),
         'pole_alpha95': r1(p['alpha95']), 'pole_k': r1(p['k']),
         'pole_n_sites': int(p['n']), 'pole_reversed_perc': rev,
         'external_database_ids': ext, 'description': desc}
    if tests:
        r.update(tests)
    r.update(GEOG)
    return r


p_cba, (d_cba, rev_cba) = pole(cba), majority_dir(cba)
p_cc, (d_cc, rev_cc) = pole(ccm), majority_dir(ccm)
p_all, (d_all, rev_all) = pole(comb), majority_dir(comb)

loc_rows = [locrow(
    ALL, 'Chuar Group combined pole', comb, p_all, d_all, rev_all,
    codes(MC_LOC_E, MC_LOC_W, 'DE-VGP:ST-F'), f'{EYS}:{WEIL}:{DEHLER}',
    liths='Sandstone:Siltstone:Mudstone:Not Specified', tests={'fold_test': 'F+'},
    desc=(
        "Combined paleomagnetic pole for the Tonian Chuar Group, eastern Grand Canyon, "
        "Arizona, from the hematite-carried, dual-polarity characteristic remanence of the "
        "Carbon Canyon Member (Galeros Formation) and the Carbon Butte and basal Awatubi "
        "members (Kwagunt Formation). The pole is the Fisher mean of the tilt-corrected "
        "virtual geomagnetic poles of the sites of the two member poles of Eyster et al. "
        "(2020, Table 1): their 17 samples, each treated as a site, and 6 sites of Weil et "
        "al. (2004) for the Carbon Butte-Awatubi members, and 14 Carbon Canyon sites of Weil "
        "et al. (2004). The members differ in age by ca. 10 Myr and their poles are distinct, "
        "so the combined pole averages apparent polar wander between them, and polarity is "
        "largely confounded with member. Each member pole passes a fold test (Eyster et al., "
        "2020). Directions and poles are not corrected for inclination shallowing. The age "
        "spans the median ages of the Dehler et al. (2023) age model for the bases of the "
        "Awatubi Member (748.35 Ma) and the Carbon Canyon Member (760.68 Ma)."))]

loc_rows.append(locrow(
    CBA, 'Combined HT Carbon Butte-Awatubi pole', cba, p_cba, d_cba, rev_cba,
    codes(MC_LOC_E, MC_LOC_W, 'DE-VGP:ST-F:ST-R'), f'{EYS}:{WEIL}:{ROONEY}:{DEHLER}',
    liths='Sandstone:Siltstone:Mudstone:Not Specified',
    tests={'fold_test': 'F+', 'reversal_test': 'Rb'}, ext='GPMDB[9981]',
    desc=(
        "Paleomagnetic pole for the Carbon Butte and basal Awatubi members of the Kwagunt "
        "Formation, Chuar Group, Nankoweap Butte and Lava Chuar areas, eastern Grand Canyon "
        "(Combined HT Carbon Butte-Awatubi pole of Eyster et al., 2020, Table 1): 17 "
        "individually oriented samples of Eyster et al. (2020), each treated as a site because "
        "each horizon may be a spot reading of the field, and 6 sites of Weil et al. (2004); "
        "Weil site AW-12-13 (alpha95 > 20) is excluded. The shallow high-temperature remanence "
        "is carried by hematite and interpreted by the authors as a detrital remanence; it "
        "passes the Tauxe and Watson (1994) bootstrap fold test (60-101% unfolding) and the "
        "Watson V reversal test, although the authors note that the small number of reversed "
        "directions may bias the tests. The authors estimate an E-I flattening factor f = 0.9 "
        "and report a corrected pole of 14.2 N, 163.8 E; the pole here is uncorrected. The "
        "sites table also holds the steep, north-directed HTs component of 10 Carbon Butte "
        "samples, which fails the fold test (-10 to 47% unfolding) and is interpreted by the "
        "authors as a chemical remanence overprint, possibly related to the Laramide orogeny "
        "or Neogene volcanism (Eyster et al., 2020, Table 1: in-situ pole 84.7 N, 96.1 E, A95 "
        "12.6), and nine Weil et al. (2004) sites that carry only a present-day-field "
        "overprint. "
        "The sampled horizons lie just below the Re-Os date on basal Awatubi Member marcasite "
        "nodules (Rooney et al., 2018); the age bracket spans the Dehler et al. (2023) "
        "age-model 95% ranges for the bases of the Carbon Butte and Awatubi members.")))

loc_rows.append(locrow(
    CCM, 'HT Carbon Canyon pole', ccm, p_cc, d_cc, rev_cc,
    codes(MC_LOC_W, 'DE-VGP:ST-F'), f'{WEIL}:{EYS}:{ROONEY}:{DEHLER}',
    tests={'fold_test': 'F+'},
    desc=(
        "Paleomagnetic pole for the Carbon Canyon Member of the Galeros Formation, Chuar "
        "Group, eastern Grand Canyon (HT Carbon Canyon pole of Eyster et al., 2020, Table 1), "
        "from 14 sites of Weil et al. (2004); Eyster et al. (2020) exclude the two Jupiter "
        "Member sites that Weil et al. (2004) included in their Galeros Formation mean. The "
        "hematite-carried characteristic remanence was isolated at unblocking temperatures up "
        "to ca. 680 C after removal of a present-day-field overprint and passes the fold test "
        "(Tauxe and Watson, 1994, bootstrap, 72-109% unfolding, Eyster et al., 2020; maximum "
        "clustering at ca. 100% unfolding, Weil et al., 2004). Eyster et al. (2020) estimate "
        "an E-I flattening factor f = 0.5 and report a corrected pole of 0.5 S, 166.0 E; the "
        "pole here is uncorrected. Weil et al. (2004) note that some individual samples carry "
        "the opposite (east, shallow) polarity. The member contains the Re-Os date of Rooney "
        "et al. (2018); the age bracket spans the Dehler et al. (2023) age-model 95% ranges "
        "for the base and top of the member.")))

# the location as a whole: every site, no pole
all_sites = list(dict.fromkeys(sites.site))
loc_all = {
    'location': ALL, 'location_type': 'Region', 'sites': ':'.join(all_sites),
    'method_codes': codes(MC_LOC_E, MC_LOC_W), 'citations': f'{EYS}:{WEIL}:{ROONEY}:{DEHLER}',
    'geologic_classes': 'Sedimentary', 'lithologies': 'Sandstone:Siltstone:Mudstone:Not Specified',
    'lat_s': round(sites.lat.min(), 3), 'lat_n': round(sites.lat.max(), 3),
    'lon_w': round(sites.lon.min(), 3), 'lon_e': round(sites.lon.max(), 3),
    'age_low': AGES[WAL][1], 'age_high': AGES[JUP][2], 'age_unit': 'Ma',
    'description': (
        "Paleomagnetic sites of the Tonian Chuar Group, eastern Grand Canyon, Arizona: the "
        "samples of Eyster et al. (2020) from the Carbon Butte and basal Awatubi members "
        "(Kwagunt Formation), each treated as a site, and the sites of Weil et al. (2004) in "
        "the Jupiter and Carbon Canyon members (Galeros Formation) and the Carbon Butte, "
        "Awatubi and Walcott members (Kwagunt Formation). Each site carries the age of its "
        "member; the age bracket spans the Dehler et al. (2023) age-model 95% ranges from the "
        "base of the Jupiter Member to the top of the Walcott Member (model ages, not dates). "
        "The pole rows of this location list their sites. Sites in no pole: the steep HTs "
        "overprint of 10 Eyster samples and the nine Weil sites with only a present-day-field "
        "overprint (described with the Carbon Butte-Awatubi pole), Weil site AW-12-13 (alpha95 "
        "> 20 deg), the two Walcott Member sites, which lie on one limb of the Chuar syncline "
        "and whose directions Eyster et al. (2020) suggest may be biased by the steep "
        "overprint, and the two Jupiter Member sites, which Weil et al. (2004) included in "
        "their Galeros Formation mean but Eyster et al. (2020) exclude from the Carbon Canyon "
        "Member pole."),
    **GEOG}
loc_rows.insert(0, loc_all)
locs = pd.DataFrame(loc_rows)

# ---------------------------------------------------------------------------
AGE_ROWS = [
    (CCM, 757.0, 3.4, 'GM-REOS-I', ROONEY,
     'Re-Os isochron (model 1) age of organic-rich carbonates of the Carbon Canyon Member, '
     'Galeros Formation (sample A1407, Carbon Canyon), 757.0 +/- 6.8 Ma (2 sigma, n = 8, '
     'MSWD 0.47; Rooney et al., 2018).'),
    (CBA, 751.0, 3.8, 'GM-REOS-I', ROONEY,
     'Re-Os isochron (model 1) age of marcasite nodules in the basal Awatubi Member, '
     'Kwagunt Formation, Nankoweap Butte, 751.0 +/- 7.6 Ma (2 sigma, n = 5, MSWD 0.37; '
     'Rooney et al., 2018); interpreted to date nodule growth during deposition.'),
    (WAL, 729.0, 0.45, 'GM-UPB', ROONEY,
     'U-Pb CA-ID-TIMS weighted mean 206Pb/238U zircon age of a tuff in the uppermost '
     'Walcott Member, Kwagunt Formation, 729.0 +/- 0.9 Ma (2 sigma, MSWD 0.86; Rooney et '
     'al., 2018).'),
    (ALL, 770.1, 0.25, 'GM-UPB', DEHLER,
     'U-Pb CA-ID-TIMS age of the youngest detrital zircon in the basal Nankoweap '
     'Formation (sample CDGC5), 770.1 +/- 0.5 Ma (2 sigma); a maximum depositional age '
     'for the Chuar Group (Dehler et al., 2023).'),
]
ages = pd.DataFrame([{
    'location': ALL, 'age': a, 'age_sigma': s, 'age_unit': 'Ma', 'method_codes': mc,
    'citations': cit, 'timescale_eon': 'Proterozoic', 'timescale_era': 'Neoproterozoic',
    'timescale_period': 'Tonian', 'description': desc}
    for loc, a, s, mc, cit, desc in AGE_ROWS])


# ---------------------------------------------------------------------------
def fmt(df):
    df = df.dropna(axis=1, how='all').copy()
    for c in df.columns:
        if df[c].dtype.kind in 'fiu' or str(df[c].dtype) == 'Int64':
            df[c] = df[c].map(lambda x: '' if pd.isna(x) else f'{float(x):.10g}')
    return df.fillna('')


def table_text(df, kind):
    return f'tab delimited\t{kind}\n' + fmt(df).to_csv(sep='\t', index=False)


for c in ['dir_alpha95', 'dir_k', 'age', 'dir_n_samples', 'dir_n_total_samples',
          'dir_n_specimens']:
    sites[c] = pd.to_numeric(sites[c])

TABLES = (('locations', locs), ('sites', sites.drop(columns=['_group'])), ('ages', ages))
for kind, df in TABLES:
    with open(os.path.join(OUT, f'{kind}.txt'), 'w', encoding='utf-8') as f:
        f.write(table_text(df, kind))
stamp = date.today().strftime('%d.%b.%Y')
upload = os.path.join(OUT, f'Eyster2020_Weil2004_Chuar_Group_{stamp}.txt')
with open(upload, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {upload}: {len(locs)} locations, {len(sites)} sites, {len(ages)} ages rows')

for r in loc_rows:
    if 'pole_lat' in r:
        print(f"  {r['result_name']:40s} pole {r['pole_lat']}/{r['pole_lon']} A95 "
              f"{r['pole_alpha95']} K {r['pole_k']} N {r['pole_n_sites']} rev% "
              f"{r['pole_reversed_perc']} | dir {r['dir_dec']}/{r['dir_inc']} k {r['dir_k']} "
              f"a95 {r['dir_alpha95']}")

# poles (lat, lon, A95, B, N samples); the combined pole is the one used in the
# Swanson-Hysell et al. (2026) compilation; the member poles reproduce Eyster et al. (2020)
# Table 1 (Carbon Butte-Awatubi 13.5/162.8 A95 3.3; Carbon Canyon -2.1/163.7 A95 8.0)
EXPECT = {'combined': (7.7, 163.1, 4.1, 37, 151), 'CBA': (13.5, 162.8, 3.4, 23, 57),
          'CCM': (-2.1, 163.7, 8.0, 14, 94)}
for key, (p, sub) in {'combined': (p_all, comb), 'CBA': (p_cba, cba),
                      'CCM': (p_cc, ccm)}.items():
    got = (round(p['inc'], 1), round(p['dec'], 1), round(p['alpha95'], 1), int(p['n']),
           int(sub['dir_n_samples'].sum()))
    assert got == EXPECT[key], (key, got)
assert (round(d_cba['dec'], 1), round(d_cba['inc'], 1), round(d_cba['k'], 1),
        round(d_cba['alpha95'], 1)) == (98.2, -21.8, 43.4, 4.6), d_cba
print('-I- poles match the compilation and Eyster et al. (2020) Table 1')

sys.path.insert(0, os.path.join(REPO, 'scripts'))
from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True)
