"""Build the MagIC contribution for the Uinta Mountain Group of Weil et al. (2006).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Weil, A. B., Geissman, J. W., & Ashby, J. M. (2006). A new paleomagnetic pole for the
Neoproterozoic Uinta Mountain supergroup, Central Rocky Mountain States, USA. Precambrian
Research 147, 234-259. doi:10.1016/j.precamres.2006.01.017

Dates: Dehler, C., et al. (2023). Precise U-Pb age models refine Neoproterozoic western
Laurentian rift initiation, correlation, and Earth system changes. Precambrian Research 396,
107156. doi:10.1016/j.precamres.2023.107156. Dehler, C. M., et al. (2010). Maximum depositional
age and provenance of the Uinta Mountain Group and Big Cottonwood Formation, northern Utah.
GSA Bulletin 122, 1686-1699. doi:10.1130/B30094.1. Rooney, A. D., et al. (2018). Coupled Re-Os
and U-Pb geochronology of the Tonian Chuar Group, Grand Canyon. GSA Bulletin 130, 1085-1098.
doi:10.1130/B31768.1.

Sources and checks
------------------
- Table 1 (east/west-directed characteristic remanence) and Table 2 (north-directed
  component) are transcribed in `Weil2006_Tables1_2_site_means.tsv` from the text layer of the
  paper and checked against the page images (pp. 242-246), as printed. Table 3 (locality means,
  sampling-area coordinates and locality VGPs, p. 248) and the locality means in the text
  (sections 4.1.1-4.1.11) are used for comparison.
- Dowd Springs fold sites DS 39-44: the tilt-corrected inclinations printed in Table 1 equal the
  alpha95 values of the same rows (11.5, 12.7, 8.5, 14.8, 12.2, 7.9). Restoring the in-situ
  directions about the tabulated bedding reproduces the printed tilt-corrected declinations,
  and the restored inclinations (5.8, -31.8, -12.6, -11.3, -3.4, -17.5) reproduce the published
  locality mean (269.0/-11.9, alpha95 11.5; asserted below). The restored inclinations are
  stored.
- Cross Mountain site CR 9 (alpha95 46.6, k 1.8) is not footnoted in Table 2, but the
  published Cross Mountain mean (352.2/60.2, alpha95 3.3) is reproduced only without it; it is
  entered with result_quality 'b'.
- Stored as printed: the IC-AW-02 dip is printed as a slash followed by 6 (bed_dip left blank);
  FG-GA-01, IC-AW-12 and most Bull Canyon sites do not restore exactly from the tabulated
  bedding (the Bull Canyon site means combine remagnetization circles and set points); some k
  and alpha95 pairs are not mutually consistent (e.g. SC 26, FG-AW-01, DS 43, UMG-GA-01,
  IC-AW-01).
- The locality mean directions are reproduced from the site data (asserted below). Table 3
  prints a Home Mountain Canyon declination of 267.7 (264.7 from its three sites) and a
  Flaming Gorge inclination of -9.8 (-11.0); Table 3 lists Bull Canyon as 10/10, whereas six
  sites are accepted in Table 1 and the locality mean uses six (Fig. 7g: 6/10).
- Coordinates: the paper gives no site coordinates. Sites of the characteristic-remanence
  localities carry the sampling-area coordinates of Table 3 (0.1 deg), except Browne Lake,
  Carter Creek and Dowd Springs: Table 3 prints 40.5 N, 250.5-250.6 E for them, whereas Fig. 2
  places them beside Sheep Creek Canyon (Table 3: 40.9 N, 250.3 E), as do the named
  features. They carry the USGS GNIS coordinates of Browne Lake (40.865 N, 250.186 E) and
  Dowd Spring (40.881 N, 250.307 E), and for Carter Creek, a stream, an approximate point on
  its reach at the Fig. 2 symbol (40.86 N, 250.27 E). The four overprint-only localities,
  which Table 3 does not list, carry approximate coordinates of the named localities.

Structure and interpretation
----------------------------
- Location 'Uinta Mountain Group' has two pole rows: (1) the pole used in the
  Swanson-Hysell et al. (2026) compilation, the Fisher mean of the tilt-corrected VGPs of the
  74 accepted Table 1 sites (sites without the rejection footnotes) of the nine localities
  (Dowd Springs split into flat-lying and folded strata), 496 samples; (2) the authors'
  preferred pole as published (0.8 N, 161.3 E, A95 4.6; GPMDB 9290), the mean of the nine
  locality VGPs of Table 3, with no site count because its confidence limit rests on the
  nine localities.
- One location per sampling locality. The nine characteristic-remanence localities carry the
  locality mean direction (tilt-corrected, majority polarity) without a pole: the localities
  are averaged into one pole by the authors and their VGPs (quoted in the descriptions) are
  not interpreted as recording different times. The four eastern localities with only the
  north-directed component carry no result row.
- The north-directed component (Table 2; in situ and tilt-corrected site rows) is interpreted
  by the authors as a young magnetization acquired after Laramide tilting; its locality means
  are given in the descriptions. Like every site, its rows carry the age of the Uinta
  Mountain Group; the interpreted age of the overprint is given in the descriptions.
- Poles: Fisher mean of polarity-unified site VGPs (calculated at the site coordinates) in the
  ca. 160 E hemisphere of the authors. Directions and poles are not corrected for inclination
  shallowing.
- No structured field-test columns: the local fold tests are not significant (Dowd Springs) or
  negative (Sheep Creek), and the paper reports no reversal test.
- Age: 759 Ma (743-767 Ma), between the CA-ID-TIMS detrital zircon maximum depositional age of
  the basal Uinta Mountain Group (766.3 +/- 0.5 Ma; Dehler et al., 2023) and the Re-Os date on
  the basal Awatubi Member of the Kwagunt Formation, Chuar Group (751.0 +/- 7.6 Ma; Rooney et
  al., 2018), with which the Red Pine Shale is correlated (all sampled localities are below the
  Red Pine Shale). The ages table holds these dates and the SHRIMP maximum depositional age of
  Dehler et al. (2010).
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

WEIL = '10.1016/j.precamres.2006.01.017'
DEHLER23 = '10.1016/j.precamres.2023.107156'
DEHLER10 = '10.1130/B30094.1'
ROONEY = '10.1130/B31768.1'
CHRM, NORTH = 'characteristic', 'north-directed overprint'
UMG = 'Uinta Mountain Group'
AGE_CHRM = (759, 743, 767)

# paper locality name -> (MagIC location, state, lat, lon); Table 3 sampling-area coordinates
# for the characteristic-remanence localities, approximate coordinates of the named
# localities for the overprint-only localities
LOCALITIES = {
    'Sheep Creek Canyon': ('Uinta Mountain Group, Sheep Creek Canyon', 'Utah', 40.9, 250.3),
    # Table 3 prints 40.5 N, 250.5-250.6 E for these three; Fig. 2 and the named features
    # place them beside Sheep Creek Canyon: USGS GNIS Browne Lake 40.865 N, 109.814 W and
    # Dowd Spring 40.881 N, 109.693 W; Carter Creek (a stream) on its reach at the Fig. 2
    # symbol (approximate)
    'Browne Lake': ('Uinta Mountain Group, Browne Lake', 'Utah', 40.865, 250.186),
    'Carter Creek': ('Uinta Mountain Group, Carter Creek', 'Utah', 40.86, 250.27),
    'Dowd Springs non-fold test rocks': ('Uinta Mountain Group, Dowd Springs flat-lying strata', 'Utah', 40.881, 250.307),
    'Dowd Springs fold test rocks': ('Uinta Mountain Group, Dowd Springs fold', 'Utah', 40.881, 250.307),
    'Flaming Gorge Section': ('Uinta Mountain Group, Flaming Gorge', 'Utah', 40.6, 250.7),
    'Bull Canyon': ('Uinta Mountain Group, Bull Canyon', 'Colorado', 40.9, 251.2),
    'Home Montain Canyon': ('Uinta Mountain Group, Home Mountain Canyon', 'Utah', 41.0, 250.8),
    'Sparks-Talamantes Creek': ('Uinta Mountain Group, Sparks-Talamantes Creek', 'Colorado', 40.9, 251.2),
    'Irish Canyon Section': ('Uinta Mountain Group, Irish Canyon', 'Colorado', 40.8617, 251.1947),
    'Cross Mountain': ('Uinta Mountain Group, Cross Mountain', 'Colorado', 40.5364, 251.6078),
    'Lone Mountain': ('Uinta Mountain Group, Lone Mountain', 'Colorado', 40.5757, 251.5049),
    'Juniper Mountain': ('Uinta Mountain Group, Juniper Mountain', 'Colorado', 40.5016, 251.8045),
}
NOT_IN_TABLE3 = {'Irish Canyon Section', 'Cross Mountain', 'Lone Mountain', 'Juniper Mountain'}
RELOCATED = {'Browne Lake': (40.5, 250.5), 'Carter Creek': (40.5, 250.6),
             'Dowd Springs non-fold test rocks': (40.5, 250.6),
             'Dowd Springs fold test rocks': (40.5, 250.6)}
# Weil et al. (2006) Table 3 and text (dec, inc, a95, VGP lat, VGP lon), tilt-corrected
PUBLISHED = {
    'Sheep Creek Canyon': (96.2, -5.6, 4.7, 6.5, 158.4), 'Browne Lake': (84.3, -8.3, 13.1, -1.6, 167.4),
    'Carter Creek': (266.0, -8.0, 5.4, -5.7, 160.1),
    'Dowd Springs non-fold test rocks': (101.8, -6.9, 12.1, 11.2, 155.6),
    'Dowd Springs fold test rocks': (269.0, -11.9, 11.5, -4.7, 156.7),
    'Flaming Gorge Section': (267.1, -9.8, 4.5, -5.4, 158.9), 'Bull Canyon': (272.0, 9.0, 5.9, 4.5, 163.3),
    'Home Montain Canyon': (267.7, 6.7, 12.5, 0.5, 164.8),
    'Sparks-Talamantes Creek': (87.9, -9.9, 8.2, 1.7, 166.3),
}
# north-directed component, in situ (text): dec, inc, a95
PUBLISHED_OVP = {
    'Sheep Creek Canyon': (40.2, 63.6, 5.5), 'Irish Canyon Section': (6.0, 60.5, 2.1),
    'Cross Mountain': (352.2, 60.2, 3.3), 'Lone Mountain': (358.3, 60.7, 2.7),
    'Juniper Mountain': (7.5, 50.6, 12.5)}
DS_FOLD = ['DS 39', 'DS 40', 'DS 41', 'DS 42', 'DS 43', 'DS 44']
BC_CIRCLES = ['BC 01', 'BC 02', 'BC 04', 'BC 05', 'BC 06', 'BC 07', 'BC 08', 'BC 10']
FLAG_REASON = {
    'a': 'rejected by the authors: direction more than two angular standard deviations from the locality mean',
    'b': 'rejected by the authors: alpha95 > 15 degrees',
}
MC_SITE = 'LP-DIR-T:LP-DIR-CHEM:DE-BFL:DE-FM'
MC_SITE_BC = 'LP-DIR-T:LP-DIR-CHEM:DE-BFL:DE-BFP-G:DE-FM-LP'
MC_LOC = 'LP-DIR-T:LP-DIR-CHEM:LP-IRM:LP-IRM-DCD:LP-IRM-3D:DE-BFL:DE-FM'
LITH = 'Sandstone:Siltstone'
GEOG = dict(continent_ocean='North America', country='United States of America',
            terranes='Laurentia')


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def r1(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 1)


T = pd.read_csv(os.path.join(HERE, 'Weil2006_Tables1_2_site_means.tsv'), sep='\t',
                comment='#', dtype=str, keep_default_na=False)
T = T[T.dec_geo != '-'].copy()           # footnote c: no interpretable results

rows = []
for _, t in T.iterrows():
    loc, state, lat, lon = LOCALITIES[t.locality]
    comp = CHRM if t.table == 'T1' else NORTH
    strike, dip = num(t.strike), num(t.dip)
    n, n0 = [int(x) for x in t.N_N0.split('/')]
    gd, gi, a95, k = num(t.dec_geo), num(t.inc_geo), num(t.alpha95), num(t.k)
    td, ti = num(t.dec_tc), num(t.inc_tc)
    quality = 'b' if t.flag in ('a', 'b') else 'g'
    notes = []
    if t.site in DS_FOLD and t.table == 'T1':
        notes.append(f'Table 1 prints a tilt-corrected inclination of {t.inc_tc}, equal to the '
                     f'alpha95 of the row; the tilt-corrected inclination stored is restored '
                     f'from the in-situ direction and bedding, which reproduces the printed '
                     f'declination and the published locality mean.')
        ti = round(pmag.dotilt(gd, gi, (strike + 90) % 360, dip)[1], 1)
    if t.site == 'IC-AW-02':
        notes.append('The bedding dip is printed as a slash followed by 6 and is not stored.')
    if t.site == 'CR 9' and t.table == 'T2':
        quality = 'b'
        notes.append('[excluded: alpha95 46.6 degrees; not footnoted in Table 2 but not in '
                     'the published Cross Mountain mean]')
    if t.flag in FLAG_REASON:
        notes.append(f'[excluded: {FLAG_REASON[t.flag]}]')
    if t.site in ('SC 13', 'DS 29', 'CR 22', 'LN 17'):
        notes.append('The site has two tabulated mean directions.')
    if strike is None:
        strike = 0.0      # flat-lying strata
    if t.table == 'T1':
        desc = ('Hematite-cemented fine- to medium-grained sandstone or siltstone bed(s) of '
                'the Neoproterozoic Uinta Mountain Group; east- or west-directed, shallow '
                'characteristic remanence (Table 1).')
    else:
        desc = ('Hematite-cemented sandstone or siltstone bed(s) of the Neoproterozoic Uinta '
                'Mountain Group carrying the north-directed, moderate to steep, normal-polarity '
                'magnetization (Table 2) that the authors interpret as a young, post-tilting '
                '(Laramide or younger) overprint.')
    if t.locality in NOT_IN_TABLE3:
        desc += ' Site coordinate is an approximate locality coordinate (not given in the paper).'
    elif t.locality in RELOCATED:
        la3, lo3 = RELOCATED[t.locality]
        desc += (' Site coordinate is that of the named locality (USGS GNIS; approximate for '
                 'Carter Creek, placed on the creek at the sampling-area symbol of Fig. 2). '
                 f'Table 3 gives {la3} N, {lo3} E, whereas Fig. 2 places the sampling area '
                 'beside Sheep Creek Canyon, ca. 0.37 deg further north.')
    else:
        desc += ' Site coordinate is the sampling-area coordinate of the locality (Table 3).'
    if notes:
        desc += ' ' + ' '.join(notes)
    age, lo, hi = AGE_CHRM
    mc = MC_SITE_BC if t.site in BC_CIRCLES and t.table == 'T1' else MC_SITE
    for tc, d, i in ((0, gd, gi), (100, td, ti)):
        plon, plat, dp, dm = pmag.dia_vgp(d, i, a95 if a95 is not None else np.nan, lat, lon)
        rows.append({
            'site': t.site, 'location': loc, 'result_type': 'i', 'result_quality': quality,
            'method_codes': mc, 'citations': WEIL, 'geologic_classes': 'Sedimentary',
            'geologic_types': 'Sediment Layer', 'lithologies': LITH,
            'bed_dip': dip, 'bed_dip_direction': (strike + 90) % 360 if dip is not None else None,
            'lat': lat, 'lon': lon, 'age': age, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
            'dir_comp_name': comp, 'dir_tilt_correction': tc, 'dir_dec': d, 'dir_inc': i,
            'dir_alpha95': a95, 'dir_k': k, 'dir_n_samples': n, 'dir_n_total_samples': n0,
            'vgp_lat': r1(plat), 'vgp_lon': r1(plon), 'vgp_dp': r1(dp), 'vgp_dm': r1(dm),
            'description': desc, '_paper_loc': t.locality})

sites = pd.DataFrame(rows)
for c in ['bed_dip', 'bed_dip_direction', 'age', 'dir_alpha95', 'dir_k', 'dir_n_samples',
          'dir_n_total_samples']:
    sites[c] = pd.to_numeric(sites[c])


def unify(pairs, ref):
    rc = pmag.dir2cart([ref[0], ref[1], 1.])
    return [p if np.dot(pmag.dir2cart([p[0], p[1], 1.]), rc) > 0 else
            [(p[0] + 180) % 360, -p[1]] for p in pairs]


def pole(sub, ref=(160, 0)):
    """Fisher mean of polarity-unified VGPs, in the hemisphere of ref (lon, lat)."""
    return pmag.fisher_mean(unify(sub[['vgp_lon', 'vgp_lat']].values.tolist(), ref))


def majority_dir(sub):
    """Fisher mean of site directions in the majority (observed) polarity, and the
    percentage of sites of the minority polarity."""
    pairs = sub[['dir_dec', 'dir_inc']].values.tolist()
    if len(pairs) == 1:
        r = sub.iloc[0]
        return {'dec': r.dir_dec, 'inc': r.dir_inc, 'alpha95': r.dir_alpha95, 'k': r.dir_k,
                'n': 1}, 0.0
    m0 = pmag.doprinc(pairs)
    rc = pmag.dir2cart([m0['dec'], m0['inc'], 1.])
    same = np.array([np.dot(pmag.dir2cart([d, i, 1.]), rc) > 0 for d, i in pairs])
    if same.sum() < len(same) / 2:
        same = ~same
        m0 = {'dec': (m0['dec'] + 180) % 360, 'inc': -m0['inc']}
    return pmag.fisher_mean(unify(pairs, (m0['dec'], m0['inc']))), \
        round(100 * (~same).sum() / len(same), 1)


good = sites[sites.result_quality == 'g']
chrm_tc = good[(good.dir_comp_name == CHRM) & (good.dir_tilt_correction == 100)]
north_is = good[(good.dir_comp_name == NORTH) & (good.dir_tilt_correction == 0)]
assert len(chrm_tc) == 74 and int(chrm_tc.dir_n_samples.sum()) == 496

# locality means reproduce Weil et al. (2006); Home Mountain Canyon: Table 3 prints D 267.7
for pl, pub in PUBLISHED.items():
    d, _ = majority_dir(chrm_tc[chrm_tc._paper_loc == pl])
    tol = 3.1 if pl == 'Home Montain Canyon' else 1.3
    assert pmag.angle([d['dec'], d['inc']], list(pub[:2]))[0] < tol, (pl, d)
d_dsf, _ = majority_dir(chrm_tc[chrm_tc._paper_loc == 'Dowd Springs fold test rocks'])
assert pmag.angle([d_dsf['dec'], d_dsf['inc']], [269.0, -11.9])[0] < 0.1 and \
    round(d_dsf['alpha95'], 1) == 11.5, d_dsf
for pl, pub in PUBLISHED_OVP.items():
    d, _ = majority_dir(north_is[north_is._paper_loc == pl])
    assert pmag.angle([d['dec'], d['inc']], list(pub[:2]))[0] < 0.2, (pl, d)


def locrow(location, sub, state, desc, rname=None, with_pole=False, age=AGE_CHRM):
    d, rev = majority_dir(sub)
    a, lo, hi = age
    r = {'location': location, 'location_type': 'Region', 'result_name': rname,
         'result_type': 'a', 'result_quality': 'g',
         'sites': ':'.join(dict.fromkeys(sub.site)),
         'method_codes': MC_LOC + (':DE-BFP-G:DE-FM-LP' if 'Bull Canyon' in location
                                   or location == UMG else ''),
         'citations': WEIL, 'geologic_classes': 'Sedimentary', 'lithologies': LITH,
         'lat_s': sub.lat.min(), 'lat_n': sub.lat.max(), 'lon_w': sub.lon.min(),
         'lon_e': sub.lon.max(), 'age': a, 'age_low': lo, 'age_high': hi, 'age_unit': 'Ma',
         'dir_tilt_correction': 100, 'dir_dec': r1(d['dec']), 'dir_inc': r1(d['inc']),
         'dir_alpha95': r1(d['alpha95']), 'dir_k': r1(d['k']), 'dir_n_sites': int(d['n']),
         'state_province': state, 'description': desc}
    if with_pole:
        p = pole(sub)
        r.update({'method_codes': r['method_codes'] + ':DE-VGP', 'pole_comp_name': CHRM,
                  'pole_lat': r1(p['inc']), 'pole_lon': r1(p['dec']),
                  'pole_alpha95': r1(p['alpha95']), 'pole_k': r1(p['k']),
                  'pole_n_sites': int(p['n']), 'pole_reversed_perc': rev})
    r.update(GEOG)
    return r


loc_rows = [locrow(
    UMG, chrm_tc, None, rname='Uinta Mountain Group pole', with_pole=True,
    desc=("Paleomagnetic pole for the Neoproterozoic (Tonian) Uinta Mountain Group, Uinta "
          "Mountains of northeastern Utah and northwestern Colorado (Weil et al., 2006): the "
          "Fisher mean of the tilt-corrected virtual geomagnetic poles of all accepted sites "
          "carrying the east- or west-directed, shallow, dual-polarity characteristic "
          "remanence at nine sampling localities (Dowd Springs split into flat-lying and "
          "folded strata). The remanence is carried by fine-grained pigmentary hematite, "
          "unblocks between ca. 660 and 680 C, and is interpreted by the authors as a near "
          "primary remanence; they note that, except at Carter Creek and Dowd Springs, each "
          "locality has a single polarity, that the direction resembles Late Cambrian to "
          "earliest Ordovician directions for North America, and that they have no direct, "
          "field-based proof of its age. Local fold tests are not significant because bedding "
          "strikes are sub-parallel to the magnetization. The authors' published pole (0.8 N, "
          "161.3 E, A95 4.6) is the mean of the nine locality VGPs. The tilt-corrected "
          "inclinations of Dowd Springs sites DS 39-44 are restored as described in the site "
          "rows. Directions are not corrected for inclination shallowing. The sites span Utah "
          "and Colorado. The age bracket is from the CA-ID-TIMS detrital zircon maximum "
          "depositional age of the basal eastern Uinta Mountain Group (Dehler et al., 2023) "
          "and the Re-Os date on the basal Awatubi Member, Kwagunt Formation, Chuar Group "
          "(Rooney et al., 2018), with which the Red Pine Shale is correlated; all sampled "
          "localities are stratigraphically below the Red Pine Shale."))]
loc_rows[0]['citations'] = f'{WEIL}:{DEHLER23}:{DEHLER10}:{ROONEY}'

LOC_DESC = {
    'Sheep Creek Canyon': 'Sheep Creek Canyon, north flank of the Uinta Mountains, Utah, where strata steepen to overturned toward the north-bounding Uinta fault; the characteristic remanence is exclusively east-directed. The locality also carries the north-directed component (20 sites; in situ 40.2/63.6, alpha95 5.5), whose dispersion increases significantly on untilting and whose northeast in-situ declinations the authors suspect may be a composite of the characteristic remanence and a young magnetization.',
    'Browne Lake': 'Browne Lake, north-central Uinta Mountains, Utah; one of four sites accepted.',
    'Carter Creek': 'Carter Creek, north-central Uinta Mountains, Utah, ca. 50 m stratigraphically below Browne Lake; some samples are of the opposite polarity.',
    'Dowd Springs non-fold test rocks': 'Dowd Springs, north-central Uinta Mountains, Utah: very gently tilted strata (five of ten sites accepted).',
    'Dowd Springs fold test rocks': 'Dowd Springs, north-central Uinta Mountains, Utah: both limbs of an upright open anticline (fold test not significant at the 90% level). Tilt-corrected inclinations restored as described in the site rows.',
    'Flaming Gorge Section': 'Flaming Gorge Dam area, north-central Uinta Mountains, Utah (16 of 18 sites).',
    'Bull Canyon': 'Bull Canyon, Browns Park, eastern Uinta Mountains, Colorado; site means other than BC 03 and BC 09 from combined remagnetization circles and set points. Six of ten sites are accepted (Table 3 lists 10/10; the locality mean uses six).',
    'Home Montain Canyon': 'Home Mountain Canyon, basal Uinta Mountain Group just west of the Utah-Colorado border.',
    'Sparks-Talamantes Creek': 'Talamantes Creek west of Sparks, Colorado, eastern Uinta Mountains (9 of 12 sites).',
}
for pl, text in LOC_DESC.items():
    loc, state = LOCALITIES[pl][:2]
    pub = PUBLISHED[pl]
    lat_txt = f'{abs(pub[3])} {"S" if pub[3] < 0 else "N"}'
    loc_rows.append(locrow(
        loc, chrm_tc[chrm_tc._paper_loc == pl], state,
        rname=f"{loc.split(', ')[1]} characteristic remanence mean direction",
        desc=(f"Mean tilt-corrected characteristic-remanence direction of the accepted sites "
              f"of the Uinta Mountain Group at {text} The authors' VGP of the locality mean "
              f"(Table 3) is {lat_txt}, {pub[4]} E; the locality VGPs are averaged in the "
              f"Uinta Mountain Group pole of the authors.")))

OVP_DESC = {
    'Irish Canyon Section': 'Irish Canyon, Browns Park, eastern Uinta Mountains, Colorado (22 of 25 sites; three without interpretable results)',
    'Cross Mountain': 'Cross Mountain, southeastern Uinta Mountains, Colorado, where the north-directed magnetization is the only component and NRM directions of incompletely demagnetized sites approximate it (18 accepted sites)',
    'Lone Mountain': 'Lone Mountain, southeastern Uinta Mountains, Colorado (12 accepted sites; limited demagnetization)',
    'Juniper Mountain': 'Juniper Mountain, easternmost Uinta Mountains, Colorado (5 accepted sites)',
}
for pl, text in OVP_DESC.items():
    loc, state, lat, lon = LOCALITIES[pl]
    sub = sites[sites._paper_loc == pl]
    dd, ii, aa = PUBLISHED_OVP[pl]
    r = {'location': loc, 'location_type': 'Region', 'sites': ':'.join(dict.fromkeys(sub.site)),
         'method_codes': MC_LOC, 'citations': WEIL, 'geologic_classes': 'Sedimentary',
         'lithologies': LITH, 'lat_s': lat, 'lat_n': lat, 'lon_w': lon, 'lon_e': lon,
         'age': AGE_CHRM[0], 'age_low': AGE_CHRM[1], 'age_high': AGE_CHRM[2], 'age_unit': 'Ma',
         'state_province': state,
         'description': (f"Sites of the Uinta Mountain Group at {text}, which carry only the "
                         f"north-directed, moderate to steep, normal-polarity magnetization "
                         f"(in-situ locality mean {dd}/{ii}, alpha95 {aa}; Weil et al., 2006). "
                         f"The authors interpret this hematite-carried magnetization as young "
                         f"and acquired after Laramide tilting, no older than mid-Cretaceous; "
                         f"the age given is that of the Uinta Mountain Group. "
                         f"Coordinates are approximate.")}
    r.update(GEOG)
    loc_rows.append(r)

# the authors' preferred pole, as published: Fisher mean of the nine locality VGPs of Table 3
# (GPMDB 9290); the confidence limit rests on the nine localities, so pole_n_sites is blank
n_rev_loc = sum(0 < v[0] < 180 for v in PUBLISHED.values())      # east-directed localities
umg_auth = {k: loc_rows[0][k] for k in (
    'location', 'location_type', 'result_type', 'result_quality', 'sites', 'method_codes',
    'citations', 'geologic_classes', 'lithologies', 'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age',
    'age_low', 'age_high', 'age_unit')}
umg_auth.update({
    'result_name': 'Uinta Mountain Group pole of Weil et al. (2006), mean of locality VGPs',
    'pole_comp_name': CHRM, 'pole_lat': 0.8, 'pole_lon': 161.3, 'pole_alpha95': 4.6,
    'pole_reversed_perc': round(100 * min(n_rev_loc, 9 - n_rev_loc) / 9, 1),
    'external_database_ids': 'GPMDB[9290]', **GEOG,
    'description': (
        "Paleomagnetic pole for the Uinta Mountain Group as published by Weil et al. (2006): "
        "the Fisher mean of the virtual geomagnetic poles of the nine locality mean "
        "directions of their Table 3 (Dowd Springs split into flat-lying and folded strata). "
        "Values are as published; the confidence limit is calculated from the nine locality "
        "VGPs, so no site count is given (the sites listed are the 74 accepted sites of the "
        "localities). Table 3 prints locality means that differ from those calculated from the "
        "site table for Home Mountain Canyon (declination 267.7, 264.7 from its sites) and "
        "Flaming Gorge (inclination -9.8, -11.0 from its sites); with the locality means "
        "calculated from the sites and the coordinates of this contribution, the mean of the "
        "locality VGPs is 0.5 N, 161.3 E, A95 4.8.")})
loc_rows.insert(1, umg_auth)
locs = pd.DataFrame(loc_rows)
sites = sites.drop(columns='_paper_loc')

AGE_ROWS = [
    (766.3, 0.25, 'GM-UPB', DEHLER23,
     'U-Pb CA-ID-TIMS detrital zircon maximum depositional age, youngest population (n = 5) '
     'of the Uinta Mountain Group, ca. 60 m above the base of the Jesse Ewing Canyon '
     'Formation, eastern Uinta Mountains (sample 17CMD-JEC-5), 766.3 +/- 0.5 Ma (Dehler et '
     'al., 2023; uncertainty treated as 2 sigma as for the other CA-ID-TIMS dates of the '
     'paper; the text also quotes 766.3 +/- 0.9 Ma for n = 3).'),
    (766.4, 4.8, 'GM-SHRIMP', DEHLER10,
     'U-Pb SHRIMP detrital zircon concordia age of the youngest grains in the Outlaw Trail '
     'formation, eastern Uinta Mountain Group (sample SCUMG-9), 766.4 +/- 4.8 Ma (Dehler et '
     'al., 2010; the concordia plot labels the uncertainty 1 sigma); a maximum depositional '
     'age.'),
    (751.0, 3.8, 'GM-REOS-I', ROONEY,
     'Re-Os isochron (model 1) age of marcasite nodules in the basal Awatubi Member, Kwagunt '
     'Formation, Chuar Group, Grand Canyon, 751.0 +/- 7.6 Ma (2 sigma; Rooney et al., 2018). '
     'Not a date on the Uinta Mountain Group: the upper Uinta Mountain Group (Red Pine '
     'Shale) is correlated with the Kwagunt Formation by vase-shaped microfossils, and the '
     'sampled localities lie below the Red Pine Shale.'),
]
ages = pd.DataFrame([{
    'location': UMG, 'age': a, 'age_sigma': s, 'age_unit': 'Ma', 'method_codes': mc,
    'citations': cit, 'timescale_eon': 'Proterozoic', 'timescale_era': 'Neoproterozoic',
    'timescale_period': 'Tonian', 'description': desc} for a, s, mc, cit, desc in AGE_ROWS])


def fmt(df):
    df = df.dropna(axis=1, how='all').copy()
    for c in df.columns:
        if df[c].dtype.kind in 'fiu':
            df[c] = df[c].map(lambda x: '' if pd.isna(x) else f'{float(x):.10g}')
    return df.fillna('')


def table_text(df, kind):
    return f'tab delimited\t{kind}\n' + fmt(df).to_csv(sep='\t', index=False)


TABLES = (('locations', locs), ('sites', sites), ('ages', ages))
for kind, df in TABLES:
    with open(os.path.join(OUT, f'{kind}.txt'), 'w', encoding='utf-8') as f:
        f.write(table_text(df, kind))
stamp = date.today().strftime('%d.%b.%Y')
upload = os.path.join(OUT, f'Weil2006_Uinta_Mountain_Group_{stamp}.txt')
with open(upload, 'w', encoding='utf-8') as f:
    f.write('>>>>>>>>>>\n'.join(table_text(df, kind) for kind, df in TABLES))
print(f'-I- wrote {upload}: {len(locs)} locations, {len(sites)} sites, {len(ages)} ages rows')

p_all = pole(chrm_tc)
r = loc_rows[0]
print(f"  {UMG} pole {r['pole_lat']}/{r['pole_lon']} A95 {r['pole_alpha95']} K {r['pole_k']} "
      f"N {r['pole_n_sites']} rev% {r['pole_reversed_perc']} | dir {r['dir_dec']}/"
      f"{r['dir_inc']} k {r['dir_k']} a95 {r['dir_alpha95']}")

# pole used in the Swanson-Hysell et al. (2026) compilation (lat, lon, A95, B, N samples)
EXPECT = (1.2, 159.9, 2.2, 74, 496)
got = (round(p_all['inc'], 1), round(p_all['dec'], 1), round(p_all['alpha95'], 1),
       int(p_all['n']), int(chrm_tc['dir_n_samples'].sum()))
assert got == EXPECT, got
print('-I- pole matches the compilation')

sys.path.insert(0, os.path.join(REPO, 'scripts'))
from validate_magic_contribution import validate_upload_file  # noqa: E402
validate_upload_file(upload, tables=['locations', 'sites', 'ages'], server=True)
