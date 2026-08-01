"""
Build a MagIC 3.0 sites.txt + locations.txt for the ca. 1108 Ma Nipigon (Logan)
sills paleomagnetic pole from the site-level compilation of the Midcontinent
Rift sills and intrusions of the Nipigon Embayment and Thunder Bay district
assembled by E. J. Iloranta (2026; `Iloranta2026_TB_sills_source.csv`).

The compilation gathers published site-mean directions from five studies:

  DuBois (1962)                 GSC Bulletin 71     10.4095/100589
  Robertson & Fahrig (1971)     Can. J. Earth Sci.  10.1139/e71-125
  Pesonen (1979)                Bull. Geol. Soc. Fi 10.17741/bgsf/51.1-2.004
  Middleton et al. (2004)       J. Geophys. Res.    10.1029/2003JB002581
  Borradaile & Middleton (2006) Precambrian Res.    10.1016/j.precamres.2005.10.007
  Piispa et al. (in prep.)      Thunder Intrusion   (no DOI yet)

VGPs in the source compilation were recomputed by Iloranta from the site mean
directions and site coordinates rather than adopted from the original papers.
That recomputation was verified here against pmagpy `pmag.dia_vgp` for all 60
entries with directions: agreement is exact to 0.003 deg in VGP latitude and
0.02 deg in VGP longitude. Reverse-polarity VGPs are reported as the antipodal
north pole, the convention used throughout this repository.

Two levels of averaging
-----------------------
Many of the same sills were sampled by more than one study, and Middleton et
al. (2004) list the same locality more than once under different demagnetization
treatments. Independent cooling units are therefore recovered in two stages:

  Stage A  Collapse the 15 entries of Middleton et al. (2004) Table 1 to their
           9 distinct localities (Disraeli, SW, Hele 1, Hele 2, Havoc, Terry
           Fox, Geikie, Rift, Fox Mountain). Table 1 reports the same locality
           under thermal (T), low-temperature (L) and combined (L & T)
           treatments; these are repeat determinations, not separate sites.
  Stage B  Merge sills sampled by more than one study, using the cross-study
           sill correlations of the source compilation ("Sill/dyke name for
           dublicates"), e.g. Red Rock North = DuBois N12-15 + N29 BCT +
           Robertson & Fahrig S13 + Pesonen R46, together with the two further
           correlations established here from the source publications
           (`ADDITIONAL_MERGES`). Three further candidates that the
           publications do not resolve are listed in `CANDIDATE_MERGES` and are
           not applied; the notebook reports the pole's sensitivity to them.

A caveat on "one VGP per sill": the Logan sills are laterally extensive sheets
that transgress from one bedding plane to another (Robertson & Fahrig, 1971),
so widely separated exposures may belong to a single cooling unit. Neither
Robertson & Fahrig (1971) nor Pesonen (1979) describe their sites, and their
site coordinates here are digitized from the site maps of those papers, so
correlation across studies rests on the sill names of the source compilation
(traced from the Ontario geological map) plus site proximity. Sills sampled at
separated exposures without a mapped correlation are therefore still counted
more than once, and the unit count of 37 is an upper bound on the number of
independent cooling units.

Within a merged unit the contributing VGPs are averaged with a Fisher mean, and
the unit mean direction is the Fisher mean of the contributing directions. Each
row of the output `sites.txt` is one such independent unit; the pole is the
Fisher mean of the unit VGPs.

Exclusions (`result_quality = 'b'`; retained in the contribution, dropped by
`pole_tools.load_magic_sites`)
------------------------------------------------------------------------------
  * Pillar Lake Lava (4 entries, Borradaile & Middleton 2006). Lavas of the
    Nipigon Embayment whose age relative to the sills is unconstrained; excluded
    from this sill pole.
  * Normal-polarity entries (Middleton et al. 2004 sites 6, 13, 15: Hele 2 T,
    Hele 2 L & T, Terry Fox T). Each is the alternative-treatment counterpart of
    a reverse-polarity determination of the same locality, and Middleton et al.
    themselves describe the Terry Fox thermal result as spurious and the Hele 2
    site as only yielding the SE-up direction when thermal demagnetization
    followed low-temperature cycling. The ca. 1108 Ma sills are reverse polarity;
    these normal directions are interpreted as incompletely removed overprints.
  * Two entries flagged "Excluded!" in the source compilation (Robertson &
    Fahrig S11, which has no direction, and S4, k = 2.4).

Notes and caveats carried forward for review
--------------------------------------------
  * `result_quality` in the source compilation is not used here as an
    acceptance criterion; its assignment does not track n, k or alpha95 in an
    obvious way and needs the compiler's clarification. It is preserved in
    `accepted_results.csv`.
  * Middleton et al. (2004) site 11 (Rift) has no reported coordinates; the
    study mean location (49.0 N, 271.0 E, their Table 1 footnote) is used.
  * The `Area` labels of the source compilation are used for the regional
    comparison but are inconsistent in places (e.g. DuBois SN5-15, at Current
    River in Thunder Bay, is labelled "Pigeon"). Regions are therefore assigned
    by grouping Area into "Nipigon Embayment" and "Thunder Bay-Logan".
  * The Seagull Intrusion (Borradaile & Middleton, 2006) is dated at
    1112 +/- 2.4 Ma (Hart & Whaley, 2005) in that paper and so may be slightly
    older than the sills; the source compilation instead notes 1107.8 +/- 1.4 Ma.
    It is retained but the sensitivity to its removal is reported in the notebook.
  * Robertson & Fahrig (1971) S5 and S6 ("Colville Lake") are flagged in the
    source compilation as possibly belonging to the younger ca. 1100.8 Ma
    McIntyre diabase. They are reverse polarity and are retained; the notebook
    reports the sensitivity to their removal.
  * The compilation assigns Pesonen R47 and R49 identical coordinates. Figure 1
    of the Geological Survey of Finland report version of Pesonen (1979) plots
    them as separate sites with R49 slightly north of R47; the coordinates
    should be separated in the source compilation.
  * The compilation merges Pesonen R45 into the Terry Fox unit on the basis of
    the Ontario geological map. Pesonen's Figure 1 places R45 near Pass Lake at
    the neck of the Sibley Peninsula, ~27 km northeast of the Current River /
    Terry Fox exposure sampled by DuBois, Robertson & Fahrig and Middleton et
    al., so the merge asserts that one sill sheet is continuous over that
    distance. Plausible for these sills, but it should be confirmed.
  * This is a sills-only pole. Pesonen's reversed Thunder Bay *dikes* (R22-R38,
    his Table 6, pole 48.0 N / 212.2 E) are a separate result and are not in the
    compilation, whereas the Lulea Working Group grand mean that this entry
    replaces did include them.

Age: the Nipigon (Logan) sills cluster near 1105-1109 Ma. Nipigon sills
1108.2 +/- 0.9 Ma (zircon upper intercept recalculated from Davis & Sutcliffe,
1985); Bleeker et al. (2020) report 1106.3 +/- 2.0 Ma for the main Logan Sill at
Mount McKay and 1105.5 +/- 3.0 Ma for the Inspiration sill. Pole age 1108 +/- 2 Ma.
MagIC data model v3.0.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = Path(__file__).parent
OUT = HERE.parent
SOURCE = HERE / 'Iloranta2026_TB_sills_source.csv'

AGE, AGE_LOW, AGE_HIGH = '1108', '1106', '1110'

CITATION = {
    'DuBois, 1962': '10.4095/100589',
    'Robertson and Fahrig, 1971': '10.1139/e71-125',
    'Pesonen, 1979': '10.17741/bgsf/51.1-2.004',
    'Middleton et al. 2004': '10.1029/2003JB002581',
    'Borradaile & Middleton 2006': '10.1016/j.precamres.2005.10.007',
    'Piispa et al. manuscript in preparation': 'Piispa et al., in preparation',
}

# Stage A: Middleton et al. (2004) Table 1 site numbers -> sampling locality.
# The source compilation's site codes M1-M15 are Middleton's Table 1 site numbers.
MIDDLETON_LOCALITY = {
    'M1': 'Disraeli', 'M2': 'Middleton SW', 'M3': 'Middleton SW',
    'M4': 'Hele 1', 'M5': 'Hele 1', 'M6': 'Hele 2', 'M7': 'Havoc',
    'M8': 'Havoc', 'M9': 'Terry Fox', 'M10': 'Geikie', 'M11': 'Rift',
    'M12': 'Fox Mountain', 'M13': 'Hele 2', 'M14': 'Hele 2', 'M15': 'Terry Fox',
}

# Middleton et al. (2004) Table 1 footnote: mean site location 49 N, 271 E.
MIDDLETON_MEAN_LATLON = (49.0, 271.0)

# Readable names for units whose source compilation code is an abbreviation.
UNIT_RENAME = {
    'SI': 'Seagull Intrusion',
    'IS': 'Inspiration Sills',
    'TI': 'Thunder Intrusion',
}

# Stage B additions established from the source publications rather than from the
# cross-study sill correlations of the source compilation. Applied.
#
#   N27-29 -> Red Rock North
#     DuBois (1962) Table XIII places N27, N28 and N29 at one cliff exposure
#     above the CPR and CNR tracks at Red Rock: "N27's and N28's from sediment
#     immediately below diabase sill ... N29's from lower part of the sill".
#     The compilation splits them into "N29 BCT" (the baked sediment, N27-N28)
#     and "N27-29" (the sill itself, N29). A sill and the sediment it baked
#     record the same cooling event and are not independent, so both belong to
#     the Red Rock North unit, to which "N29 BCT" was already assigned.
#
#   S7 -> Mt Mackay
#     Robertson & Fahrig (1971) give no site descriptions, but their site S7
#     plots (as digitized from their Figure 1 in the source compilation) 0.33 km
#     from DuBois's Mount Mackay site, the type exposure of the lower Mount
#     McKay sill at Fort William dated by Bleeker et al. (2020) at
#     1106.3 +/- 2.0 Ma. DuBois (1962, Table XV) obtained badly scattered
#     directions there, attributed them to a superimposed random component
#     (possibly lightning, given the hill's elevation) and resorted to
#     great-circle fits; S7 (k = 92, alpha95 = 10) is the better determination
#     of the same sill.
ADDITIONAL_MERGES = {
    'N27-29': 'Red Rock North',
    'S7': 'Mt Mackay',
}

# Further same-sill correlations suggested by site proximity, but NOT established
# by the publications and therefore NOT applied. The notebook reports the
# sensitivity of the pole to each. Both need the Ontario geological map (or the
# compiler's field knowledge) to resolve, because the Logan sills are laterally
# extensive sheets and Robertson & Fahrig (1971) describe none of their sites.
#
#   S8 -> N19-21          Robertson & Fahrig S8 lies 3.5 km from DuBois's N19-21
#                         sill on the Highway 11 scarp south of Orient Bay.
#   N16-18 -> Doghead Mountain
#                         DuBois locates N16 "2 miles southwest of Ozone
#                         Station" and N9-N11 "just south of Ozone Station". The
#                         source compilation's note reads this as 2 miles
#                         southwest of N1m instead, which appears to be a
#                         misreading and would misplace the site.
#
# A third candidate, merging Pesonen R49 into Kama Hill (R47 + Robertson &
# Fahrig S14), is REFUTED and deliberately absent. R47 and R49 carry identical
# coordinates and identical notes in the source compilation, but the Geological
# Survey of Finland report version of Pesonen (1979) (Report 1942,
# Q 20/27.2/79/1) states that the 270 hand samples came "from 18 reversed and 22
# normal dikes and from 10 reversed sills", and its Figure 1 plots R47 and R49
# as two separate boxed sill sites on the north shore of Nipigon Bay, R49 lying
# slightly north of R47. Sites R40-R49 are therefore ten distinct sills by the
# author's own count, and the shared coordinates are a digitizing artifact that
# should be corrected in the source compilation.
CANDIDATE_MERGES = {
    'S8': 'N19-21',
    'N16-18': 'Doghead Mountain',
}

# Excluded from the pole but retained in the contribution (result_quality 'b').
PILLAR_LAKE = {'PL_ST', 'PL_SAF', 'PL_BHT', 'PL_BHAF'}
NORMAL_POLARITY = {'M6', 'M13', 'M15'}

# MagIC method codes per source study. The codes carried in the source
# compilation mix lab-treatment codes with analysis codes and do not distinguish
# blanket cleaning from stepwise demagnetization, so they are reassigned here
# from the methods described in each paper:
#   DuBois (1962)      blanket AF cleaning (~20 mT) or NRM; no stepwise demag,
#                      no PCA.
#   Robertson & Fahrig blanket AF cleaning at 40 mT with thermal treatment to
#     (1971)           550 C used to demonstrate stability.
#   Pesonen (1979)     stepwise AF in 50-100 Oe steps to 1000 Oe with no blanket
#                      cleaning, plus thermal cleaning.
#   Middleton et al.   stepwise thermal and/or low-temperature cycling with PCA
#     (2004)           best-fit lines.
#   Borradaile &       low-temperature cycling followed by incremental thermal
#     Middleton (2006) (6-12 steps) and/or AF (12-20 steps) with PCA.
METHOD_CODES = {
    'DuBois, 1962': 'DE-BLANKET:LP-DIR-AF-BLANKET',
    'Robertson and Fahrig, 1971': 'DE-BLANKET:LP-DIR-AF-BLANKET:LP-DIR-T',
    'Pesonen, 1979': 'LP-DIR-AF:LP-DIR-T:DE-BFL',
    'Middleton et al. 2004': 'LP-DIR-T:LP-DIR-LT:DE-BFL',
    'Borradaile & Middleton 2006': 'LP-DIR-AF:LP-DIR-T:LP-DIR-LT:DE-BFL',
    'Piispa et al. manuscript in preparation': 'LP-DIR-AF:LP-DIR-T:DE-BFL',
}

AREA_TO_REGION = {
    'Nipigon': 'Nipigon Embayment',
    'Osler': 'Nipigon Embayment',
    'Logan': 'Thunder Bay-Logan',
    'Thunder Bay': 'Thunder Bay-Logan',
    'Thunder Bay?': 'Thunder Bay-Logan',
    'Pigeon': 'Thunder Bay-Logan',
    'Sibley': 'Thunder Bay-Logan',
    'Sibley tip': 'Thunder Bay-Logan',
}

SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'dir_comp_name', 'dir_dec', 'dir_inc', 'dir_polarity',
    'dir_k', 'dir_alpha95', 'dir_n_samples',
    'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description',
]
LOC_COLS = [
    'location', 'location_type', 'result_name', 'result_type', 'result_quality',
    'method_codes', 'citations', 'geologic_classes', 'lithologies',
    'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k',
    'pole_n_sites', 'sites', 'description',
]


def normalize_sill_name(name):
    """Normalizes the cross-study sill names of the source compilation.

    Strips the trailing query marks and the appended geochronology note so that
    e.g. "Terry Fox?" and "Terry Fox", or "Mt Mackay 1106.3 B20" and
    "Mt Mackay", resolve to the same unit.

    Args:
        name (str or float): Raw value of the "Sill/dyke name for dublicates"
            column; NaN where the compilation gives no cross-study correlation.

    Returns:
        str or None: The normalized sill name, or None if no name was given.
    """
    if not isinstance(name, str) or not name.strip():
        return None
    name = name.strip().rstrip('?').strip()
    for suffix in (' 1106.3 B20',):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name.strip()


def assign_units(df):
    """Assigns each accepted published result to an independent cooling unit.

    Applies Stage A (collapse of the repeat demagnetization treatments of
    Middleton et al., 2004 to their sampling localities) followed by Stage B
    (merge of sills sampled by more than one study, using the cross-study sill
    correlations of the source compilation). A cross-study sill name given for
    any member of a Middleton locality is propagated to the whole locality, so
    that e.g. Middleton's Hele 1 (sites 4 and 5) joins Pesonen's R44 as the
    single unit "Steward Lake".

    Args:
        df (pd.DataFrame): Accepted rows of the source compilation.

    Returns:
        tuple[pd.Series, pd.Series]: ``(provisional, unit)``, both indexed like
        ``df``. ``provisional`` is the Stage A grouping (the sampling locality,
        so that repeat demagnetization treatments of the same specimens share a
        value); ``unit`` is the final independent cooling unit after Stage B.
    """
    # Stage A: Middleton site codes -> sampling locality; everything else keeps
    # its own site code as the provisional unit.
    provisional = df.apply(
        lambda r: MIDDLETON_LOCALITY.get(r['site'], r['site'])
        if r['Reference'] == 'Middleton et al. 2004' else r['site'],
        axis=1)

    sill = df['Sill/dyke name for dublicates'].map(normalize_sill_name)

    # Propagate a cross-study sill name to every row of the provisional unit.
    propagated = {}
    for prov, group in sill.groupby(provisional):
        named = group.dropna().unique()
        if len(named) > 1:
            raise ValueError(f'conflicting sill names for {prov}: {named}')
        if len(named) == 1:
            propagated[prov] = named[0]

    # Stage B: the unit is the cross-study sill name where one exists, plus the
    # same-sill correlations established here from the source publications.
    unit = provisional.map(lambda p: propagated.get(p, p))
    unit = unit.map(lambda u: ADDITIONAL_MERGES.get(u, u))
    return provisional, unit.map(lambda u: UNIT_RENAME.get(u, u))


def load_accepted():
    """Loads the source compilation and flags the rows excluded from the pole.

    Returns:
        pd.DataFrame: All 62 source rows with added ``unit``, ``region``,
        ``exclusion`` and ``accepted`` columns. ``exclusion`` is an empty
        string for accepted rows and otherwise states why the row is excluded.
    """
    df = pd.read_csv(SOURCE)

    excluded = pd.Series('', index=df.index)
    excluded[df['site'].isin(PILLAR_LAKE)] = (
        'Pillar Lake Lava; age relative to the sills unconstrained')
    excluded[df['site'].isin(NORMAL_POLARITY)] = (
        'normal-polarity counterpart of a reverse determination of the same '
        'Middleton et al. (2004) locality; interpreted as an incompletely '
        'removed overprint')
    excluded[df['Notes'].fillna('').str.contains('Excluded')] = (
        'flagged "Excluded!" in the source compilation')
    df['exclusion'] = excluded
    df['accepted'] = excluded == ''

    # Middleton site 11 (Rift): no reported coordinates; use the study mean and
    # compute the VGP the source compilation could not.
    missing = df['lat'].isna() & (df['Reference'] == 'Middleton et al. 2004')
    df.loc[missing, ['lat', 'lon']] = MIDDLETON_MEAN_LATLON
    for i in df.index[missing & df['dir_dec'].notna()]:
        r = df.loc[i]
        plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'],
                                          r['dir_alpha95'], r['lat'], r['lon'])
        if plat < 0:                       # report reverse VGPs as north poles
            plat, plon = -plat, (plon + 180) % 360
        df.loc[i, ['vgp_lat', 'vgp_lon', 'vdp_dp', 'vgp_dm']] = plat, plon % 360, dp, dm

    df['locality'], df['unit'] = assign_units(df)
    # not applied; carried so the notebook can report the sensitivity
    df['candidate_merge'] = df['unit'].map(CANDIDATE_MERGES).fillna('')
    df['region'] = df['Area'].map(AREA_TO_REGION)
    return df


def unit_means(accepted):
    """Averages the accepted published results into independent unit means.

    Contributing VGPs are combined with a Fisher mean and, separately, the
    contributing directions are combined with a Fisher mean. Where a unit has a
    single contributing result its published direction, k and alpha95 are
    carried through unchanged; where results are merged, k and alpha95 describe
    the dispersion among the merged determinations rather than within a site.

    Args:
        accepted (pd.DataFrame): Accepted rows with a ``unit`` column.

    Returns:
        pd.DataFrame: One row per unit, sorted by region then unit name.
    """
    rows = []
    for unit, g in accepted.groupby('unit'):
        g = g.sort_values('site')
        vgp = ipmag.fisher_mean(dec=list(g['vgp_lon']), inc=list(g['vgp_lat']))
        direction = ipmag.fisher_mean(dec=list(g['dir_dec']), inc=list(g['dir_inc']))
        single = len(g) == 1
        rows.append({
            'unit': unit,
            'region': g['region'].mode().iat[0],
            'area': ':'.join(sorted(set(g['Area']))),
            'lat': g['lat'].mean(),
            'lon': g['lon'].mean(),
            'dir_dec': direction['dec'],
            'dir_inc': direction['inc'],
            'dir_k': g['dir_k'].iat[0] if single else direction['k'],
            'dir_alpha95': g['dir_alpha95'].iat[0] if single else direction['alpha95'],
            # Repeat demagnetization treatments of a Middleton et al. (2004)
            # locality re-measure the same specimens, so samples are counted
            # once per locality (the largest treatment) rather than summed.
            'dir_n_samples': int(g.groupby('locality')['dir_n_samples'].max().sum()),
            'n_results': len(g),
            'vgp_lat': vgp['inc'],
            'vgp_lon': vgp['dec'],
            'citations': ':'.join(sorted({CITATION[r] for r in g['Reference']})),
            'lithologies': 'Diabase' if 'Diabase' in set(g['lithologies']) else
                           sorted(set(g['lithologies']))[0],
            'geologic_types': 'Sill' if 'Sill' in set(g['geologic_types']) else
                              sorted(set(g['geologic_types']))[0],
            'method_codes': ':'.join(sorted({c for ref in g['Reference']
                                             for c in METHOD_CODES[ref].split(':')})),
            'source_sites': ':'.join(g['site']),
            'source_refs': '; '.join(sorted(set(g['Reference']))),
        })
    units = pd.DataFrame(rows)

    # dp/dm of the unit mean VGP, from the unit mean direction at the unit
    # mean location.
    dp, dm = [], []
    for _, r in units.iterrows():
        _, _, r_dp, r_dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'],
                                        r['dir_alpha95'], r['lat'], r['lon'])
        dp.append(r_dp)
        dm.append(r_dm)
    units['vgp_dp'] = dp
    units['vgp_dm'] = dm
    return units.sort_values(['region', 'unit']).reset_index(drop=True)


def main():
    df = load_accepted()
    accepted = df[df['accepted']].copy()
    units = unit_means(accepted)

    # ---- accepted_results.csv: the published site-level results behind the units
    keep = ['site', 'Reference', 'unit', 'candidate_merge', 'locality', 'region', 'Area', 'lat', 'lon',
            'dir_dec', 'dir_inc', 'dir_alpha95', 'dir_k', 'dir_n_samples',
            'Pol', 'vgp_lat', 'vgp_lon', 'result_quality', 'Notes',
            'accepted', 'exclusion']
    df[keep].to_csv(HERE / 'accepted_results.csv', index=False)
    print(f'Wrote {HERE / "accepted_results.csv"}: {len(df)} source rows '
          f'({int(df["accepted"].sum())} accepted, '
          f'{int((~df["accepted"]).sum())} excluded)')

    # ---- sites.txt: one row per independent unit
    rows = []
    for _, u in units.iterrows():
        merged = (f'{u["n_results"]} determinations merged '
                  f'({u["source_sites"]}; {u["source_refs"]})'
                  if u['n_results'] > 1 else
                  f'{u["source_sites"]} ({u["source_refs"]})')
        rows.append({
            'site': u['unit'], 'location': u['region'], 'result_type': 'i',
            'result_quality': 'g', 'method_codes': u['method_codes'],
            'citations': u['citations'], 'geologic_classes': 'Igneous',
            'geologic_types': u['geologic_types'], 'lithologies': u['lithologies'],
            'lat': f'{u["lat"]:.4f}', 'lon': f'{u["lon"]:.4f}',
            'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '0', 'dir_comp_name': 'ChRM',
            'dir_dec': f'{u["dir_dec"]:.1f}', 'dir_inc': f'{u["dir_inc"]:.1f}',
            'dir_polarity': 'r',
            'dir_k': f'{u["dir_k"]:.1f}', 'dir_alpha95': f'{u["dir_alpha95"]:.1f}',
            'dir_n_samples': str(u['dir_n_samples']),
            'vgp_lat': f'{u["vgp_lat"]:.1f}', 'vgp_lon': f'{u["vgp_lon"]:.1f}',
            'vgp_dp': f'{u["vgp_dp"]:.1f}', 'vgp_dm': f'{u["vgp_dm"]:.1f}',
            'description': f'{u["area"]}; {merged}',
        })

    # excluded source results are retained with result_quality 'b'
    for _, r in df[~df['accepted']].iterrows():
        has_dir = not np.isnan(r['dir_dec'])
        rows.append({
            'site': r['site'], 'location': r['region'], 'result_type': 'i',
            'result_quality': 'b', 'method_codes': METHOD_CODES[r['Reference']],
            'citations': CITATION[r['Reference']], 'geologic_classes': 'Igneous',
            'geologic_types': r['geologic_types'], 'lithologies': r['lithologies'],
            'lat': f'{r["lat"]:.4f}' if not np.isnan(r['lat']) else '',
            'lon': f'{r["lon"]:.4f}' if not np.isnan(r['lon']) else '',
            'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '0', 'dir_comp_name': 'ChRM',
            'dir_dec': f'{r["dir_dec"]:.1f}' if has_dir else '',
            'dir_inc': f'{r["dir_inc"]:.1f}' if has_dir else '',
            'dir_polarity': 'r' if r['Pol'] == 'R' else 'n' if r['Pol'] == 'N' else '',
            'dir_k': f'{r["dir_k"]:.1f}' if has_dir else '',
            'dir_alpha95': f'{r["dir_alpha95"]:.1f}' if has_dir else '',
            'dir_n_samples': str(int(r['dir_n_samples'])),
            'vgp_lat': f'{r["vgp_lat"]:.1f}' if has_dir else '',
            'vgp_lon': f'{r["vgp_lon"]:.1f}' if has_dir else '',
            'vgp_dp': f'{r["vdp_dp"]:.1f}' if has_dir else '',
            'vgp_dm': f'{r["vgp_dm"]:.1f}' if has_dir else '',
            'description': f'EXCLUDED FROM POLE: {r["exclusion"]} '
                           f'({r["site"]}; {r["Reference"]})',
        })

    sites_path = OUT / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in SITE_COLS) + '\n')
    n_good = sum(r['result_quality'] == 'g' for r in rows)
    print(f'Wrote {sites_path}: {n_good} unit means entering the pole '
          f'+ {len(rows) - n_good} excluded source results')

    # ---- locations.txt: the pole and the two regional sub-poles
    def fisher_pole(sub):
        block = ipmag.make_di_block(list(sub['vgp_lon']), list(sub['vgp_lat']))
        return ipmag.fisher_mean(di_block=pmag.flip(block, combine=True))

    loc_rows = []
    groups = [('Nipigon sills', units,
               'Nipigon (Logan) sills ca. 1108 Ma pole')]
    for region in ['Nipigon Embayment', 'Thunder Bay-Logan']:
        groups.append((region, units[units['region'] == region],
                       f'{region} sills ca. 1108 Ma pole'))
    for name, sub, result_name in groups:
        pole = fisher_pole(sub)
        loc_rows.append({
            'location': name, 'location_type': 'Region',
            'result_name': result_name, 'result_type': 'a', 'result_quality': 'g',
            'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DA-DIR-GEO:DE-VGP',
            'citations': ':'.join(sorted({c for cits in sub['citations']
                                          for c in cits.split(':')})),
            'geologic_classes': 'Igneous', 'lithologies': 'Diabase',
            'lat_s': f'{sub["lat"].min():.4f}', 'lat_n': f'{sub["lat"].max():.4f}',
            'lon_w': f'{sub["lon"].min():.4f}', 'lon_e': f'{sub["lon"].max():.4f}',
            'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
            'dir_tilt_correction': '0',
            'pole_lat': f'{pole["inc"]:.1f}', 'pole_lon': f'{pole["dec"]:.1f}',
            'pole_alpha95': f'{pole["alpha95"]:.1f}', 'pole_k': f'{pole["k"]:.1f}',
            'pole_n_sites': str(pole['n']),
            'sites': ':'.join(sub['unit']),
            'description': (
                f'Fisher mean of the {pole["n"]} independent unit VGPs. '
                'Geographic coordinates (intrusive, no tilt correction). '
                'Compiled from DuBois (1962), Robertson & Fahrig (1971), '
                'Pesonen (1979), Middleton et al. (2004), Borradaile & '
                'Middleton (2006) and Piispa et al. (in prep.) by '
                'E. J. Iloranta (2026).'),
        })

    loc_path = OUT / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(LOC_COLS) + '\n')
        for r in loc_rows:
            f.write('\t'.join(str(r[c]) for c in LOC_COLS) + '\n')
    for r in loc_rows:
        print(f'  {r["location"]:20s} pole {r["pole_lat"]}N {r["pole_lon"]}E '
              f'A95={r["pole_alpha95"]} K={r["pole_k"]} N={r["pole_n_sites"]}')
    print(f'Wrote {loc_path}')


if __name__ == '__main__':
    main()
