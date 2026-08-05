"""Combine per-notebook Nordic summary CSVs into a single table.

Each pole notebook writes a single-row CSV (named after the notebook) into this
folder via ``pole_tools.save_nordic_summary``, using the exact Nordic Workshop
compilation columns (``pole_tools.NORDIC_COLUMNS``). This script concatenates
those rows into ``nordic_summaries_combined.csv`` so the result can be pasted
directly into the Nordic format.

It is done with the ``csv`` module rather than pandas so the exact header is
preserved — the Nordic format intentionally repeats some column labels
(a second ``f``/``INCf``/... block and a duplicate ``ROCKNAME``), which pandas
would rename. Rows are sorted by ``nominal age``.

After concatenation, any **empty** cell in a summary row is back-filled from the
matching row of the existing compilation (``data/Laurentia_poles.csv``), matched
by ``ROCKNAME``. This carries over values the per-notebook summaries leave blank
(e.g. the legacy Van der Voo ``Q`` criteria) without ever overwriting a value
the notebook computed — recreated quantities (pole position, A95, R-scores)
always win; only blanks are filled.

The combined table is then rendered to ``pole_table.tex``, the LaTeX
``longtable`` of the compilation used in the manuscript. Pole references are
emitted as ``\\citet`` keys in the manuscript's ``SurnameYYYYa`` convention.
Each pole credits every study that contributed sites to it, taken from
:data:`POLE_REFERENCES`; a pole not listed there is keyed from its own
``POLE AUTHORS`` and ``YEAR`` by :func:`citekey`. The keys are written out
rather than looked up, so the table renders identically wherever it is built;
whether every key resolves is settled by the manuscript's BibTeX file at
compile time.

Usage:
    python combine_nordic_summaries.py
    python combine_nordic_summaries.py --no-tex
"""

import argparse
import csv
import glob
import math
import os
import re
import unicodedata
from collections import defaultdict

SUMMARY_DIR = os.path.dirname(os.path.abspath(__file__))
COMBINED_FILENAME = 'nordic_summaries_combined.csv'
# the existing compilation (same Nordic columns) used to back-fill blank cells
COMPILATION_PATH = os.path.join(SUMMARY_DIR, os.pardir, 'Laurentia_poles.csv')
# Columns never back-filled from the compilation: the recreated poles are
# VGP-Fisher-mean poles with a circular A95, so their oval semi-axes DP/DM (and
# the flattening-block DPf/DMf) are intentionally blank and must stay blank
# rather than inherit the compilation's oval values.
NO_BACKFILL_COLUMNS = frozenset({'DP', 'DM', 'DPf', 'DMf'})

# Nordic grades kept in the compilation. Grade A and B poles are compiled;
# Grade C marks a pole assessed and judged not reliable enough to guide
# reconstruction. A Grade C pole, like one named in EXCLUDED_POLES below, keeps
# its notebook as the record of the assessment but is dropped from the combined
# table and from everything built on it -- the compilation page, the pole map,
# the paleolatitude figures, and the manuscript pole table.
COMPILED_GRADES = frozenset({'A', 'B'})

# Poles excluded from the compilation for reasons other than their grade,
# ROCKNAME -> reason. Their notebooks are retained as the record of the work.
EXCLUDED_POLES = {
    'Dubawnt Group': (
        'out of scope: the compilation covers Laurentia from ca. 1800 Ma, and '
        'the Baker Lake Group strata sampled by Park et al. (1973) are ca. '
        '1830-1790 Ma (Rainbird et al., 2003, 2006). The ca. 1756 Ma Wharton '
        'Group pole of Raub et al. (2026) carries the upper Dubawnt '
        'Supergroup in the compilation instead. There is a second, '
        'independent reason for caution: the Baker Lake pole sits only 1.6 '
        'deg from the Wharton pole despite an age gap of ~40 Myr, and up to '
        '77 Myr against the older bound of its bracket. Either apparent polar '
        'wander was near-stationary across that interval, or the Baker Lake '
        'remanence was reset during Wharton time -- a possibility the Park et '
        'al. (1973) data carry no field test to exclude. Raub et al. (2026) '
        'describe the two as complementary rather than raising '
        'remagnetization, but the coincidence is worrisome'),
}

# --- LaTeX pole table -------------------------------------------------------
TEX_FILENAME = 'pole_table.tex'
# Locality whose implied paleolatitude is tabulated (Duluth, Minnesota),
# longitude in 0-360 deg east. Matches scripts/build_pole_map.py.
DULUTH_LAT = 46.7867
DULUTH_LON = 267.8995
# Terranes that have rifted from Laurentia are rotated back into the Laurentia
# reference frame before the Duluth paleolatitude is computed, so that column is
# comparable across the whole table. Values mirror pole_tools.TERRANE_EULER_POLES
# (Greenland: Roest & Srivastava, 1989; Scotland: Torsvik & Cocks, 2017), and are
# [pole latitude, pole longitude, rotation angle] in degrees. The Plat/Plon
# columns remain present-day positions; only the paleolatitude is rotated.
TERRANE_EULER_POLES = {
    'Laurentia-Greenland': [67.5, -118.5, -13.8],
    'Laurentia-Greenland-Nain': [67.5, -118.5, -13.8],
    'Laurentia-Scotland': [78.6, 161.9, -31.0],
    'Laurentia-Svalbard': [-81.0, 125.0, 68.0],
}
# Subtle ROCKNAME changes between the per-notebook summaries and the legacy
# compilation (Laurentia_poles.csv): spelling, abbreviation, a "MEAN " prefix,
# Dike/Dyke, or a reordered/expanded label. Each entry maps the summary ROCKNAME
# to its compilation ROCKNAME so blank cells still back-fill. Every mapping below
# was confirmed to be the SAME pole (pole position within a few degrees). NOT
# included (genuinely different / superseding / new poles, even where the name is
# similar): "Cardenas Basalts" vs "Cardenas Basalts and Intrusions" (lavas-only
# recompilation, ~16 deg away); "Central Arizona diabases" vs "... -N" (~15 deg
# away); "NW Ontario Lamprophyre Dykes" vs "... and Abitibi Dykes" (compilation
# pole also includes the Abitibi dykes); "Chuar Group (combined)" (new combined
# pole; legacy was "Kwagunt Formation"); "Lower Freda Formation" / "Upper Freda
# Formation" (new split of the legacy "Freda Sandstone").
ROCKNAME_ALIASES = {
    'Giant Gabbro Dikes': 'Giant Gabbro Dykes',
    'Hviddal': 'Hviddal Giant Dyke',
    'McNamara': 'McNamara Formation',
    'Midsommersoe Dolerites': 'Midsommersoe Dolerite',
    'Mean Rocky Mountain intrusions': 'MEAN Rocky Mountain intrusions',
    'Mistastin Batholith': 'Mistastin Pluton',
    'Narsaqq': 'Narssaq Gabbro',
    'Nipigon sills and lavas': 'MEAN Nipigon sills and lavas',
    # This compilation uses the modern lithostratigraphic name; the legacy
    # compilation and the GPMDB carry the unit as the "Nonesuch Shale".
    'Nonesuch Formation': 'Nonesuch Shale',
    'North Qoroq Intrusion': 'North Qoroq Intr.',
    'South Qoroq Intrusion': 'South Qoroq Intr.',
    'Snowslip': 'Snowslip Formation',
    'Spokane': 'Spokane Formation',          # same formation; recomputed pole ~8 deg off the legacy one
    'St. Francois Mountains Acidic Rocks': 'St.Francois Mountains Acidic Rocks',
    'Stoer Group': 'MEAN Stoer Group',
    'Torridon Group': 'MEAN Torridon Group',
    'Western Channel diabase': 'Western Channel Diabase',
    'Pilcher, Garnet Range, Libby': 'Pilcher, Garnet Range and Libby Formations',
    'North Shore Volcanic Group -N (combined)': 'North Shore lavas -N',
    'Sudbury Dike Swarm': 'Sudbury Dykes Combined',
    'Osler Volcanic Group reverse lower': 'Lower Osler volcanics -R',
    'Osler Volcanic Group reverse middle': 'Middle Osler volcanics -R',
    'Osler Volcanic Group reverse upper': 'Upper Osler volcanics -R',
}

# --- pole reference -> BibTeX key(s) ----------------------------------------
# The reference column credits every study that contributed sites to a pole,
# not only the paper the pole is named for: these are recalculated poles, and
# the site-level data behind them frequently come from several studies (the
# Nipigon sills pole, for instance, pools sites from five). The lists below were
# derived from the ``citations`` column of the MagIC site files each notebook
# reads -- ``sites.txt`` for most poles, but the study-specific files for the
# Sudbury, Mackenzie and Mistastin poles, which their notebooks combine -- and
# restricted to the sites each notebook actually uses (its quality and
# polarity-zone filters), then resolved to keys against the manuscript bib by
# DOI. Five references cited by the site data were missing from the bib and were
# added with it: Larochelle (1967), Murthy et al. (1968), Evans et al. (1975),
# Stupavsky & Symons (1982) and Fahrig (1986).
#
# Keys follow the manuscript's SurnameYYYY[a-z] convention. A ROCKNAME absent
# from this mapping falls back to citekey() on its own POLE AUTHORS/YEAR, which
# gives the "a" suffix and so needs checking against the bib; the table build
# reports any entry here that matches no row.
POLE_REFERENCES = {
    'Sept-Iles Layered Intrusion': ['Tanczyk1987a'],
    'Catoctin Basalts': ['Meert1994a'],
    'Callander Alkaline Complex': ['Symons1991a'],
    'Baie des Moutons complex A': ['McCausland2011a'],
    'Baie des Moutons complex B': ['McCausland2011a'],
    'Long Range Dykes': ['Murthy1992a'],
    'Franklin event grand mean': ['Denyszyn2009b'],
    'Chuar Group (combined)': ['Weil2004a', 'Eyster2020a'],
    'Uinta Mountain Group': ['Weil2006b'],
    'Gunbarrel LIP': ['Harlan1997a', 'Harlan2003a', 'Ding2025a'],
    'Adirondack metamorphic anorthosite': ['Brown2012a'],
    'Torridon Group': ['Evans2021a'],
    'Jacobsville Formation': ['Zhang2024a'],
    'Upper Freda Formation': ['Fuentes2025a'],
    'Lower Freda Formation': ['Henry1977a'],
    'Nonesuch Formation': ['Slotznick2024a'],
    'Cardenas Basalts': ['Zhang2024b'],
    'Michipicoten Island Formation': ['Palmer1987a', 'Fairchild2017a'],
    'Lake Shore Traps': ['Diehl1994a', 'Kulakov2013a'],
    'Schroeder Lutsen Basalts': ['Tauxe2009a', 'Fairchild2017a'],
    'Portage Lake Volcanics': [
        'Books1972a', 'Hnat2006a', 'Foucher2018a', 'Swanson-Hysell2019a'
    ],
    'Uppermost Mamainse Point volcanics -N': ['Swanson-Hysell2014a'],
    'North Shore Volcanic Group -N (combined)': [
        'Books1972a', 'Tauxe2009a', 'Swanson-Hysell2019a'
    ],
    'Central Arizona diabases': ['Harlan1993a', 'Donadini2011b'],
    'Mamainse Point volcanics -C (lower N, upper R)': ['Swanson-Hysell2014a'],
    'Lower Mamainse Point volcanics -R2': ['Swanson-Hysell2014a'],
    'Osler Volcanic Group reverse upper': [
        'Halls1974a', 'Swanson-Hysell2014b', 'Swanson-Hysell2019a'
    ],
    'Coldwell Complex': ['Kulakov2014a'],
    'Osler Volcanic Group reverse middle': ['Swanson-Hysell2014b'],
    'Nipigon sills': [
        'Dubois1962a', 'Robertson1971a', 'Pesonen1979a',
        'Middleton2004a', 'Borradaile2006a'
    ],
    'Osler Volcanic Group reverse lower': ['Swanson-Hysell2014b'],
    'Lowermost Mamainse Point volcanics -R1': ['Swanson-Hysell2014a'],
    'Abitibi Dykes': ['Ernst1993a', 'Halls2005a'],
    'NW Ontario Lamprophyre Dykes': ['Ernst1993a', 'Piispa2018a'],
    'NE-SW Trending Dyke Swarm': ['Piper1992a'],
    'Giant Gabbro Dikes': ['Piper1977a'],
    'South Qoroq Intrusion': ['Piper1992a'],
    'Hviddal': ['Piper1977a'],
    'Narsaqq': ['Piper1977a'],
    'Stoer Group': ['Evans2021a'],
    'Sudbury Dike Swarm': [
        'Larochelle1967a', 'Palmer1977a', 'Stupavsky1982b'
    ],
    'Mackenzie dykes recompiled': [
        'Fahrig1969a', 'Robertson1969a', 'Irving1972c', 'Park1974a',
        'Fahrig1986a'
    ],
    'West Gardar Dolerite Dykes': ['Piper1977b'],
    'West Gardar Lamprophyre Dykes': ['Piper1977b'],
    'Kungnat Ring Dyke': ['Piper1977b'],
    'North Qoroq Intrusion': ['Piper1992a'],
    'Nain Anorthosite': ['Murthy1978a'],
    'Midsommersoe Dolerites': ['Marcussen1983a'],
    'Victoria Fjord dolerite dykes': ['Abrahamsen1987a'],
    'Zig-Zag Dal Basalts': ['Marcussen1983a'],
    'Pilcher, Garnet Range, Libby': ['Elston2002a'],
    'McNamara': ['Elston2002a'],
    'Purcell Lava': ['Evans1975a', 'Elston2002a'],
    'Mean Rocky Mountain intrusions': ['Harlan1994a', 'Harlan1998a'],
    'Mistastin Batholith': ['Fahrig1976a', 'Herve2015a'],
    'Snowslip': ['Elston2002a'],
    'Spokane': ['Elston2002a'],
    'St. Francois Mountains Acidic Rocks': ['Meert2002b', 'Bray2021a'],
    'Michikamau Intrusion Combined': ['Murthy1968a', 'Emslie1976a'],
    'Western Channel diabase': ['Irving1972b'],
    'Melville Bugt diabase dykes': ['Halls2011a'],
    'Cleaver Dykes': ['Harlan2003a', 'Irving2004a', 'Ootes2015a'],
    'Wharton Group': ['Raub2026a'],
    'NE trending ECMB Diabase Dykes': ['Swanson-Hysell2021b'],
}

# --- key-pole flag ----------------------------------------------------------
# The key-pole flag is a SEPARATE axis from the R-criteria columns. The R
# scores follow Meert et al. (2020) as closely as possible
# (resources/Meert2020_R_criteria.md); the key-pole flag follows Buchan (2013)
# as closely as possible. The two use different vocabularies and different
# thresholds, so the flag is NOT derived from R1/R2/R4/R5 -- it is a curated
# assignment, recorded below one pole at a time.
#
# Buchan's criteria (Buchan, 2013, section 2, and the scope statement in
# section 6) are:
#
#   B1  The age of the paleopole is precisely determined -- U-Pb, or
#       occasionally Ar-Ar; typically within +/-10 Myr, +/-20 Myr at most.
#       DEVIATION: we hold to +/-15 Myr, the tighter standard of Meert R1,
#       rather than Buchan's +/-20 Myr ceiling. Buchan (2013, section 2)
#       anticipates this, noting that "it should in future be possible to
#       significantly tighten this latter uncertainty". This is the one place
#       where the compilation departs from him, and it is deliberate.
#   B2  The paleopole is of good quality -- the primary remanence is properly
#       isolated by stepwise AF or thermal demagnetization, and secular
#       variation is largely averaged out. This is a qualitative judgment and
#       is deliberately weaker than the statistical gates of Meert R2.
#   B3  A positive field test establishes the remanence as primary. Buchan
#       admits exactly the seven tests of BUCHAN_TESTS below; the tests of his
#       section 4 (partial baked contact, antipodal reversals, reversals in a
#       single stratigraphic section, tectonic fold test, conglomerate test on
#       a younger conglomerate) do NOT establish a remanence as primary.
#   B4  The pole is from the interior of a large craton. Buchan excludes
#       cratonic margins and microcontinents from his analysis outright rather
#       than scoring them, so this acts as a scope filter.
#
# Buchan's field-test letters, from the notes to his Table 1. Note that these
# collide with resources/field_test_codes.md ('c' is an intraformational
# conglomerate test here but an inverse baked-contact test there), which is a
# further reason to keep the two systems apart.
BUCHAN_TESTS = {
    'b': 'baked contact test',
    'b(p)': 'baked contact profile test',
    'c': 'intraformational conglomerate test',
    'p': 'polarity correlation test',
    's': 'secular variation correlation test',
    'd': 'remanence direction correlation test',
    'x': 'regional consistency test',
    'ef': 'fold test contemporaneous with emplacement',
}

# Key-pole assignments, ROCKNAME -> (is_key, test, basis).
#
#   is_key  True or False.
#   test    Buchan test letter(s) from BUCHAN_TESTS supporting a True call;
#           '' when is_key is False.
#   basis   Where the call comes from. 'Buchan (2013) Table 1' / 'Table 2'
#           means his own published call, adopted verbatim. Anything else is
#           our assignment applying his criteria, and should say which
#           criterion decided it.
#
# A pole with no entry here is UNASSESSED: it renders as '--' in the table and
# is reported by the script until a call is recorded.
KEY_POLES = {
    # -- Buchan (2013) Table 1, Laurentia and northeastern Greenland blocks,
    # adopted verbatim. Where our compilation recalculates the pole from site
    # level, his judgment on the rock unit and its field test carries over.
    'Long Range Dykes': (
        True, 'b', 'Buchan (2013) Table 1. NOTE: flagged to stay faithful to '
                   'his call, but this is the weakest key pole in the '
                   'compilation and warrants discussion in the text -- his '
                   'pole rests on 5 dykes with A95 18 deg, and the pole '
                   'recreated here has A95 20.8 deg, so secular variation is '
                   'unlikely to be averaged. It also sets the effective floor '
                   'on how permissive his quality criterion is'),
    'Franklin event grand mean': (True, 'x', 'Buchan (2013) Table 1'),
    'Gunbarrel LIP': (True, 'x', 'Buchan (2013) Table 1'),
    'Lake Shore Traps': (True, 'b', 'Buchan (2013) Table 1'),
    'Portage Lake Volcanics': (True, 's', 'Buchan (2013) Table 1'),
    'Osler Volcanic Group reverse upper': (True, 'ef', 'Buchan (2013) Table 1'),
    'Nipigon sills and lavas': (True, 'b', 'Buchan (2013) Table 1 (Logan sills)'),
    'Abitibi Dykes': (True, 'b', 'Buchan (2013) Table 1'),
    'Sudbury Dike Swarm': (True, 'b(p)', 'Buchan (2013) Table 1'),
    'Mackenzie dykes recompiled': (True, 'b', 'Buchan (2013) Table 1'),
    'Midsommersoe Dolerites': (True, 'b', 'Buchan (2013) Table 1'),
    'Mistastin Batholith': (
        True, 'x', 'Buchan (2013) Table 1 -- regional consistency with the '
                   'ca. 1.46-1.42 Ga suite, not the impact-crater test'),
    'Mean Rocky Mountain intrusions': (
        True, 'x', 'Buchan (2013) Table 1 (Laramie complex and Sherman '
                   'Granite) -- regional consistency'),
    'St. Francois Mountains Acidic Rocks': (
        True, 'ef,x', 'Buchan (2013) Table 1'),
    'Michikamau Intrusion Combined': (True, 'b', 'Buchan (2013) Table 1'),
    'Western Channel diabase': (True, 's,d', 'Buchan (2013) Table 1'),
    'Cleaver Dykes': (True, 'b', 'Buchan (2013) Table 1'),
    # -- Buchan (2013) Table 2: listed as a non-key pole used in his figures.
    'Melville Bugt diabase dykes': (
        False, '', 'Buchan (2013) Table 2 -- listed as non-key, no field test'),

    # -- Midcontinent Rift and coeval magmatism, assessed here against the
    # Buchan criteria. Two arguments do most of the work in this setting.
    #
    # Intraformational conglomerate test (c): at Mamainse Point conglomerate
    # horizons lie between the -R2 and -C zones and between the -C and -N
    # zones. Following Buchan section 3.3, a conglomerate test constrains the
    # unit the clasts were derived from, i.e. the zone immediately below, so
    # the lower conglomerate tests -R2 and the Great Conglomerate tests -C.
    #
    # Regional consistency (x) and remanence direction correlation (d): the
    # rift record is a rapid apparent polar wander track sampled at localities
    # hundreds of kilometres apart. Coeval poles from separate successions
    # agree while poles of different age differ systematically -- a pattern a
    # regional overprint cannot produce. Buchan applies the regional
    # consistency test qualitatively, and it is applied qualitatively here:
    # requiring overlap of A95 circles would perversely penalise the most
    # precisely determined poles, since genuine plate motion separates poles
    # only a few million years apart.
    'Cardenas Basalts': (
        True, 'x', 'regional consistency with the coeval Lake Shore Traps and '
                   'Michipicoten Island poles ~2400 km away'),
    'Michipicoten Island Formation': (
        True, 'x', 'regional consistency with the coeval Cardenas Basalts pole'),
    'Schroeder Lutsen Basalts': (
        True, 'x', 'regional consistency with the coeval Lake Shore Traps and '
                   'Portage Lake Volcanics poles'),
    'Uppermost Mamainse Point volcanics -N': (
        True, 'd,x', 'remanence direction correlation between the Mamainse '
                     'Point polarity zones; regional consistency with the '
                     'coeval Portage Lake Volcanics pole'),
    'North Shore Volcanic Group -N (combined)': (
        True, 'x', 'regional consistency with the coeval Portage Lake '
                   'Volcanics and uppermost Mamainse Point poles'),
    'Mamainse Point volcanics -C (lower N, upper R)': (
        True, 'c', 'Great Conglomerate above the zone carries clasts derived '
                   'from it (Swanson-Hysell et al., 2009)'),
    'Lower Mamainse Point volcanics -R2': (
        True, 'c', 'conglomerate above the zone carries clasts derived from '
                   'it (Swanson-Hysell et al., 2009)'),
    'Coldwell Complex': (
        True, 'x', 'regional consistency with the coeval upper Osler and '
                   'Nipigon (Logan) sill poles'),
    'Osler Volcanic Group reverse middle': (
        True, 'x', 'regional consistency with the coeval Nipigon (Logan) sill '
                   'and lowermost Mamainse Point poles'),
    'Nipigon sills': (
        True, 'b', 'positive baked-contact test: baked Sibley Group and Rove '
                   'Formation adjacent to the sills carry the sill direction '
                   'while unbaked host rock does not (Pesonen, 1979, Table 4). '
                   'This pole is also the coeval anchor for the '
                   'regional-consistency calls on the Coldwell, middle and '
                   'lower Osler, and lowermost Mamainse Point poles'),
    'Osler Volcanic Group reverse lower': (
        True, 'x', 'regional consistency with the coeval Nipigon (Logan) sill '
                   'pole'),
    'Lowermost Mamainse Point volcanics -R1': (
        True, 'd,x', 'remanence direction correlation between the Mamainse '
                     'Point polarity zones; regional consistency with the '
                     'coeval Nipigon (Logan) sill pole'),
    'NW Ontario Lamprophyre Dykes': (
        True, 'b', 'positive baked-contact test on dykes D1 and MM2 (Piispa '
                   'et al., 2018)'),
    'Central Arizona diabases': (
        False, '', 'age uncertain to +/-16 Myr, outside the +/-15 Myr standard '
                   'adopted here'),

    # -- Northeastern Greenland, ca. 1382 Ma. Buchan's key pole "Midsommerso
    # sills, dykes and related volcanics" is a grand mean of exactly these
    # three units -- his N of "10s, 9s, 19s" matches the site counts of the
    # Midsommerso dolerites (10), Victoria Fjord dykes (9) and Zig-Zag Dal
    # basalts (19) compiled separately here, and his footnote j names all
    # three. The baked-contact test earning his designation is the one in the
    # Victoria Fjord dykes (Abrahamsen and Van der Voo, 1987), which he cites
    # as the field-test reference for the grand mean. His call is therefore
    # carried to all three components.
    'Victoria Fjord dolerite dykes': (
        True, 'b', 'component of the Buchan (2013) Table 1 Midsommerso grand '
                   'mean; carries the positive baked-contact test he cites '
                   'for it. NOTE: these dykes have no direct radiometric date '
                   '-- their age rests on the antipodality of their directions '
                   'to the dated Midsommerso dolerites, which Buchan accepted '
                   'in folding them into one 1382 +/- 2 Ma pole'),
    'Zig-Zag Dal Basalts': (
        True, 'b,x', 'component of the Buchan (2013) Table 1 Midsommerso grand '
                     'mean; same 1382 +/- 2 Ma event as the baked-contact-'
                     'tested Victoria Fjord dykes'),

    # -- Gardar province, southern Greenland. Buchan lists none of these
    # despite the Piper (1977, 1992) data long predating his review. Two
    # considerations apply across the province: most Gardar ages trace to
    # unpublished Heaman determinations summarised by Upton (2013), and the
    # Gardar is a continental rift on the southern margin of the Greenland
    # shield, whereas Buchan (section 6) excludes cratonic margins from his
    # analysis outright. Each pole below also fails a specific criterion.
    'NE-SW Trending Dyke Swarm': (
        False, '', 'no field test; age bracketed only by unpublished dates'),
    'Giant Gabbro Dikes': (False, '', 'no field test'),
    'South Qoroq Intrusion': (
        False, '', 'has a baked-contact test, but its age is one of the '
                   'unpublished Heaman determinations and it sits in the '
                   'Gardar rift province. Its A95 of 14.9 deg over 9 sites is '
                   'NOT the basis for exclusion: Buchan flagged the Long '
                   'Range dykes at A95 18 deg on 5 dykes, so his quality '
                   'criterion admits determinations this imprecise'),
    'Hviddal': (False, '', 'no field test; B=7 sites'),
    'Narsaqq': (False, '', 'no field test; B=4 sites, N=22'),
    'West Gardar Dolerite Dykes': (
        False, '', 'age rests on Rb-Sr biotite dates from an unpublished '
                   'study; Buchan admits only U-Pb and occasionally Ar-Ar'),
    'West Gardar Lamprophyre Dykes': (
        False, '', 'age rests on Rb-Sr biotite dates from an unpublished '
                   'study; Buchan admits only U-Pb and occasionally Ar-Ar'),
    'Kungnat Ring Dyke': (
        False, '', 'no field test; B=4 sites with K=455 record a spot reading '
                   'rather than averaged secular variation'),
    'North Qoroq Intrusion': (
        False, '', 'only an inverse baked-contact test; the U-Pb date is on '
                   'the adjacent Motzfeldt complex'),

    # -- Belt-Purcell basin. Every pole rests on a tectonic fold test and
    # nothing else. Belt-Purcell folding is Cretaceous-Paleogene, roughly
    # 1.4 Gyr after deposition, and Buchan section 4.4 admits a fold test only
    # where folding is contemporaneous with deposition or emplacement.
    #
    # The regional consistency test does not substitute. That test is not
    # restricted to a single igneous province -- Buchan section 3.7 extends it
    # to "widely scattered but discrete units of similar age", and he applies
    # it that way in Table 1, granting it to the Harp Lake, Laramie-Sherman
    # and Mistastin complexes, discrete intrusions ~3500 km apart spanning
    # ca. 1476-1420 Ma. What it does require is units in "the interior of a
    # craton with no evidence of subsequent regional metamorphism", and that
    # is where the Belt-Purcell poles fail: these rocks were transported in
    # the Cordilleran fold-and-thrust belt. Four of the five also fail Meert
    # R5, i.e. they are not in structural coherence with the craton.
    'Pilcher, Garnet Range, Libby': (
        False, '', 'tectonic fold test only, folding ca. 1.4 Gyr after '
                   'deposition; age also uncertain to +/-22 Myr'),
    'McNamara': (
        False, '', 'tectonic fold test only, folding ca. 1.4 Gyr after '
                   'deposition'),
    'Purcell Lava': (
        False, '', 'tectonic fold test only, folding ca. 1.4 Gyr after '
                   'deposition, and the determination is too poorly resolved '
                   'to argue primary remanence: the remanence is largely '
                   'hematite-held, with pigmentary and specular hematite '
                   'alongside magnetite and pyrrhotite, and site directions '
                   'are scattered (K=17.9 over B=15 sites). Its pole does lie '
                   '4.0 deg from the coeval Rocky Mountain intrusions pole '
                   '1209 km away, a similarity worth noting, but that '
                   'agreement cannot carry a pole whose primary remanence is '
                   'not established'),
    'Snowslip': (
        False, '', 'tectonic fold test only, folding ca. 1.4 Gyr after '
                   'deposition'),
    'Spokane': (
        False, '', 'tectonic fold test only, folding ca. 1.4 Gyr after '
                   'deposition'),

    # -- Scotland. Both are Lulea Working Group grand means of six or seven
    # studies, and both fail on age, on field test, and on scope.
    'Torridon Group': (
        False, '', 'age uncertain to +/-75 Myr; tectonic fold test only; from '
                   'a detached fragment of the Laurentian margin, outside the '
                   'craton interiors Buchan considers'),
    'Stoer Group': (
        False, '', 'age uncertain to +/-70 Myr; tectonic fold test only; from '
                   'a detached fragment of the Laurentian margin, outside the '
                   'craton interiors Buchan considers'),

    # -- Neoproterozoic. All seven fail, four of them plainly. Note that
    # Buchan's quality criterion cannot be the discriminator anywhere in this
    # block: his Long Range dykes key pole has A95 18 deg on 5 dykes, so site
    # counts and confidence circles of this order do not disqualify.
    'Uinta Mountain Group': (
        False, '', 'no field test; the age rests on a 766.3 Ma detrital-zircon '
                   'maximum against a correlation-based minimum of ca. 730 Ma'),
    'Baie des Moutons complex A': (
        False, '', 'baked-contact tests inconclusive (no stable remanence in '
                   'the host syenite); Buchan section 4.1 holds that an '
                   'incomplete baked contact test yields no information on '
                   'whether a remanence is primary'),
    'Baie des Moutons complex B': (
        False, '', 'baked-contact tests inconclusive, as for complex A; the '
                   'determination is also poorly resolved (B=6, K=10.5)'),
    'Sept-Iles Layered Intrusion': (
        False, '', 'only an inverse baked-contact test, and the direction '
                   'overlaps the Ordovician-Silurian segment of the Laurentia '
                   'path, the expected time of Taconic-Salinic '
                   'remagnetization, which no available test excludes'),
    'Callander Alkaline Complex': (
        False, '', 'has a positive baked-contact test, but the age rests on '
                   'K-Ar biotite dates and an unpublished Pb-Pb determination '
                   'whose value differs between abstract and presentation; '
                   'Buchan admits only U-Pb and occasionally Ar-Ar'),
    'Catoctin Basalts': (
        False, '', 'not flagged, to stay faithful to Buchan, who did not list '
                   'it although Meert et al. (1994) long predates his review. '
                   'It does have a positive baked-contact test and a U-Pb '
                   'age; the likely reason for his omission is that the '
                   'Catoctin lies in the Appalachian Blue Ridge, deformed '
                   'well after magnetization, outside the craton interiors he '
                   'considers'),
    'Chuar Group (combined)': (
        False, '', 'the fold test is not demonstrably syn-depositional, so it '
                   'is a plain tectonic fold test and does not establish a '
                   'primary remanence under Buchan section 4.4'),

    # -- Keweenawan sedimentary rocks. Two considerations run across the block.
    #
    # Field test: the intraformational conglomerate test of Swanson-Hysell,
    # Fairchild and Slotznick (2019b) is on fluvial intraclasts in the lower
    # Freda Formation. Those intraclasts are of the same lithofacies as the
    # fluvial rocks that carry the Nonesuch pole and the lower Freda pole, and
    # as the upper Freda, so the test constrains the remanence of all three.
    #
    # Inclination shallowing: Buchan section 2 requires it to be corrected, and
    # in section 6 he admits sedimentary key poles where shallowing is absent
    # or, for the Lower Akitkan sediments, where "the remanence inclination is
    # very low" so its effect is small. That is the situation here. These units
    # have inclinations of 0.7 to -13.3 deg, so although the flattening factors
    # are substantial (f = 0.43 to 0.73), the poles move only 0.1 to 3.7 deg
    # between the tabulated and inclination-corrected positions. Shallowing is
    # both well quantified and minor in its effect on pole position.
    'Nonesuch Formation': (
        True, 'c', 'intraformational conglomerate test on fluvial intraclasts '
                   'in the conformably overlying lower Freda Formation '
                   '(Swanson-Hysell, Fairchild & Slotznick, 2019b), of the '
                   'same lithofacies as, and in close stratigraphic proximity '
                   'to, the fluvial rocks that carry this pole; inclination '
                   'shallowing is quantified and displaces the pole by only '
                   '3.7 deg'),
    'Lower Freda Formation': (
        True, 'c', 'the intraformational conglomerate test of Swanson-Hysell, '
                   'Fairchild & Slotznick (2019b) is on intraclasts within '
                   'this formation; inclination shallowing displaces the pole '
                   'by 0.1 deg'),
    'Jacobsville Formation': (
        True, 'c', 'two positive intraformational conglomerate tests within '
                   'the formation (Agate Falls, Dover Creek); inclination '
                   'shallowing displaces the pole by 3.5 deg. NOTE: the age '
                   'is equivocal against the +/-15 Myr standard, resting on a '
                   'detrital-zircon maximum depositional age of 992.51 +/- '
                   '0.64 Ma rather than a date on the unit, and the reversal '
                   'test fails'),
    'Upper Freda Formation': (
        False, '', 'the lower Freda intraformational conglomerate test applies '
                   'to this formation as well, but the age is too loose: the '
                   'Bayesian thermal-subsidence model gives +18.2/-13.3 Myr, '
                   'outside the +/-15 Myr standard adopted here'),

    # -- Grenville-age and Paleoproterozoic remainder.
    'Adirondack metamorphic anorthosite': (
        False, '', 'a cooling (exhumation) pole with no field test, and the '
                   'age of magnetization is set by a thermochronologic model '
                   'rather than dated directly (+/-23 Myr)'),
    'Nain Anorthosite': (
        False, '', 'no field test, and the age is uncertain to +/-23 Myr'),
    'Wharton Group': (
        True, 'c', 'robustly positive, definitively intraformational '
                   'conglomerate test on rhyolite clasts within the Pitz '
                   'Formation section at McRae Lake (Raub et al., 2026); the '
                   'pole is dated directly by a 1756.4 +/- 0.9 Ma U-Pb zircon '
                   'age on one of its contributing flows. Raub et al. reach '
                   'the same conclusion independently, writing that the pole '
                   '"constitutes a key pole by the criteria outlined in '
                   'Buchan et al. (2000) and Buchan (2013)"'),
    'NE trending ECMB Diabase Dykes': (
        True, 'x', 'regional consistency with the Wharton Group pole 2027 km '
                   'away and the Cleaver dykes pole 2659 km away, both of '
                   'which are key on a direct field test; three discrete '
                   'igneous units spanning 1779-1740 Ma across the craton '
                   'interior. NOTE: this application of the regional '
                   'consistency test is provisional and is the least settled '
                   'flag in the compilation. Apparent polar wander is slow '
                   'across this interval -- Raub et al. (2026) call it "a '
                   'short arc" -- so agreement between poles tens of Myr '
                   'apart is only weakly diagnostic of primary remanence, '
                   'unlike in the rapidly moving Midcontinent Rift where '
                   'coeval poles a few Myr apart routinely fail to overlap. '
                   'The agreement with Wharton is marginal (15.6 deg against '
                   'a 15.9 deg summed A95); the agreement with Cleaver is '
                   'comfortable (8.5 deg against 11.5 deg). The pole has no '
                   'primary field test of its own, only an inverse baked '
                   'contact against a ca. 1096 Ma dyke. Worth putting to '
                   'Buchan directly before publication'),
}


def _is_empty(value):
    """True for cells that carry no information (blank or a NaN literal)."""
    return value is None or str(value).strip() in ('', 'nan', 'NaN', 'NA')


def _build_column_map(summary_header, compilation_header):
    """Map each summary column index to a compilation column index.

    Columns are matched by (name, occurrence), so the intentionally repeated
    Nordic labels (the second ``f``/``INCf``/... block, the duplicate
    ``ROCKNAME``) line up with the corresponding repeat in the compilation,
    regardless of column order or any extra trailing columns the compilation
    file may carry.

    Returns:
        dict[int, int]: summary column index -> compilation column index.
    """
    comp_occurrences = defaultdict(list)
    for j, name in enumerate(compilation_header):
        comp_occurrences[name].append(j)
    seen = defaultdict(int)
    column_map = {}
    for i, name in enumerate(summary_header):
        k = seen[name]
        seen[name] += 1
        if name in comp_occurrences and k < len(comp_occurrences[name]):
            column_map[i] = comp_occurrences[name][k]
    return column_map


def backfill_from_compilation(header, rows, compilation_path=COMPILATION_PATH):
    """Fill empty cells in ``rows`` from the matching compilation row.

    Rows are matched to the compilation by ``ROCKNAME``. For each matched row,
    every empty summary cell whose corresponding compilation cell has a value is
    filled from the compilation; non-empty summary cells are left untouched.

    Args:
        header (list[str]): The Nordic column header (``NORDIC_COLUMNS``).
        rows (list[list[str]]): Summary rows (mutated in place).
        compilation_path (str): Path to ``Laurentia_poles.csv``.

    Returns:
        tuple[int, list[str]]: (number of cells filled, ROCKNAMEs with no
        compilation match).
    """
    if not os.path.exists(compilation_path):
        print(f'-W- compilation not found at {compilation_path}; '
              'skipping back-fill of blank cells')
        return 0, []
    with open(compilation_path, encoding='utf-8-sig', newline='') as fh:
        records = list(csv.reader(fh))
    if len(records) < 2:
        return 0, []
    comp_header, comp_rows = records[0], records[1:]
    column_map = _build_column_map(header, comp_header)

    rock_idx = header.index('ROCKNAME')
    comp_rock_idx = comp_header.index('ROCKNAME')
    comp_by_rock = {}
    for cr in comp_rows:
        if len(cr) > comp_rock_idx and cr[comp_rock_idx].strip():
            comp_by_rock.setdefault(cr[comp_rock_idx].strip(), cr)

    filled = 0
    unmatched = []
    for row in rows:
        rockname = row[rock_idx].strip() if len(row) > rock_idx else ''
        comp_row = comp_by_rock.get(rockname)
        if comp_row is None and rockname in ROCKNAME_ALIASES:
            comp_row = comp_by_rock.get(ROCKNAME_ALIASES[rockname])
        if comp_row is None:
            unmatched.append(rockname)
            continue
        for i, j in column_map.items():
            if header[i] in NO_BACKFILL_COLUMNS:
                continue
            if (i < len(row) and j < len(comp_row)
                    and _is_empty(row[i]) and not _is_empty(comp_row[j])):
                row[i] = comp_row[j]
                filled += 1
    return filled, unmatched


def combine_summaries(summary_dir=SUMMARY_DIR, combined_filename=COMBINED_FILENAME):
    """Concatenate all per-notebook summary CSVs in a directory into one CSV.

    Args:
        summary_dir (str): Directory holding the per-notebook summary CSVs.
        combined_filename (str): Name of the combined CSV to write into
            ``summary_dir``. Excluded from the inputs if present.

    Returns:
        str: Path to the combined CSV written to disk.
    """
    combined_path = os.path.join(summary_dir, combined_filename)
    # per-pole summary CSVs are named after their notebooks, with a numeric age
    # prefix (e.g. 1086_Lake_Shore_Traps.csv). Only those are combined — this
    # excludes the combined output and any reference file (e.g.
    # Iloranta_Laurentia_preworkshop.csv) kept in this folder.
    csv_paths = sorted(
        p for p in glob.glob(os.path.join(summary_dir, '*.csv'))
        if os.path.basename(p)[:1].isdigit()
        and os.path.abspath(p) != os.path.abspath(combined_path)
    )
    if not csv_paths:
        raise FileNotFoundError(
            f'No per-notebook summary CSVs found in {summary_dir}. '
            'Run the pole notebooks to generate them first.'
        )

    header = None
    rows = []
    for path in csv_paths:
        with open(path, encoding='utf-8-sig', newline='') as fh:
            records = list(csv.reader(fh))
        if len(records) < 2:
            continue
        if header is None:
            header = records[0]
        elif records[0] != header:
            raise ValueError(
                f'Column header of {os.path.basename(path)} does not match the '
                'Nordic columns of the other summaries; re-run that notebook.'
            )
        rows.extend(records[1:])

    # sort by nominal age (column label is unique, so .index is unambiguous)
    age_idx = header.index('nominal age')

    def age_key(row):
        try:
            return float(row[age_idx])
        except (ValueError, IndexError):
            return float('inf')

    rows.sort(key=age_key)

    # drop poles graded outside COMPILED_GRADES (their notebooks are kept as
    # the record of the assessment, but the pole is not compiled)
    grade_idx, rock_idx = header.index('Grade'), header.index('ROCKNAME')

    def _drop_reason(row):
        if row[rock_idx].strip() in EXCLUDED_POLES:
            return EXCLUDED_POLES[row[rock_idx].strip()]
        if row[grade_idx].strip() not in COMPILED_GRADES:
            return f'grade {row[grade_idx].strip()}'
        return None

    excluded = [(r[rock_idx], _drop_reason(r)) for r in rows if _drop_reason(r)]
    rows = [r for r in rows if not _drop_reason(r)]

    # back-fill blank cells from the compilation (e.g. the legacy Q criteria)
    filled, unmatched = backfill_from_compilation(header, rows)

    # Site longitude to 0-360 deg E. Notebooks record the sampling locality in
    # whichever convention their source used, so the summaries arrive with a mix
    # of signed and 0-360 values. The Nordic layout is 0-360 deg E, and the
    # manuscript table normalises on the way out, so without this the same pole
    # reads -88.2 in the compilation CSV and 271.8 in the table.
    slong_idx = header.index('SLONG')
    wrapped = 0
    for row in rows:
        raw = row[slong_idx].strip()
        if not raw:
            continue
        try:
            lon = float(raw)
        except ValueError:
            continue
        if lon < 0 or lon >= 360:
            row[slong_idx] = f'{lon % 360.0:.2f}'.rstrip('0').rstrip('.')
            wrapped += 1

    with open(combined_path, 'w', encoding='utf-8-sig', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    print(f'Combined {len(rows)} summaries into {combined_path}')
    if excluded:
        print(f'Excluded {len(excluded)} pole(s) from the compilation:')
        for name, reason in excluded:
            print(f'    {name}: {reason}')
    if filled:
        print(f'Back-filled {filled} blank cell(s) from '
              f'{os.path.basename(COMPILATION_PATH)}')
    if wrapped:
        print(f'Wrapped {wrapped} site longitude(s) into 0-360 deg E')
    if unmatched:
        print(f'-W- {len(unmatched)} row(s) had no ROCKNAME match in '
              f'{os.path.basename(COMPILATION_PATH)}: '
              f'{", ".join(sorted(set(unmatched)))}')
    return combined_path


# ---------------------------------------------------------------------------
# LaTeX pole table
# ---------------------------------------------------------------------------

_LATEX_SPECIALS = {
    '\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
    '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}',
    '~': r'\textasciitilde{}', '^': r'\textasciicircum{}',
}


def latex_escape(text):
    """Escape the LaTeX special characters in a table cell.

    Args:
        text (str): Raw cell text.

    Returns:
        str: The text with LaTeX special characters escaped.
    """
    return ''.join(_LATEX_SPECIALS.get(ch, ch) for ch in str(text))


def citekey(authors, year):
    """Build the manuscript citation key for a pole reference.

    Keys follow the manuscript's ``SurnameYYYYa`` convention: the first
    author's surname folded to ASCII, keeping the hyphens of compound surnames
    (``Swanson-Hysell2019a``), the four-digit year, and the ``a``
    disambiguation suffix. Rows whose reference is not a single publication --
    combined,
    recompiled, or working-group poles -- are keyed by hand in
    :data:`POLE_REFERENCES` instead.

    Args:
        authors (str): Summary ``POLE AUTHORS``.
        year (str): Summary ``YEAR`` (may list several years, in which case
            the first is used).

    Returns:
        str | None: The citation key, or ``None`` when the row carries neither
        an author nor a year to build one from.
    """
    surname = str(authors).split(',')[0].split(';')[0].split(' and ')[0]
    surname = unicodedata.normalize('NFKD', surname.strip())
    surname = ''.join(c for c in surname if not unicodedata.combining(c))
    surname = re.sub(r'[^A-Za-z-]', '', surname).strip('-')
    years = re.findall(r'\d{4}', str(year))
    if not surname or not years:
        return None
    return f'{surname}{years[0]}a'


def pole_citation(rockname, authors, year):
    """The ``\\citet`` command for a row's pole reference.

    Args:
        rockname (str): Summary ``ROCKNAME``, looked up in
            :data:`POLE_REFERENCES` first, which is where poles built from more
            than one study get their full credit list.
        authors (str): Summary ``POLE AUTHORS``.
        year (str): Summary ``YEAR``.

    Returns:
        tuple[str, bool]: The LaTeX citation and whether it came from
        :data:`POLE_REFERENCES`. A row that yields no key at all is
        rendered as a dash.
    """
    keys = POLE_REFERENCES.get(rockname)
    if keys:
        return r'\citet{' + ','.join(keys) + '}', True
    key = citekey(authors, year)
    return (r'\citet{' + key + '}') if key else '--', False


def rotate_pole(euler, pole_lat, pole_lon):
    """Rotate a pole about an Euler pole (Rodrigues' rotation formula).

    Verified to agree with ``pmagpy.pmag.pt_rot`` to within 1e-14 degrees; it
    is reimplemented here so this script keeps its standard-library-only
    footprint.

    Args:
        euler (list): ``[pole latitude, pole longitude, rotation angle]``, deg.
        pole_lat, pole_lon (float): The pole to rotate, in degrees.

    Returns:
        tuple[float, float]: Rotated ``(latitude, longitude)``, longitude in
        0-360 degrees east.
    """
    elat, elon, angle = euler
    el, eo, a = math.radians(elat), math.radians(elon), math.radians(angle)
    k = (math.cos(el) * math.cos(eo), math.cos(el) * math.sin(eo), math.sin(el))
    p, po = math.radians(pole_lat), math.radians(pole_lon)
    v = (math.cos(p) * math.cos(po), math.cos(p) * math.sin(po), math.sin(p))
    k_dot_v = sum(x * y for x, y in zip(k, v))
    k_cross_v = (k[1] * v[2] - k[2] * v[1],
                 k[2] * v[0] - k[0] * v[2],
                 k[0] * v[1] - k[1] * v[0])
    r = [v[i] * math.cos(a) + k_cross_v[i] * math.sin(a)
         + k[i] * k_dot_v * (1 - math.cos(a)) for i in range(3)]
    return (math.degrees(math.asin(max(-1.0, min(1.0, r[2])))),
            math.degrees(math.atan2(r[1], r[0])) % 360.0)


def paleolatitude(site_lat, site_lon, pole_lat, pole_lon):
    """Paleolatitude of a locality implied by a paleomagnetic pole (GAD).

    The angular distance between the locality and the pole is the paleo
    co-latitude, so the paleolatitude is ``90 - distance``.

    Args:
        site_lat, site_lon (float): Locality latitude and longitude (deg,
            longitude in degrees east).
        pole_lat, pole_lon (float): Pole latitude and longitude (deg).

    Returns:
        float: Paleolatitude of the locality in degrees.
    """
    slat, slon = math.radians(site_lat), math.radians(site_lon)
    plat, plon = math.radians(pole_lat), math.radians(pole_lon)
    cos_p = (math.sin(slat) * math.sin(plat)
             + math.cos(slat) * math.cos(plat) * math.cos(plon - slon))
    return 90.0 - math.degrees(math.acos(max(-1.0, min(1.0, cos_p))))


def bearing(lat1, lon1, lat2, lon2):
    """Initial great-circle bearing from one point toward another.

    Args:
        lat1, lon1 (float): Start point (degrees, longitude degrees east).
        lat2, lon2 (float): End point (degrees, longitude degrees east).

    Returns:
        float: Bearing in degrees east of north, 0-360.
    """
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlon = math.radians(lon2 - lon1)
    return math.degrees(math.atan2(
        math.sin(dlon) * math.cos(phi2),
        math.cos(phi1) * math.sin(phi2)
        - math.sin(phi1) * math.cos(phi2) * math.cos(dlon))) % 360.0


def ellipse_radius_toward(plat, plon, zeta, zdec, zinc, eta, site_lat,
                          site_lon):
    """Radius of a Kent 95% confidence ellipse in the direction of a site.

    The paleolatitude a pole implies for a site is ``90`` minus their angular
    distance, so what limits the paleolatitude is the extent of the confidence
    region along the pole-to-site great circle. For an ellipse with semi-axes
    ``zeta`` (major) and ``eta`` (minor) that extent is

        r(theta) = ((cos(theta) / zeta)^2 + (sin(theta) / eta)^2)^(-1/2),

    with ``theta`` the angle between the major axis and the bearing to the
    site. This is the quantity the paleolatitude figures plot as the vertical
    error bar for the Kent poles (``paleolat_uncertainty`` in
    ``scripts/build_paleolatitude_figure.py``); the equal-area circular
    approximation ``sqrt(zeta * eta)`` carried in the ``A95`` column would
    understate it, because the flattening-factor uncertainty elongates the
    ellipse along the site-to-pole direction.

    Args:
        plat, plon (float): Pole position (degrees, longitude degrees east).
        zeta (float): Semi-angle of the major axis, in degrees.
        zdec, zinc (float): Longitude and latitude of the major-axis
            direction, in the same frame as the pole.
        eta (float): Semi-angle of the minor axis, in degrees.
        site_lat, site_lon (float): The site, in degrees.

    Returns:
        float: The ellipse radius toward the site, in degrees.
    """
    theta = math.radians(bearing(plat, plon, site_lat, site_lon)
                         - bearing(plat, plon, zinc, zdec))
    return 1.0 / math.hypot(math.cos(theta) / zeta, math.sin(theta) / eta)


def key_pole_status(rockname):
    """Look up a pole's Buchan (2013) key-pole assignment.

    The flag is curated in :data:`KEY_POLES` rather than derived from the
    R-criteria columns -- see the commentary there for why the two axes are
    kept separate.

    Args:
        rockname (str): Summary ``ROCKNAME``.

    Returns:
        tuple[bool | None, str]: The flag (``None`` when the pole has not yet
        been assessed) and a short description of the basis for it.
    """
    if rockname not in KEY_POLES:
        return None, 'not yet assessed against the Buchan (2013) criteria'
    is_key, test, basis = KEY_POLES[rockname]
    if not is_key:
        return False, basis
    named = ', '.join(BUCHAN_TESTS.get(t, t) for t in test.split(',') if t)
    return True, f'{basis}; {named}' if named else basis


def _r_criteria(scores):
    """The seven R-criteria of Meert et al. (2020) as one compact string.

    Each criterion is rendered as 0 or 1 in R1--R7 order, e.g. ``1111011``.
    R4 is stored as the field-test letter code(s) rather than a score, so it
    contributes 1 when a test is recorded and 0 when it is not, matching how
    ``pole_tools.make_nordic_summary`` sums it into ``Rsum``. A criterion that
    was never scored is rendered as ``?``.

    Args:
        scores (list[str]): The R1--R7 cell values, in order.

    Returns:
        str: A seven-character string.
    """
    out = []
    for i, value in enumerate(scores):
        value = str(value).strip()
        if i == 3:  # R4 holds letter codes; any recorded test scores 1
            out.append('0' if value in ('', '0') else '1')
        elif value in ('0', '1'):
            out.append(value)
        else:
            out.append('?')
    return ''.join(out)


def _number(value, decimals=1):
    """Format a numeric cell to fixed decimals, or an en dash when missing.

    Negative values are written with a math-mode minus (``$-$``) rather than a
    hyphen, which LaTeX would set as a short dash in a text-mode table cell.
    """
    try:
        text = f'{float(value):.{decimals}f}'
    except (TypeError, ValueError):
        return '--'
    return text.replace('-', '$-$', 1) if text.startswith('-') else text


def _east_longitude(value):
    """Format a longitude in the 0-360 degrees east convention."""
    try:
        return f'{float(value) % 360.0:.1f}'
    except (TypeError, ValueError):
        return '--'


def _age_cell(nominal, lomagage, himagage):
    """Nominal age with its magnetization-age bounds as scripts.

    The upper bound is set as a superscript and the lower bound as a
    subscript, so the whole age fits on one line at the reduced script size
    (e.g. ``1109`` with ``1114`` above and ``1104`` below).
    """
    age = _number(nominal, 0)
    lo, hi = _number(lomagage, 0), _number(himagage, 0)
    if '--' in (lo, hi):
        return age
    return f'{age}$^{{{hi}}}_{{{lo}}}$'


def _paleolat_cell(value, error):
    """Paleolatitude with its 95% bounds as sub- and superscripts.

    The bounds are the paleolatitude plus and minus the pole's 95% confidence
    radius in the direction of the site (``A95`` for a circular confidence,
    :func:`ellipse_radius_toward` for a Kent ellipse), clipped at the poles.
    They are set inside one math group, so the minus signs of negative bounds
    are true minus signs rather than hyphens.

    Args:
        value (float): The paleolatitude, in degrees.
        error (float | None): Its 95% uncertainty in degrees; ``None`` leaves
            the bounds off.

    Returns:
        str: The formatted cell.
    """
    text = _number(value, 1)
    if text == '--' or error is None:
        return text
    hi = min(90.0, float(value) + float(error))
    lo = max(-90.0, float(value) - float(error))
    return f'{text}$^{{{hi:.1f}}}_{{{lo:.1f}}}$'


TEX_PREAMBLE_NOTE = r"""% Auto-generated by data/nordic_summaries/combine_nordic_summaries.py
% Source: data/nordic_summaries/nordic_summaries_combined.csv
% Do not edit by hand -- rerun the script instead.
%
% Requires: longtable, booktabs, amssymb (\checkmark), natbib (\citet).
"""

TEX_CAPTION = (
    r'Compiled paleomagnetic poles for Laurentia spanning ca.\ '
    r'{oldest}--{youngest}~Ma, beginning after the amalgamation of Laurentia '
    r'and running through the remainder of the Proterozoic, listed from '
    r'youngest to oldest. Each pole is recalculated from '
    r"site-level data in the accompanying compilation. ``Age (Ma)'' gives the "
    r'nominal age of magnetization, with its upper bound as a superscript '
    r'and its lower bound as a subscript. '
    r"``Rating'' is the Nordic reliability grade (A/B). ``Key pole'' "
    r'($\checkmark$) marks poles meeting the key-pole criteria of '
    r'\citet{Buchan2013a}: a precisely determined age (U--Pb, or occasionally '
    r'Ar--Ar, typically within $\pm$10~Myr), a good-quality determination in '
    r'which the primary remanence is isolated by stepwise demagnetization and '
    r'secular variation is largely averaged, and a positive field test '
    r'establishing that the remanence is primary rather than solely older '
    r'than a subsequent geologic event. This flag is assessed independently '
    r'of the Nordic grade and of the R-criteria scores of \citet{Meert2020a}, '
    r'which are reported in the accompanying compilation; a blank marks a '
    r'pole assessed against these criteria and found not to meet them, and a '
    r'dash marks one not yet assessed. Longitudes run $0$--$360^{\circ}$E. '
    r"``Duluth paleolat'' is "
    r'the paleolatitude of Duluth, Minnesota ($46.8^{\circ}$N, '
    r'$267.9^{\circ}$E) implied by the pole, with the bounds set by $A_{95}$ '
    r'given as a superscript and a subscript. Site and pole positions are '
    r'present-day; the Greenland and Scotland poles are rotated back into '
    r'Laurentia coordinates before their Duluth paleolatitude is computed, so '
    r'that column is comparable across the whole table.')

# The three text columns and the reference column are ``p`` columns whose
# cells start with \raggedright, so unit names and citations are not stretched
# to fill the measure. Rows are therefore ended with \tabularnewline rather
# than \\, which \raggedright would otherwise clobber in the final column; this
# keeps the table free of any package beyond longtable and booktabs.
# The Duluth column is 1.6cm because its widest cell -- a negative paleolatitude
# with both bounds negative, e.g. $-$41.4$^{-28.5}_{-54.3}$ -- sets 41.3pt and
# cannot be broken across lines.
# 'Key pole' is a centred column with no width of its own, so its header used
# to set it: the single line "Key pole" made it ~1.5cm wide to hold a
# checkmark. Stacking the header on two lines (see TEX_HEADER_ROW) narrows it
# to the width of "Key", and the space goes to the reference column, whose
# multi-author \citet keys are what most need the room. The Duluth column
# gains a little too, clearing the negative-with-both-bounds-negative cells
# that were fractionally overfull at 1.6cm.
# Column-unit labels. Stacking the unit under the column name with \shortstack
# keeps the column no wider than the longer of the two lines (a `c`/`r` column
# is sized by its header here), and lets both tables state their conventions in
# the headers rather than in the caption.
DEG = r'$^{\circ}$'
DEG_E = DEG + 'E'
DEG_N = DEG + 'N'


def unit_header(label, unit):
    """Bold column label with its unit stacked beneath it.

    Args:
        label (str): Column name, e.g. ``'Plon'``.
        unit (str): Unit, e.g. :data:`DEG_E`; wrapped in parentheses.

    Returns:
        str: A ``\shortstack`` cell for the header row.
    """
    return r'\shortstack{\textbf{' + label + r'}\\\textbf{(' + unit + r')}}'


TEX_COLUMN_SPEC = r'p{1.8cm}p{3.1cm}p{1.5cm}cccrrrrrp{1.7cm}p{3.2cm}'
RAGGED = r'\raggedright '
ROW_END = r' \tabularnewline'

TEX_HEADER_ROW = ' & '.join([
    r'\textbf{Terrane}', r'\textbf{Unit}', r'\textbf{Age (Ma)}',
    r'\textbf{Rating}', r'\textbf{R1--R7}',
    r'\shortstack{\textbf{Key}\\\textbf{pole}}',
    unit_header('Site lon', DEG_E), unit_header('Site lat', DEG_N),
    unit_header('Plon', DEG_E), unit_header('Plat', DEG_N),
    unit_header('A95', DEG),
    r'\textbf{Duluth}\newline\textbf{paleolat}\newline\textbf{(' + DEG_N + r')}',
    r'\textbf{Pole reference}',
]) + ROW_END


# A hairline between body rows. The table is 13 columns wide with many two-line
# cells, so the eye has a long way to travel from unit name to reference; a rule
# lighter than \midrule keeps rows tied together without the heavy banding a
# full-weight rule would give across 64 rows. Set to None for plain booktabs
# spacing (no rules between rows).
ROW_RULE = r'\midrule[0.1pt]'


def _interleave_rules(body):
    """Body rows separated by ROW_RULE, or unchanged if ROW_RULE is None."""
    if not ROW_RULE or not body:
        return body
    out = []
    for row in body[:-1]:
        out.extend([row, ROW_RULE])
    out.append(body[-1])
    return out


def build_pole_table_tex(header, rows):
    """Render the combined summary rows as a LaTeX ``longtable``.

    Args:
        header (list[str]): The Nordic column header.
        rows (list[list[str]]): Summary rows, in the order they are tabulated.

    Returns:
        tuple[str, dict]: The LaTeX source, and a report with the
        ``key_poles`` (list of ROCKNAMEs flagged),
        ``unassessed`` (ROCKNAMEs with no entry in :data:`KEY_POLES`),
        ``unkeyed`` (ROCKNAMEs with no author/year to build a key from), and
        ``unused_overrides`` (:data:`POLE_REFERENCES` entries matching no
        row, which are stale and should be renamed or dropped).
    """
    col = {name: header.index(name) for name in
           ('Terrane', 'ROCKNAME', 'nominal age', 'lomagage', 'himagage',
            'Grade', 'SLAT', 'SLONG', 'PLAT', 'PLONG', 'A95',
            'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7',
            'POLE AUTHORS', 'YEAR', 'TITLE')}

    def cell(row, name):
        i = col[name]
        return row[i].strip() if i < len(row) else ''

    body, key_poles, unassessed, unkeyed = [], [], [], []
    overrides_used = set()
    basis_notes, caveats = [], []
    for row in rows:
        rockname = cell(row, 'ROCKNAME')
        is_key, basis = key_pole_status(rockname)
        if is_key:
            key_poles.append(rockname)
            basis_notes.append(f'%   {rockname}: {basis}')
        if 'NOTE:' in basis:
            caveats.append(f'%   {rockname}: {basis.split("NOTE:", 1)[1].strip()}')
        elif is_key is None:
            unassessed.append(rockname)

        reference, from_override = pole_citation(
            rockname, cell(row, 'POLE AUTHORS'), cell(row, 'YEAR'))
        if from_override:
            overrides_used.add(rockname)
        elif reference == '--':
            unkeyed.append(rockname)

        try:
            plat, plon = float(cell(row, 'PLAT')), float(cell(row, 'PLONG'))
            euler = TERRANE_EULER_POLES.get(cell(row, 'Terrane'))
            if euler is not None:   # rift back into Laurentia coordinates
                plat, plon = rotate_pole(euler, plat, plon)
            try:
                a95 = float(cell(row, 'A95'))
            except ValueError:  # a pole tabulated without a confidence radius
                a95 = None
            duluth = _paleolat_cell(
                paleolatitude(DULUTH_LAT, DULUTH_LON, plat, plon), a95)
        except ValueError:  # a row without a pole position
            duluth = '--'

        body.append(' & '.join([
            RAGGED + latex_escape(cell(row, 'Terrane')),
            RAGGED + latex_escape(rockname),
            RAGGED + _age_cell(cell(row, 'nominal age'),
                               cell(row, 'lomagage'), cell(row, 'himagage')),
            latex_escape(cell(row, 'Grade')) or '--',
            _r_criteria([cell(row, f'R{i}') for i in range(1, 8)]),
            {True: r'$\checkmark$', False: '', None: '--'}[is_key],
            _east_longitude(cell(row, 'SLONG')),
            _number(cell(row, 'SLAT')),
            _east_longitude(cell(row, 'PLONG')),
            _number(cell(row, 'PLAT')),
            _number(cell(row, 'A95')),
            r'\hfill ' + duluth,
            RAGGED + reference,
        ]) + ROW_END)

    notes = TEX_PREAMBLE_NOTE
    if basis_notes:
        notes += '%\n% Basis for each key-pole flag:\n'
        notes += '\n'.join(basis_notes) + '\n'
    if caveats:
        notes += ('%\n% CAVEATS -- flags carried for fidelity to Buchan (2013) '
                  'that warrant\n%           discussion in the text:\n')
        notes += '\n'.join(caveats) + '\n'

    # the caption states the interval the table covers; take it from the data
    # rather than restating it by hand as poles are added at either end
    ages = [float(cell(row, 'nominal age')) for row in rows
            if cell(row, 'nominal age')]
    caption = (TEX_CAPTION
               .replace('{oldest}', f'{max(ages):.0f}' if ages else '?')
               .replace('{youngest}', f'{min(ages):.0f}' if ages else '?'))

    tex = '\n'.join([
        notes.rstrip('\n'),
        r'\begingroup',
        r'\footnotesize',
        r'\setlength{\tabcolsep}{3pt}',
        # The age, paleolatitude and f cells carry sub- and superscripts that
        # reach beyond the normal row box, so rows set at the default height
        # read as crowded. Stretching them is local to this \begingroup.
        r'\renewcommand{\arraystretch}{1.25}',
        r'\begin{longtable}{' + TEX_COLUMN_SPEC + '}',
        r'\caption{' + caption + r'}\label{tab:poles} \\',
        r'\toprule',
        TEX_HEADER_ROW,
        r'\midrule',
        r'\endfirsthead',
        r'\toprule',
        TEX_HEADER_ROW,
        r'\midrule',
        r'\endhead',
        r'\midrule',
        r'\multicolumn{13}{r}{\textit{continued on next page}}' + ROW_END,
        r'\endfoot',
        r'\bottomrule',
        r'\endlastfoot',
        *_interleave_rules(body),
        r'\end{longtable}',
        r'\endgroup',
        '',
    ])
    return tex, {'key_poles': key_poles, 'unassessed': unassessed,
                 'unkeyed': unkeyed,
                 'unused_overrides': sorted(set(POLE_REFERENCES)
                                            - overrides_used)}


def write_pole_table(header, rows, summary_dir=SUMMARY_DIR):
    """Write ``pole_table.tex`` next to the combined CSV and report on it.

    Args:
        header (list[str]): The Nordic column header.
        rows (list[list[str]]): Combined summary rows, already sorted.
        summary_dir (str): Directory the table is written into.

    Returns:
        str: Path to the LaTeX table written to disk.
    """
    tex, report = build_pole_table_tex(header, rows)
    tex_path = os.path.join(summary_dir, TEX_FILENAME)
    with open(tex_path, 'w', encoding='utf-8') as fh:
        fh.write(tex)
    print(f'Wrote {len(rows)}-row LaTeX pole table to {tex_path}')
    print(f'  key poles: {len(report["key_poles"])} of {len(rows)}')
    if report['unassessed']:
        print(f'  -I- {len(report["unassessed"])} pole(s) not yet assessed '
              'against the Buchan (2013) criteria (shown as -- in the table); '
              'record a call for each in KEY_POLES:')
        for rockname in report['unassessed']:
            print(f'      {rockname}')
    if report['unkeyed']:
        print(f'  -W- {len(report["unkeyed"])} pole(s) carry no author/year to '
              'build a citation key from (shown as -- in the table): '
              + ', '.join(report['unkeyed']))
    if report['unused_overrides']:
        print(f'  -W- {len(report["unused_overrides"])} POLE_REFERENCES '
              'entry(ies) match no row and are stale: '
              + ', '.join(report['unused_overrides']))
    return tex_path


def main(argv=None):
    """Combine the summaries and render the manuscript pole table."""
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--no-tex', action='store_true',
                        help='only write the combined CSV')
    args = parser.parse_args(argv)

    combined_path = combine_summaries()
    if args.no_tex:
        return combined_path
    with open(combined_path, encoding='utf-8-sig', newline='') as fh:
        records = list(csv.reader(fh))
    write_pole_table(records[0], records[1:])
    return combined_path


if __name__ == '__main__':
    main()
