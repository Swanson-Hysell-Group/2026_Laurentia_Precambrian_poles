"""Combine the per-notebook Kent summaries of the sedimentary poles.

The main compilation (``nordic_summaries_combined.csv`` and ``pole_table.tex``)
reports every pole **as measured** — the uncorrected, f = 1 position with a
circular A95. Detrital remanence in sedimentary rocks is shallowed during
compaction, so for a sedimentary unit that position implies a paleolatitude that
is a minimum. This script assembles the companion table that carries the
inclination-shallowing-corrected position of each sedimentary pole, summarized
as a Kent (1982) distribution whose 95% confidence ellipse propagates the
uncertainty in the flattening factor f following Pierce et al. (2022).

Each sedimentary pole's notebook writes a single-row CSV into ``kent/`` via
``pole_tools.save_kent_summary`` (columns ``pole_tools.KENT_COLUMNS``), by one of
two routes:

- **quantified** — the study determined f from the directional distribution
  itself, by elongation/inclination (E/I; Tauxe & Kent, 2004) or by SVEI (Tauxe
  et al., 2024), and propagated it to the pole with ``ipmag.find_ei_kent`` /
  ``ipmag.find_svei_kent``.
- **compilation** — no f could be determined from the data (a published mean
  pole, or too few directions), so f is resampled from the compiled distribution
  of measured flattening factors of Pierce et al. (2022) with
  ``ipmag.find_compilation_kent``.

This script concatenates those rows into ``kent_poles_combined.csv`` (sorted by
nominal age) and renders ``kent_pole_table.tex``, the manuscript's Kent-pole
table. Every unit in :data:`SEDIMENTARY_POLES` must have a Kent summary and a
row in the main combined table, and no other unit may appear in ``kent/``; the
script fails loudly otherwise, so a sedimentary pole cannot be silently added to
the compilation without a Kent treatment.

The Kent positions written here are what the pole map, the paleolatitude
figures, and the APWP figure plot for these units in place of the uncorrected
positions (see ``scripts/build_pole_map.py``).

Usage:
    python build_kent_poles.py
    python build_kent_poles.py --no-tex
"""

import argparse
import csv
import os
import sys

SUMMARY_DIR = os.path.dirname(os.path.abspath(__file__))
KENT_DIR = os.path.join(SUMMARY_DIR, 'kent')
KENT_COMBINED_FILENAME = 'kent_poles_combined.csv'
KENT_TEX_FILENAME = 'kent_pole_table.tex'

sys.path.insert(0, SUMMARY_DIR)
from combine_nordic_summaries import (  # noqa: E402
    COMBINED_FILENAME, DULUTH_LAT, DULUTH_LON, TERRANE_EULER_POLES,
    DEG, DEG_E, DEG_N, _age_cell, _east_longitude, _number, _paleolat_cell,
    ellipse_radius_toward, latex_escape, paleolatitude, rotate_pole,
    unit_header)

# Every sedimentary pole in the compilation, ROCKNAME -> the notebook that
# assesses it (and so writes its Kent summary into kent/). A unit belongs here
# when its remanence is a detrital or early-diagenetic sedimentary remanence
# whose inclination is subject to compaction shallowing — the Chuar, Uinta
# Mountain, Torridon, Stoer, Oronto (Nonesuch, Freda) and Jacobsville clastic
# successions, and the Belt-Purcell redbeds. Volcanic units interbedded with
# sediments (e.g. the Purcell Lava, the Cardenas Basalts) are not sedimentary
# poles and are absent.
SEDIMENTARY_POLES = {
    'Chuar Group (combined)': '755_Chuar_Group',
    'Uinta Mountain Group': '759_Uinta_Mountain_Group',
    'Torridon Group': '975_Torridon',
    'Jacobsville Formation': '990_Jacobsville',
    'Upper Freda Formation': '1045_Upper_Freda',
    'Lower Freda Formation': '1075_Lower_Freda',
    'Nonesuch Formation': '1078_Nonesuch',
    'Stoer Group': '1199_Stoer',
    'Pilcher, Garnet Range, Libby': '1385_Garnet_Range',
    'McNamara': '1392_McNamara',
    'Snowslip': '1449_Snowslip',
    'Spokane': '1458_Spokane',
}


def read_rows(path):
    """Read a CSV written by ``save_kent_summary`` as (header, rows)."""
    with open(path, encoding='utf-8-sig', newline='') as fh:
        records = list(csv.reader(fh))
    return records[0], records[1:]


def combine_kent_summaries(kent_dir=KENT_DIR, summary_dir=SUMMARY_DIR):
    """Concatenate the per-notebook Kent summaries into one CSV.

    Args:
        kent_dir (str): Directory holding the per-notebook Kent summary CSVs.
        summary_dir (str): Directory the combined CSV is written into.

    Returns:
        tuple[str, list[str], list[list[str]]]: Path written, header, and the
        rows in the order they were written (sorted by nominal age).

    Raises:
        SystemExit: If a unit in :data:`SEDIMENTARY_POLES` has no Kent summary,
            if ``kent/`` holds a unit that is not in it, or if a Kent row does
            not correspond to a row of the main combined table.
    """
    if not os.path.isdir(kent_dir):
        raise SystemExit(f'-E- no Kent summaries found in {kent_dir}; run the '
                         'sedimentary pole notebooks first')

    header, rows, seen = None, [], {}
    for path in sorted(os.listdir(kent_dir)):
        if not path.endswith('.csv'):
            continue
        file_header, file_rows = read_rows(os.path.join(kent_dir, path))
        if header is None:
            header = file_header
        elif file_header != header:
            raise SystemExit(f'-E- {path} has a different header than the other '
                             'Kent summaries; regenerate it from the notebook')
        for row in file_rows:
            rockname = row[header.index('ROCKNAME')].strip()
            if rockname in seen:
                raise SystemExit(f'-E- {rockname} appears in both '
                                 f'{seen[rockname]} and {path}')
            seen[rockname] = path
            rows.append(row)

    missing = [name for name in SEDIMENTARY_POLES if name not in seen]
    if missing:
        raise SystemExit(
            '-E- no Kent summary for: ' + ', '.join(missing)
            + '\n    run the notebook(s) '
            + ', '.join(SEDIMENTARY_POLES[name] for name in missing)
            + ' to write them into kent/')
    extra = [name for name in seen if name not in SEDIMENTARY_POLES]
    if extra:
        raise SystemExit(
            '-E- Kent summaries for units that are not in SEDIMENTARY_POLES: '
            + ', '.join(extra) + '\n    add them there (or remove the CSV) so '
            'the two tables stay in step')

    # Every Kent pole must correspond to a pole in the main compilation; a
    # sedimentary unit dropped from the compilation (Grade C, or excluded)
    # must not linger here.
    combined_path = os.path.join(summary_dir, COMBINED_FILENAME)
    main_header, main_rows = read_rows(combined_path)
    main_names = {row[main_header.index('ROCKNAME')].strip() for row in main_rows}
    orphaned = sorted(set(seen) - main_names)
    if orphaned:
        raise SystemExit(
            '-E- Kent summaries for units absent from '
            f'{COMBINED_FILENAME}: ' + ', '.join(orphaned))

    def age(row):
        try:
            return float(row[header.index('nominal age')])
        except ValueError:
            return float('inf')

    rows.sort(key=age)
    out_path = os.path.join(summary_dir, KENT_COMBINED_FILENAME)
    with open(out_path, 'w', encoding='utf-8-sig', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    print(f'Wrote {len(rows)} Kent pole rows to {out_path}')
    return out_path, header, rows


# --- LaTeX Kent pole table --------------------------------------------------
# Citation keys used in the caption. These are methodological citations, the
# same on every row, so they are written out here rather than built from a row's
# authors as the per-pole references in pole_table.tex are. Like those, they
# follow the manuscript's AuthorYYYYa key convention and must exist in its bib.
CAPTION_CITEKEYS = ('Kent1982a', 'Pierce2022a', 'Tauxe2004a', 'Tauxe2024a')

TEX_PREAMBLE_NOTE = r"""% Auto-generated by data/nordic_summaries/build_kent_poles.py
% Source: data/nordic_summaries/kent_poles_combined.csv
% Do not edit by hand -- rerun the script instead.
%
% Requires: longtable, booktabs, natbib (\citet, \citep).
% Caption cites: {citekeys}
"""

# The caption as authored for the manuscript. It leans on the methods section,
# which already introduces the Pierce et al. (2022) approach and the orientation
# of the Kent ellipse, and so restates only what a reader needs to parse the
# columns -- including that the compilation f is a property of that compilation
# rather than a determination for the individual unit.
TEX_CAPTION = (
    r'Inclination-shallowing-corrected mean poles for the sedimentary units of '
    r'the compilation; Table~\ref{tab:poles} gives the same poles as measured. '
    r'Each is the mean of a \citet{Kent1982a} distribution whose 95\% '
    r'confidence ellipse propagates the flattening-factor ($f$) uncertainty '
    r'following \citet{Pierce2022a}. The ellipse is given as the mean pole, and '
    r'the direction (longitude, latitude) and semi-angle of each of its axes, '
    r'$\zeta_{95}$ for the major and $\eta_{95}$ for the minor. '
    r"``Correction'' gives the source of $f$: E/I \citep{Tauxe2004a} and SVEI "
    r'\citep{Tauxe2024a} determine it from the directional distribution of the '
    r'unit itself, whereas the compilation approach resamples flattening factors '
    r'as compiled by \citet{Pierce2022a} (the median and 95$\%$ bounds are what '
    r'is shown for those $f$ values). Scotland poles (Torridon, Stoer) are '
    r'rotated into Laurentia coordinates before the paleolatitude of Duluth, '
    r'Minnesota ($46.8^{\circ}$N, $267.9^{\circ}$E) is computed. That '
    r'paleolatitude is repeated for the uncorrected pole of '
    r'Table~\ref{tab:poles} alongside the corrected one; their bounds follow '
    r'from $A_{95}$ and from the Kent ellipse along the Duluth-to-pole great '
    r'circle.')

# Columns: the four identifying/left-hand columns, then the Kent mean pole, then
# each ellipse axis as (longitude, latitude, semi-angle), then the two Duluth
# paleolatitudes. The mean plus both axes fully specify the ellipse, so the table
# reproduces it without reference to anything else. Neither the uncorrected pole
# nor the terrane is repeated here, since Table 1 carries both; the terrane is
# still read from the row to rotate the Scotland poles for the Duluth column. The f column sets its 95% bounds as
# scripts on a single line (see _f_cell), so it needs no paragraph column.
TEX_COLUMN_SPEC = (r'p{4.0cm}p{1.4cm}p{2.0cm}c'
                   + 'rr' + 'rrr' + 'rrr' + 'rr')
N_COLUMNS = 14
RAGGED = r'\raggedright '
ROW_END = r' \tabularnewline'

# Two-row header: a spanning row grouping the pole and ellipse-axis blocks, then
# the column labels. \cmidrule(lr) keeps the rules clear of the neighboring
# groups.
TEX_HEADER_ROWS = '\n'.join([
    (r'& & & & \multicolumn{2}{c}{\textbf{Kent pole}} & '
     r'\multicolumn{3}{c}{\textbf{Major axis}} & '
     r'\multicolumn{3}{c}{\textbf{Minor axis}} & '
     r'\multicolumn{2}{c}{\textbf{Duluth paleolat}}' + ROW_END),
    (r'\cmidrule(lr){5-6}\cmidrule(lr){7-9}\cmidrule(lr){10-12}'
     r'\cmidrule(lr){13-14}'),
    (' & '.join([
        r'\textbf{Unit}', r'\textbf{Age (Ma)}', r'\textbf{Correction}',
        r'\textbf{$f$}',
        unit_header('Plon', DEG_E), unit_header('Plat', DEG_N),
        unit_header('lon', DEG_E), unit_header('lat', DEG_N),
        unit_header(r'$\zeta_{95}$', DEG),
        unit_header('lon', DEG_E), unit_header('lat', DEG_N),
        unit_header(r'$\eta_{95}$', DEG),
        unit_header('uncorr.', DEG_N), unit_header('corr.', DEG_N),
     ]) + ROW_END),
])

# Short label for the correction method in the table; the full citation is in
# the caption rather than repeated on every row.
METHOD_LABELS = {
    'E/I (Tauxe & Kent, 2004)': 'E/I',
    'SVEI (Tauxe et al., 2024)': 'SVEI',
    'compilation f (Pierce et al., 2022)': 'compilation',
}


def _f_cell(f, f_low, f_high):
    """Flattening factor with its 95% bounds as scripts, matching the age cell.

    The upper bound is set as a superscript and the lower bound as a subscript
    so the whole entry fits on one line (e.g. ``0.62`` with ``0.91`` above and
    ``0.43`` below), as ``combine_nordic_summaries._age_cell`` does for ages.
    Returns a dash when no flattening factor is recorded.
    """
    if not str(f).strip():
        return '--'
    lo, hi = _number(f_low, 2), _number(f_high, 2)
    if '--' in (lo, hi):
        return _number(f, 2)
    return f'{_number(f, 2)}$^{{{hi}}}_{{{lo}}}$'


def _duluth(terrane, plat, plon, a95='', ellipse=None):
    """Duluth paleolatitude with its 95% bounds, rotated into Laurentia first.

    Args:
        terrane (str): Terrane label, used to look up the rotation.
        plat, plon (str | float): The pole, in present-day coordinates.
        a95 (str | float): Circular 95% confidence radius; sets the bounds for
            the uncorrected pole.
        ellipse (tuple | None): ``(zeta, zdec, zinc, eta)`` of the Kent 95%
            ellipse. When given, the bounds come from the ellipse projected
            onto the pole-to-Duluth great circle rather than from ``a95``,
            since the flattening-factor uncertainty elongates it along that
            direction.

    Returns:
        str: The formatted table cell.
    """
    try:
        plat, plon = float(plat), float(plon)
    except (TypeError, ValueError):
        return '--'
    euler = TERRANE_EULER_POLES.get(terrane)
    if euler is not None:
        plat, plon = rotate_pole(euler, plat, plon)

    error = None
    if ellipse is not None:
        try:
            zeta, zdec, zinc, eta = (float(v) for v in ellipse)
        except (TypeError, ValueError):
            zeta = eta = 0.0
        if zeta > 0 and eta > 0:
            if euler is not None:  # a rigid rotation carries the ellipse too
                zinc, zdec = rotate_pole(euler, zinc, zdec)
            error = ellipse_radius_toward(plat, plon, zeta, zdec, zinc, eta,
                                          DULUTH_LAT, DULUTH_LON)
    if error is None:
        try:
            error = float(a95)
        except (TypeError, ValueError):
            error = None
    return _paleolat_cell(paleolatitude(DULUTH_LAT, DULUTH_LON, plat, plon),
                          error)


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


def build_kent_table_tex(header, rows):
    """Render the combined Kent rows as a LaTeX ``longtable``.

    Args:
        header (list[str]): The Kent column header.
        rows (list[list[str]]): Kent rows, in the order they are tabulated.

    Returns:
        str: The LaTeX source.
    """
    col = {name: header.index(name) for name in header}

    def cell(row, name):
        i = col[name]
        return row[i].strip() if i < len(row) else ''

    body = []
    for row in rows:
        terrane = cell(row, 'Terrane')
        method = cell(row, 'f method')
        body.append(' & '.join([
            RAGGED + latex_escape(cell(row, 'ROCKNAME')),
            RAGGED + _age_cell(cell(row, 'nominal age'), cell(row, 'lomagage'),
                               cell(row, 'himagage')),
            RAGGED + latex_escape(METHOD_LABELS.get(method, method)),
            _f_cell(cell(row, 'f'), cell(row, 'f low'), cell(row, 'f high')),
            # Kent mean pole, then each 95% ellipse axis as its direction
            # (longitude, latitude) and semi-angle -- together these reproduce
            # the ellipse. The uncorrected pole is in Table 1, not repeated.
            _east_longitude(cell(row, 'PLONG')),
            _number(cell(row, 'PLAT')),
            _east_longitude(cell(row, 'Zdec')),
            _number(cell(row, 'Zinc')),
            _number(cell(row, 'Zeta')),
            _east_longitude(cell(row, 'Edec')),
            _number(cell(row, 'Einc')),
            _number(cell(row, 'Eta')),
            _duluth(terrane, cell(row, 'PLAT uncorrected'),
                    cell(row, 'PLONG uncorrected'),
                    a95=cell(row, 'A95 uncorrected')),
            _duluth(terrane, cell(row, 'PLAT'), cell(row, 'PLONG'),
                    a95=cell(row, 'A95'),
                    ellipse=(cell(row, 'Zeta'), cell(row, 'Zdec'),
                             cell(row, 'Zinc'), cell(row, 'Eta'))),
        ]) + ROW_END)

    ncol = N_COLUMNS
    return '\n'.join([
        TEX_PREAMBLE_NOTE.format(citekeys=', '.join(CAPTION_CITEKEYS)).rstrip('\n'),
        r'\begingroup',
        r'\footnotesize',
        r'\setlength{\tabcolsep}{3pt}',
        # The age, paleolatitude and f cells carry sub- and superscripts that
        # reach beyond the normal row box, so rows set at the default height
        # read as crowded. Stretching them is local to this \begingroup.
        r'\renewcommand{\arraystretch}{1.25}',
        r'\begin{longtable}{' + TEX_COLUMN_SPEC + '}',
        r'\caption{' + TEX_CAPTION + r'}\label{tab:kent_poles} \\',
        r'\toprule',
        TEX_HEADER_ROWS,
        r'\midrule',
        r'\endfirsthead',
        r'\toprule',
        TEX_HEADER_ROWS,
        r'\midrule',
        r'\endhead',
        r'\midrule',
        r'\multicolumn{' + str(ncol) + r'}{r}{\textit{continued on next page}}'
        + ROW_END,
        r'\endfoot',
        r'\bottomrule',
        r'\endlastfoot',
        *_interleave_rules(body),
        r'\end{longtable}',
        r'\endgroup',
        '',
    ])


def write_kent_table(header, rows, summary_dir=SUMMARY_DIR):
    """Write ``kent_pole_table.tex`` next to the combined Kent CSV."""
    tex_path = os.path.join(summary_dir, KENT_TEX_FILENAME)
    with open(tex_path, 'w', encoding='utf-8') as fh:
        fh.write(build_kent_table_tex(header, rows))
    quantified = sum(1 for row in rows
                     if 'compilation' not in row[header.index('f method')])
    print(f'Wrote {len(rows)}-row LaTeX Kent pole table to {tex_path}')
    print(f'  {quantified} of {len(rows)} poles carry a study determination of '
          'f (E/I or SVEI); the rest resample the Pierce et al. (2022) '
          'compilation')
    return tex_path


def main(argv=None):
    """Combine the Kent summaries and render the Kent pole table."""
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--no-tex', action='store_true',
                        help='only write the combined CSV')
    args = parser.parse_args(argv)

    combined_path, header, rows = combine_kent_summaries()
    if not args.no_tex:
        write_kent_table(header, rows)
    return combined_path


if __name__ == '__main__':
    main()
