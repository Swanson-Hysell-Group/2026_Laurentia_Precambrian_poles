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

Usage:
    python combine_nordic_summaries.py
"""

import csv
import glob
import os
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

    # back-fill blank cells from the compilation (e.g. the legacy Q criteria)
    filled, unmatched = backfill_from_compilation(header, rows)

    with open(combined_path, 'w', encoding='utf-8-sig', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    print(f'Combined {len(rows)} summaries into {combined_path}')
    if filled:
        print(f'Back-filled {filled} blank cell(s) from '
              f'{os.path.basename(COMPILATION_PATH)}')
    if unmatched:
        print(f'-W- {len(unmatched)} row(s) had no ROCKNAME match in '
              f'{os.path.basename(COMPILATION_PATH)}: '
              f'{", ".join(sorted(set(unmatched)))}')
    return combined_path


if __name__ == '__main__':
    combine_summaries()
