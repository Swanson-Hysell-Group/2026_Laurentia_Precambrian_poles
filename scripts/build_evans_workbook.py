"""Export the compilation as an Excel workbook in the Nordic-workshop format.

The workshop compilation is maintained by Dave Evans as a single spreadsheet with
a fixed 71-column layout. ``data/nordic_summaries`` was built to that layout, so
this script is a faithful transcription rather than a reshaping: the header row
is copied verbatim from the reference workbook and each pole is written beneath
it in the same column order.

The workbook is a complete Laurentia record, assembled from two sources:

- Poles at or below ``PUBLISH_AGE_MAX`` come from this project's site-level
  summaries. These are the poles recomputed from site VGPs and documented in a
  notebook apiece, and are the set the accompanying manuscript tabulates.
- Older poles come from ``data/older_Laurentia_poles.csv``, which
  ``extract_older_poles.py`` lifts from the workshop workbook in the same
  layout. They have not been recreated at the site level here, so they are
  carried across unchanged -- same layout, same values, no transcription step to
  get wrong.

The two ranges are disjoint at ``PUBLISH_AGE_MAX``, and any older row whose
ROCKNAME matches a site-level pole is dropped as superseded, so no pole is
listed twice.

Two properties of the reference layout are deliberately preserved even though
they would be unusual in a file written from scratch:

- Column names repeat. ``INCf``/``PLATf``/``PLONf``/``DPf``/``DMf``/``A95f``
  appear twice, once for each fixed flattening-factor scenario, and three of the
  Q-criterion columns are headed with the integers 24, 10 and 16. Columns are
  therefore addressed positionally throughout, never by name.
- ``age diff`` holds a formula (``=BI{row}-BH{row}``, i.e. himagage minus
  lomagage) rather than a value, matching the reference workbook so the column
  keeps recomputing if ages are edited in place. Rows carried over from the
  reference workbook have their formula rewritten for the new row number.

Note on flattening factors: this layout carries fixed *f* corrections only. The
inclination-shallowing treatment used in this repository and its manuscript --
Kent distributions that propagate the uncertainty in *f* (see
``build_kent_poles.py``) -- has no home here, and the ``INCf`` blocks are the
fixed-*f* values as computed for the compilation. The two are not
interchangeable, and the Kent results are reported separately.

    python scripts/build_evans_workbook.py
"""

import csv
import os
import re

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SUMMARY_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                           "nordic_summaries_combined.csv")
BUILD_POLE_MAP = os.path.join(HERE, "build_pole_map.py")
OUT_XLSX = os.path.join(ROOT, "data", "nordic_summaries",
                        "Laurentia_compilation_Nordic_format.xlsx")

# The pre-compilation poles, committed in the Nordic layout by
# ``extract_older_poles.py``. Both halves of the compilation are therefore in
# the repository and the build needs no file outside it.
OLDER_CSV = os.path.join(ROOT, "data", "older_Laurentia_poles.csv")

SHEET_TITLE = "Laurentia 2026"
# Matches the reference workbook, which freezes above the Q-criterion block so
# the identifying columns stay visible while scrolling through the scores.
FREEZE_PANES = "AP2"

# 1-based column indices in the fixed layout.
COL_TERRANE = 1
COL_ROCKNAME = 2
COL_AGE_DIFF = 59     # BG, written as a formula
COL_NOMINAL_AGE = 58  # BF
COL_LOMAGAGE = 60     # BH
COL_HIMAGAGE = 61     # BI


def publish_age_max():
    """The age at which the site-level compilation hands over to older poles.

    Read from ``build_pole_map.py`` rather than duplicated, so the workbook and
    the compilation page cannot drift apart. That module is not imported because
    it pulls in the plotting stack for a single constant.

    Returns:
        float: The cutoff in Ma.
    """
    with open(BUILD_POLE_MAP, encoding="utf-8") as fh:
        m = re.search(r"^PUBLISH_AGE_MAX\s*=\s*([0-9.]+)", fh.read(), re.M)
    if not m:
        raise SystemExit("PUBLISH_AGE_MAX not found in build_pole_map.py")
    return float(m.group(1))


def coerce(value):
    """Convert a CSV field to the type the spreadsheet should carry.

    Numeric fields become int or float so they sort and plot as numbers; fields
    such as ``%REV`` that legitimately hold text ("0or100", "MIXED") are left as
    strings, and empty fields become ``None`` so the cell is genuinely blank
    rather than an empty string.

    Text is returned unchanged rather than stripped. Some of the carried-over
    prose fields end in a space or a newline in the workshop workbook, and the
    point of carrying those rows is that they come across exactly as they are.

    Args:
        value (str): Raw CSV field.

    Returns:
        int | float | str | None: The coerced value.
    """
    stripped = value.strip()
    if not stripped:
        return None
    if re.fullmatch(r"-?\d+", stripped):
        return int(stripped)
    if re.fullmatch(r"-?\d*\.\d+([eE][-+]?\d+)?", stripped):
        return float(stripped)
    return value


def read_layout(path):
    """Read one half of the compilation.

    Args:
        path (str): CSV in the Nordic 71-column layout.

    Returns:
        tuple[list, list[list]]: The header and the coerced data rows.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        # The header is coerced too: three of the Q-criterion columns are headed
        # with the integers 24, 10 and 16 in the reference workbook, and writing
        # them back as text would not reproduce its layout.
        header = [coerce(v) for v in next(reader)]
        rows = [[coerce(v) for v in r] for r in reader
                if any(f.strip() for f in r)]
    return header, rows


def main():
    cutoff = publish_age_max()
    header, our_rows = read_layout(SUMMARY_CSV)
    older_header, carried = read_layout(OLDER_CSV)

    if [str(h) for h in header] != [str(h) for h in older_header]:
        raise SystemExit(
            "The two halves of the compilation no longer share a layout. The "
            "export is positional, so it must not be written until they agree.")

    # The halves are meant to be disjoint at the cutoff and free of duplicates;
    # both are enforced here rather than assumed, since a silent overlap would
    # double-list a pole.
    recreated = {str(r[COL_ROCKNAME - 1]).strip() for r in our_rows}
    for row in carried:
        name = str(row[COL_ROCKNAME - 1]).strip()
        if name in recreated:
            raise SystemExit(f"{name!r} appears in both halves of the "
                             "compilation; regenerate with extract_older_poles.py")
        try:
            if float(row[COL_NOMINAL_AGE - 1]) <= cutoff:
                raise SystemExit(f"{name!r} is at or below the {cutoff:.0f} Ma "
                                 "cutoff but sits in the older-pole file")
        except (TypeError, ValueError):
            raise SystemExit(f"{name!r} in the older-pole file has no age")

    def age_of(row):
        try:
            return float(row[COL_NOMINAL_AGE - 1])
        except (TypeError, ValueError):
            return float("inf")

    rows = sorted(our_rows + carried, key=age_of)

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = SHEET_TITLE
    ws.append(header)
    for i, row in enumerate(rows):
        if len(row) != len(header):
            raise SystemExit(f"row {i + 2} has {len(row)} fields, "
                             f"expected {len(header)}")
        excel_row = i + 2
        values = list(row)
        # age diff is a formula in the reference layout, not a stored value.
        # Both halves store it as a number, so it is written as a formula here,
        # for the row it actually lands on.
        values[COL_AGE_DIFF - 1] = (
            f"=BI{excel_row}-BH{excel_row}"
            if values[COL_LOMAGAGE - 1] is not None
            and values[COL_HIMAGAGE - 1] is not None else None)
        ws.append(values)

    ws.freeze_panes = FREEZE_PANES
    wb.save(OUT_XLSX)
    print(f"Site-level poles (<= {cutoff:.0f} Ma): {len(our_rows)}")
    print(f"Older poles from {os.path.basename(OLDER_CSV)} (> {cutoff:.0f} Ma): "
          f"{len(carried)}")
    print(f"Wrote {os.path.relpath(OUT_XLSX, ROOT)} "
          f"({len(rows)} poles x {len(header)} columns)")


if __name__ == "__main__":
    main()
