"""Extract the pre-compilation Laurentia poles into a CSV in the Nordic layout.

The compilation is assembled from two files:

- ``data/nordic_summaries/nordic_summaries_combined.csv`` -- poles at or below
  ``PUBLISH_AGE_MAX``, recomputed from site-level data by this project and
  documented in a notebook apiece.
- ``data/older_Laurentia_poles.csv`` -- this file's output. Poles older than
  that cutoff, which have *not* been recreated at the site level here and are
  carried over from the Nordic workshop compilation unchanged.

Both carry the same 71 columns in the same order, so concatenating them gives
the full record without any reshaping.

The rows are taken from the workshop workbook rather than from this project's
own copy of the earlier compilation (``data/Laurentia_poles.csv``, a 68-column
layout) because the workbook is already in the target layout: no column mapping
step means no column mapping to get wrong. This script is therefore run when the
workshop compilation is updated, not on every build; its output is committed.

    python scripts/extract_older_poles.py [path/to/workshop_workbook.xlsx]
"""

import csv
import os
import re
import sys

import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SUMMARY_CSV = os.path.join(ROOT, "data", "nordic_summaries",
                           "nordic_summaries_combined.csv")
BUILD_POLE_MAP = os.path.join(HERE, "build_pole_map.py")
OUT_CSV = os.path.join(ROOT, "data", "older_Laurentia_poles.csv")

DEFAULT_WORKBOOK = os.path.expanduser(
    "~/Downloads/Iloranta-global-2026-06-24.xlsx")

# Sampling localities absent from the workshop workbook, ROCKNAME -> (lat, lon
# in degrees east). A pole without one is still tabulated but cannot be placed
# on the interactive map, so the coordinate is supplied here rather than left
# blank. Each entry records where the value comes from.
#
#   Elbow Creek dikes -- mean of the thirteen component-A dikes of Ding et al.
#   (2024, Nature Communications 15:10814, doi:10.1038/s41467-024-55117-w)
#   Supplementary Table S1, the same thirteen that make their pole. They span
#   45.378-45.474 deg N and 109.878-110.139 deg W in the Stillwater Complex of
#   the Beartooth Mountains, Montana, giving a mean of 45.42 deg N, 109.96 deg W
#   (250.04 deg E). As a check, the VGP of their tilt-corrected mean direction
#   (D = 145.8, I = 59.6) at this site is 1.2 N / 275.4 E against their published
#   2.0 N / 275.3 E -- the 0.8 deg offset expected between the VGP of the mean
#   direction and the mean of the thirteen site VGPs that the pole actually is.
SITE_COORD_OVERRIDES = {
    "Elbow Creek dikes": (45.42, 250.04),
}

# 1-based column indices in the fixed layout.
COL_TERRANE = 1
COL_ROCKNAME = 2
COL_AGE_DIFF = 59
COL_NOMINAL_AGE = 58
COL_LOMAGAGE = 60
COL_HIMAGAGE = 61
COL_SLAT = 8
COL_SLONG = 9


def publish_age_max():
    """The age at which the site-level compilation hands over to older poles.

    Read from ``build_pole_map.py`` rather than duplicated, so the two cannot
    drift apart. That module is not imported because it pulls in the plotting
    stack for a single constant.

    Returns:
        float: The cutoff in Ma.
    """
    with open(BUILD_POLE_MAP, encoding="utf-8") as fh:
        m = re.search(r"^PUBLISH_AGE_MAX\s*=\s*([0-9.]+)", fh.read(), re.M)
    if not m:
        raise SystemExit("PUBLISH_AGE_MAX not found in build_pole_map.py")
    return float(m.group(1))


def cell(value):
    """Render a spreadsheet value for CSV, leaving blanks genuinely empty."""
    return "" if value is None else str(value)


def main():
    workbook = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_WORKBOOK
    if not (os.path.exists(workbook) and os.path.getsize(workbook) > 0):
        raise SystemExit(
            f"Workshop workbook not readable at {workbook}. Pass a path as the "
            "first argument, and if it is a cloud placeholder make it available "
            "offline first.")

    ws = openpyxl.load_workbook(workbook).active
    ncol = ws.max_column
    header = [ws.cell(row=1, column=c).value for c in range(1, ncol + 1)]

    with open(SUMMARY_CSV, newline="", encoding="utf-8-sig") as fh:
        our_header = next(csv.reader(fh))
        recreated = {r[COL_ROCKNAME - 1].strip()
                     for r in csv.reader(fh) if len(r) > 1}
    if [str(h) for h in header] != [str(h) for h in our_header]:
        raise SystemExit(
            "The workbook layout does not match the site-level summaries. The "
            "two halves of the compilation must share a layout to concatenate.")

    cutoff = publish_age_max()
    kept, superseded, filled = [], [], []
    for r in range(2, ws.max_row + 1):
        row = [ws.cell(row=r, column=c).value for c in range(1, ncol + 1)]
        if not str(row[COL_TERRANE - 1]).startswith("Laurentia"):
            continue
        try:
            age = float(row[COL_NOMINAL_AGE - 1])
        except (TypeError, ValueError):
            continue
        if age <= cutoff:
            continue
        name = str(row[COL_ROCKNAME - 1]).strip()
        if name in recreated:
            superseded.append(name)
            continue
        override = SITE_COORD_OVERRIDES.get(name)
        if override is not None:
            if row[COL_SLAT - 1] not in (None, "") \
                    and row[COL_SLONG - 1] not in (None, ""):
                raise SystemExit(
                    f"{name!r} now carries a sampling locality in the workbook; "
                    "drop its SITE_COORD_OVERRIDES entry rather than overwriting "
                    "the workbook value.")
            row[COL_SLAT - 1], row[COL_SLONG - 1] = override
            filled.append(name)
        # 'age diff' is a formula in the workbook, pointing at its row there.
        # Store the value instead, so the CSV stands on its own; the workbook
        # export rewrites it as a formula for the row it lands on.
        try:
            row[COL_AGE_DIFF - 1] = (int(round(float(row[COL_HIMAGAGE - 1])))
                                     - int(round(float(row[COL_LOMAGAGE - 1]))))
        except (TypeError, ValueError):
            row[COL_AGE_DIFF - 1] = None
        kept.append((age, row))

    kept.sort(key=lambda pair: pair[0])
    with open(OUT_CSV, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows([cell(v) for v in row] for _, row in kept)

    print(f"Source: {os.path.basename(workbook)}")
    if superseded:
        print("  dropped as superseded by a site-level recreation: "
              + ", ".join(sorted(superseded)))
    if filled:
        print("  sampling locality supplied from the literature: "
              + ", ".join(sorted(filled)))
    print(f"Wrote {os.path.relpath(OUT_CSV, ROOT)} "
          f"({len(kept)} poles > {cutoff:.0f} Ma x {len(header)} columns)")


if __name__ == "__main__":
    main()
