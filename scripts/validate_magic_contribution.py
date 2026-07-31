"""Validate a combined MagIC upload file before submission.

A combined upload file is the ``.txt`` with the individual tables separated by
``>>>>>>>>>>`` on their own line (what the build scripts write and what is
uploaded to MagIC). This runs two layers of checks:

1. **PmagPy data-model validation** (``pmagpy.validate_upload3``) on every table:
   column names, controlled vocabularies, value types, and the data-model's
   conditional-requirement rules. Empty numeric cells are filled with ``''``
   first, because ``validate_upload3.test_type`` raises on ``NaN`` for
   integer-typed columns.

2. **Server-side rules the offline data model does NOT encode** — so PmagPy
   passes them but the MagIC online uploader rejects them. Currently:

   - **age completeness:** every row that has an ``age`` must also carry either
     ``age_sigma`` or *both* ``age_low`` and ``age_high``. The data-model
     validations only require that *some* age representation exists (``age`` or
     ``age_low``/``age_high``) plus ``age_unit``; they never require an
     uncertainty, so a bare ``age`` slips past PmagPy but is rejected on upload.

The two validators disagree in *both* directions, so neither is a perfect proxy
for the MagIC online uploader (which is authoritative): the online uploader is
stricter on the age rule above, but more lenient elsewhere — e.g. it accepts the
``standard='300'`` measurement values in the published ECMB contribution that
``validate_upload3`` flags as a controlled-vocabulary error. Treat the PmagPy
results as advisory. Pass ``tables=[...]`` to restrict validation to the tables a
build actually authors (skip pass-through tables copied verbatim from an
already-published, already-online-validated contribution) to avoid such false
positives.

Use programmatically (e.g. at the end of a build script)::

    from validate_magic_contribution import validate_upload_file
    if not validate_upload_file(path):
        raise SystemExit('MagIC validation failed')

or as a CLI::

    python validate_magic_contribution.py <upload_file.txt>

Exit status is 0 if the contribution passes, 1 otherwise.
"""

import os
import shutil
import sys
import tempfile
import warnings

SEP = '>>>>>>>>>>'


def _split_to_dir(path):
    """Write each table block of a combined upload file to its own file."""
    work = tempfile.mkdtemp(prefix='magic_val_')
    with open(path) as f:
        content = f.read()
    for chunk in content.split(SEP):
        chunk = chunk.strip('\n')
        if not chunk:
            continue
        marker, _, body = chunk.partition('\n')
        table = marker.split('\t')[-1].strip()
        with open(os.path.join(work, f'{table}.txt'), 'w') as out:
            out.write(f'tab delimited\t{table}\n')
            out.write(body if body.endswith('\n') else body + '\n')
    return work


def _bare_age_rows(df):
    """Rows with an ``age`` but neither ``age_sigma`` nor age_low+age_high."""
    if 'age' not in df.columns:
        return []

    def val(row, col):
        return str(row[col]).strip() if col in df.columns else ''

    bad = []
    for idx, row in df.iterrows():
        if val(row, 'age') == '':
            continue
        has_sigma = val(row, 'age_sigma') != ''
        has_range = val(row, 'age_low') != '' and val(row, 'age_high') != ''
        if not (has_sigma or has_range):
            bad.append(val(row, 'site') or val(row, 'location') or f'row {idx}')
    return bad


def validate_upload_file(path, tables=None, verbose=False):
    """Validate a combined MagIC upload file. Return True if it passes.

    Args:
        path: combined upload file (tables separated by ``>>>>>>>>>>``).
        tables: optional list of table types to validate (e.g.
            ``['locations', 'sites']``). Use it to restrict validation to the
            tables a build authors and skip pass-through tables copied verbatim
            from an already-published contribution. ``None`` validates all.
        verbose: pass through to ``validate_upload3`` for per-row detail.

    Returns:
        True if the age-completeness check passes on every checked table. The
        PmagPy data-model result is reported but advisory (it both misses and
        over-flags relative to the MagIC online uploader), so it does not by
        itself flip the return value.
    """
    warnings.filterwarnings('ignore')
    from pmagpy import contribution_builder as cb
    from pmagpy import validate_upload3 as v

    work = _split_to_dir(path)
    ok = True
    try:
        con = cb.Contribution(work)
        dtypes = [t for t in con.tables if tables is None or t in tables]
        for dtype in dtypes:
            table = con.tables[dtype]
            # avoid the NaN->int crash in validate_upload3.test_type
            table.df = table.df.fillna('')
            failed = v.validate_table(con, dtype, verbose=verbose,
                                      output_dir=work)
            if failed:
                print(f'-W- {dtype}: PmagPy data-model validation reported '
                      f'issues (advisory; verify against the MagIC uploader)')
            # server-side rule PmagPy does not enforce -> treat as a real error
            bare = _bare_age_rows(table.df)
            if bare:
                ok = False
                print(f'-E- {dtype}: {len(bare)} row(s) have an age with no '
                      f'age_sigma and no age_low/age_high (MagIC rejects a bare '
                      f'age): {bare}')
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print('-I- age-completeness check PASSED' if ok
          else '-E- age-completeness check FAILED')
    return ok


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('usage: python validate_magic_contribution.py <upload_file.txt>')
    sys.exit(0 if validate_upload_file(sys.argv[1], verbose=True) else 1)
