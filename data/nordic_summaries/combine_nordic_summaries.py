"""Combine per-notebook Nordic summary CSVs into a single table.

Each pole notebook writes a single-row CSV (named after the notebook) into
this folder via ``pole_tools.save_nordic_summary``. This script reads all of
those per-pole CSVs and concatenates them into one combined CSV
(``nordic_summaries_combined.csv``), with one row per pole sorted by the
source filename.

Usage:
    python combine_nordic_summaries.py
"""

import glob
import os

import pandas as pd

# directory containing this script (and the per-notebook summary CSVs)
SUMMARY_DIR = os.path.dirname(os.path.abspath(__file__))
COMBINED_FILENAME = 'nordic_summaries_combined.csv'


def combine_summaries(summary_dir=SUMMARY_DIR, combined_filename=COMBINED_FILENAME):
    """Combine all per-notebook summary CSVs in a directory into one CSV.

    Args:
        summary_dir (str): Directory holding the per-notebook summary CSVs.
        combined_filename (str): Name of the combined CSV to write into
            ``summary_dir``. This file is excluded from the inputs if present.

    Returns:
        pandas.DataFrame: The combined summary table (also written to disk).
    """
    csv_paths = sorted(glob.glob(os.path.join(summary_dir, '*.csv')))
    combined_path = os.path.join(summary_dir, combined_filename)
    # do not fold a previously written combined file back into the inputs
    csv_paths = [p for p in csv_paths if os.path.abspath(p) != os.path.abspath(combined_path)]

    if not csv_paths:
        raise FileNotFoundError(
            f'No per-notebook summary CSVs found in {summary_dir}. '
            'Run the pole notebooks to generate them first.'
        )

    frames = []
    for path in csv_paths:
        df = pd.read_csv(path)
        # record provenance so each row is traceable to its source notebook
        df.insert(0, 'source_file', os.path.splitext(os.path.basename(path))[0])
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(combined_path, index=False)
    print(f'Combined {len(csv_paths)} summaries into {combined_path}')
    return combined


if __name__ == '__main__':
    combine_summaries()
