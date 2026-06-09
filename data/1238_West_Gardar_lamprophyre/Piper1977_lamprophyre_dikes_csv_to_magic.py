from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='West Gardar Lamprophyre Dykes',
    result_name='West Gardar Lamprophyre Dykes site means',
    citations='Piper1977Lamprophyre',
    geologic_classes='Igneous',
    geologic_types='Volcanic Dike',
    lithologies='Lamprophyre',
    method_codes='DE-FM',
    description='Site-mean directions from Piper and Stearn (1977).',
    sources=[
        CsvInput(
            filename='Piper1977_lamprophyre_dikes.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a95',
            k_col='k',
            n_col='n',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)