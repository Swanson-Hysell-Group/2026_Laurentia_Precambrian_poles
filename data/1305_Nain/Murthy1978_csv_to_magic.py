from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Nain Anorthosite',
    result_name='Nain Anorthosite site means',
    citations='Murthy1978',
    geologic_classes='Igneous',
    geologic_types='Pluton',
    lithologies='Anorthosite',
    method_codes='DE-FM',
    description='Site-mean directions from Murthy (1978).',
    age='1305',
    age_low='1283',
    age_high='1327',
    sources=[
        CsvInput(
            filename='Murthy1978.csv',
            site_col='Unnamed: 0',
            dec_col='dec',
            inc_col='inc',
            n_col='Unnamed: 1',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)