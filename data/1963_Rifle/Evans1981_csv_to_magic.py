from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Rifle Formation',
    result_name='Rifle Formation site means',
    citations='Evans1981',
    geologic_classes='Sedimentary',
    geologic_types='Formation',
    lithologies='Sandstone',
    method_codes='DE-FM',
    description='Site-mean directions from Evans and Hoye (1981).',
    age='1963',
    age_low='1957',
    age_high='1969',
    sources=[
        CsvInput(
            filename='Evans1981.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)