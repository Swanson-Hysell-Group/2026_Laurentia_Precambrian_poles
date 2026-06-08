from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Callander Alkaline Complex',
    result_name='Callander Alkaline Complex site means',
    citations='10.1139/e91-033',
    geologic_classes='Igneous',
    geologic_types='Pluton',
    lithologies='Alkaline Rock',
    method_codes='DE-FM',
    description='Site-mean directions from Symons and Chiasson (1991).',
    age='575',
    age_low='570',
    age_high='580',
    sources=[
        CsvInput(
            filename='Symons1991.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a_95',
            k_col='k',
            description_col='rock_type',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)