from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Boot-Phantom plutons',
    result_name='Boot-Phantom site means',
    citations='Symons1999',
    geologic_classes='Igneous',
    geologic_types='Pluton',
    lithologies='Granitoid',
    method_codes='DE-FM',
    description='Site-mean directions from Symons et al. (1999).',
    age='1838',
    age_low='1835',
    age_high='1841',
    sources=[
        CsvInput(
            filename='Symons1999.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a95',
            k_col='k',
            description_col='rock type',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)