from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Takiyuak Formation',
    result_name='Takiyuak Formation site means',
    citations='Irving1979Takiyuak',
    geologic_classes='Sedimentary',
    geologic_types='Formation',
    lithologies='Sandstone:Siltstone',
    method_codes='DE-FM',
    description='Site-mean directions from Irving and McGlynn (1979).',
    age='1876',
    age_low='1866',
    age_high='1886',
    sources=[
        CsvInput(
            filename='Irving1979_Takiyuak.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a95',
            k_col='k',
            n_col='N',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)