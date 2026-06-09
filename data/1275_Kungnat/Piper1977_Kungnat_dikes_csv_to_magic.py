from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Kungnat Ring Dyke',
    result_name='Kungnat Ring Dyke site means',
    citations='Piper1977Kungnat',
    geologic_classes='Igneous',
    geologic_types='Volcanic Dike',
    lithologies='Diabase',
    method_codes='DE-FM',
    description='Site-mean directions from Piper and Stearn (1977).',
    age='1275',
    age_low='1273',
    age_high='1277',
    sources=[
        CsvInput(
            filename='Piper1977_Kungnat_dikes.csv',
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