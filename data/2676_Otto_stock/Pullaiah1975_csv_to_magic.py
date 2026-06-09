from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Otto Stock dike swarm',
    result_name='Otto Stock site means',
    citations='Pullaiah1975',
    geologic_classes='Igneous',
    geologic_types='Volcanic Dike',
    lithologies='Diabase',
    method_codes='DE-FM',
    description='Site-mean directions from Pullaiah et al. (1975).',
    age='2676',
    age_low='2672',
    age_high='2680',
    sources=[
        CsvInput(
            filename='Pullaiah1975.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a95',
            k_col='k',
            n_col='n5',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)