from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Sparrow Dykes',
    result_name='Sparrow Dykes site means',
    citations='McGlynn1974',
    geologic_classes='Igneous',
    geologic_types='Volcanic Dike',
    lithologies='Diabase',
    method_codes='DE-FM',
    description='Site-mean directions from McGlynn et al. (1974).',
    age='1827',
    age_low='1823',
    age_high='1831',
    sources=[
        CsvInput(
            filename='McGlynn1974.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            k_col='k',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)