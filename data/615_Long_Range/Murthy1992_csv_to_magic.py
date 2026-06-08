from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Long Range Dykes',
    result_name='Long Range Dykes site means',
    citations='Murthy1992',
    geologic_classes='Igneous',
    geologic_types='Volcanic Dike',
    lithologies='Diabase',
    method_codes='DE-FM',
    description='Site-mean directions from Murthy et al. (1992).',
    age='615',
    age_low='613',
    age_high='617',
    sources=[
        CsvInput(
            filename='Murthy1992.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a95',
            k_col='k',
            n_col='N',
            lat_col='site_lat',
            lon_col='site_lon',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)