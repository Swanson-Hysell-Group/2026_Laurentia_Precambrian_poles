from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Defeat Suite',
    result_name='Defeat Suite site means',
    citations='Mitchell2014',
    geologic_classes='Igneous',
    geologic_types='Pluton',
    lithologies='Granitoid',
    method_codes='DE-FM',
    description='Site-mean directions from Mitchell et al. (2014).',
    age='2625',
    age_low='2620',
    age_high='2630',
    sources=[
        CsvInput(
            filename='Mitchell2014.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a95',
            k_col='k',
            n_col='n',
            lat_col='site_lat',
            lon_col='site_lon',
            vgp_lat_col='Plat',
            vgp_lon_col='Plon',
            vgp_tilt_correction=0,
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)