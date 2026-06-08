from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Franklin LIP',
    result_name='Franklin LIP legacy VGP compilation',
    citations='Denyszyn2009',
    geologic_classes='Igneous',
    geologic_types='Volcanic Dike:Sill:Flow',
    lithologies='Diabase:Basalt',
    method_codes='DE-FM',
    description='Legacy VGP compilation assembled for the Franklin grand mean.',
    age='719',
    age_low='718',
    age_high='720',
    sources=[
        CsvInput(
            filename='Denyszyn2009.csv',
            site_col='site',
            default_tilt_correction=100,
            vgp_tilt_correction=100,
            vgp_lat_col='vgp_lat_rotated',
            vgp_lon_col='vgp_lon_rotated',
            description_col='note',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)