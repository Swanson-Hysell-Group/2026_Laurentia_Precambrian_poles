from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Rocky Mountain intrusions',
    result_name='Rocky Mountain intrusion site means',
    citations='Harlan1994:Harlan1998',
    geologic_classes='Igneous',
    geologic_types='Pluton',
    lithologies='Anorthosite:Syenite',
    method_codes='DE-FM',
    description='Site-mean in-situ and tilt-corrected directions from Harlan et al. (1994, 1998).',
    age='1430',
    age_low='1415',
    age_high='1445',
    sources=[
        CsvInput(
            filename='Harlan1994a.csv',
            site_col='site',
            dec_geo_col='dec',
            inc_geo_col='inc',
            dec_tc_col='dec_tc',
            inc_tc_col='inc_tc',
            a95_col='alpha95',
            k_col='k',
            n_col='N',
            lat_col='lat',
            lon_col='lon',
            vgp_lat_col='Plat',
            vgp_lon_col='Plon',
            vgp_tilt_correction=100,
            description_col='rock type',
            extra_cols={
                'strike': 'strike',
                'dip': 'dip',
                'n': 'n',
                'N': 'N',
            },
        ),
        CsvInput(
            filename='Harlan1998a.csv',
            site_col='site',
            dec_geo_col='dec',
            inc_geo_col='inc',
            dec_tc_col='dec_tc',
            inc_tc_col='inc_tc',
            a95_col='alpha95',
            k_col='k',
            n_col='N',
            lat_col='lat',
            lon_col='lon',
            vgp_lat_col='Plat',
            vgp_lon_col='Plon',
            vgp_tilt_correction=100,
            description_col='rock type',
            extra_cols={
                'n': 'n',
                'N': 'N',
            },
        ),
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)