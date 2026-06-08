from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Uinta Mountain Group',
    result_name='Uinta Mountain Group site means',
    citations='Weil2006',
    geologic_classes='Sedimentary',
    geologic_types='Formation',
    lithologies='Sandstone:Siltstone',
    method_codes='DE-FM',
    description='Site-mean in-situ and tilt-corrected directions from Weil et al. (2006).',
    age='750',
    age_low='717',
    age_high='771',
    sources=[
        CsvInput(
            filename='Weil2006.csv',
            site_col='site',
            dec_geo_col='dec_is',
            inc_geo_col='inc_is',
            dec_tc_col='dec_tc',
            inc_tc_col='inc_tc',
            a95_col='alpha95',
            n_col='n/no',
            description_col='section',
            extra_cols={
                'filtered_out_by_direction': 'filtered_out_by_direction',
                'strike': 'strike',
                'dip': 'dip',
                'section': 'section',
            },
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)