from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[2] / 'scripts'))

from legacy_csv_to_magic import CsvInput, PoleMagicConfig, convert_csvs_to_magic


CONFIG = PoleMagicConfig(
    location='Narssaq Gabbro',
    result_name='Narssaq Gabbro site means',
    citations='Piper1977Narssaq',
    geologic_classes='Igneous',
    geologic_types='Pluton',
    lithologies='Gabbro',
    method_codes='DE-FM',
    description='Site-mean directions from Piper (1977).',
    age='1184',
    age_low='1179',
    age_high='1189',
    sources=[
        CsvInput(
            filename='Piper1977_Narssaq.csv',
            site_col='site',
            dec_col='dec',
            inc_col='inc',
            a95_col='a95',
            k_col='k',
            n_col='N',
            vgp_lat_col='Plat',
            vgp_lon_col='Plon',
        )
    ],
)


if __name__ == '__main__':
    convert_csvs_to_magic(CONFIG, Path(__file__).parent)