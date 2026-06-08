from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag


SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_sigma', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95',
    'dir_n_samples', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description',
]

LOCATION_COLS = [
    'location', 'location_type', 'result_name', 'result_type', 'result_quality',
    'method_codes', 'citations', 'geologic_classes', 'geologic_types',
    'lithologies', 'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age', 'age_sigma',
    'age_low', 'age_high', 'age_unit', 'pole_lat', 'pole_lon', 'pole_alpha95',
    'pole_k', 'pole_n_sites', 'description',
]


@dataclass
class CsvInput:
    filename: str
    site_col: str
    dec_col: str | None = None
    inc_col: str | None = None
    dec_geo_col: str | None = None
    inc_geo_col: str | None = None
    dec_tc_col: str | None = None
    inc_tc_col: str | None = None
    k_col: str | None = None
    a95_col: str | None = None
    n_col: str | None = None
    lat_col: str | None = None
    lon_col: str | None = None
    vgp_lat_col: str | None = None
    vgp_lon_col: str | None = None
    vgp_dp_col: str | None = None
    vgp_dm_col: str | None = None
    description_col: str | None = None
    citations_col: str | None = None
    default_tilt_correction: int = 0
    vgp_tilt_correction: int | None = None
    read_csv_kwargs: dict[str, Any] = field(default_factory=dict)
    column_names: list[str] | None = None
    extra_cols: dict[str, str] = field(default_factory=dict)


@dataclass
class PoleMagicConfig:
    location: str
    result_name: str
    citations: str = ''
    geologic_classes: str = ''
    geologic_types: str = ''
    lithologies: str = ''
    method_codes: str = 'DE-FM'
    description: str = ''
    age: str = ''
    age_sigma: str = ''
    age_low: str = ''
    age_high: str = ''
    age_unit: str = 'Ma'
    result_quality: str = 'g'
    result_type: str = 'i'
    location_type: str = 'i'
    sources: list[CsvInput] = field(default_factory=list)


def _value(row: pd.Series, column: str | None) -> Any:
    if column is None or column not in row.index:
        return None
    return row[column]


def _stringify(value: Any) -> str:
    if value is None or pd.isna(value):
        return ''
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _as_float(value: Any) -> float | None:
    if value is None or pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _lon_to_360(value: Any) -> str:
    lon = _as_float(value)
    if lon is None:
        return ''
    return f'{lon % 360:.3f}'


def _signed_lon(value: str) -> float | None:
    lon = _as_float(value)
    if lon is None:
        return None
    return ((lon + 180) % 360) - 180


def _direction_sets(spec: CsvInput) -> list[tuple[int, str | None, str | None]]:
    sets: list[tuple[int, str | None, str | None]] = []
    if spec.dec_geo_col or spec.inc_geo_col:
        sets.append((0, spec.dec_geo_col, spec.inc_geo_col))
    if spec.dec_tc_col or spec.inc_tc_col:
        sets.append((100, spec.dec_tc_col, spec.inc_tc_col))
    if not sets:
        sets.append((spec.default_tilt_correction, spec.dec_col, spec.inc_col))
    return sets


def _compute_vgp(dec: Any, inc: Any, a95: Any, site_lat: Any, site_lon: Any) -> tuple[str, str, str, str]:
    dec_f = _as_float(dec)
    inc_f = _as_float(inc)
    lat_f = _as_float(site_lat)
    lon_f = _as_float(site_lon)
    if None in (dec_f, inc_f, lat_f, lon_f):
        return '', '', '', ''

    alpha95 = _as_float(a95)
    pole_lon, pole_lat, dp, dm = pmag.dia_vgp(dec_f, inc_f, alpha95 or 0.0, lat_f, lon_f)
    if alpha95 is None:
        return f'{pole_lat:.3f}', f'{pole_lon % 360:.3f}', '', ''
    return f'{pole_lat:.3f}', f'{pole_lon % 360:.3f}', f'{dp:.3f}', f'{dm:.3f}'


def _build_site_rows(config: PoleMagicConfig, script_dir: Path) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    for spec in config.sources:
        df = pd.read_csv(script_dir / spec.filename, **spec.read_csv_kwargs)
        if spec.column_names is not None:
            df.columns = spec.column_names

        for _, row in df.iterrows():
            site_name = _stringify(_value(row, spec.site_col))
            if not site_name:
                continue

            for tilt_code, dec_col, inc_col in _direction_sets(spec):
                dec_val = _value(row, dec_col)
                inc_val = _value(row, inc_col)
                has_direction = pd.notna(dec_val) or pd.notna(inc_val)

                use_vgp = spec.vgp_tilt_correction is None or spec.vgp_tilt_correction == tilt_code
                vgp_lat = _stringify(_value(row, spec.vgp_lat_col)) if use_vgp else ''
                vgp_lon = _lon_to_360(_value(row, spec.vgp_lon_col)) if use_vgp else ''
                vgp_dp = _stringify(_value(row, spec.vgp_dp_col)) if use_vgp else ''
                vgp_dm = _stringify(_value(row, spec.vgp_dm_col)) if use_vgp else ''

                if use_vgp and not (vgp_lat and vgp_lon):
                    vgp_lat, vgp_lon, vgp_dp, vgp_dm = _compute_vgp(
                        dec_val,
                        inc_val,
                        _value(row, spec.a95_col),
                        _value(row, spec.lat_col),
                        _value(row, spec.lon_col),
                    )

                if not has_direction and not (vgp_lat and vgp_lon):
                    continue

                rows_out.append({
                    'site': site_name,
                    'location': config.location,
                    'result_type': config.result_type,
                    'result_quality': config.result_quality,
                    'method_codes': config.method_codes,
                    'citations': _stringify(_value(row, spec.citations_col)) or config.citations,
                    'geologic_classes': config.geologic_classes,
                    'geologic_types': config.geologic_types,
                    'lithologies': config.lithologies,
                    'lat': _stringify(_value(row, spec.lat_col)),
                    'lon': _lon_to_360(_value(row, spec.lon_col)),
                    'age': _stringify(config.age),
                    'age_sigma': _stringify(config.age_sigma),
                    'age_low': _stringify(config.age_low),
                    'age_high': _stringify(config.age_high),
                    'age_unit': config.age_unit,
                    'dir_tilt_correction': str(tilt_code),
                    'dir_dec': _stringify(dec_val),
                    'dir_inc': _stringify(inc_val),
                    'dir_k': _stringify(_value(row, spec.k_col)),
                    'dir_alpha95': _stringify(_value(row, spec.a95_col)),
                    'dir_n_samples': _stringify(_value(row, spec.n_col)),
                    'vgp_lat': vgp_lat,
                    'vgp_lon': vgp_lon,
                    'vgp_dp': vgp_dp,
                    'vgp_dm': vgp_dm,
                    'description': _stringify(_value(row, spec.description_col)) or config.description,
                })
                for out_col, source_col in spec.extra_cols.items():
                    rows_out[-1][out_col] = _stringify(_value(row, source_col))
    return rows_out


def _build_location_row(config: PoleMagicConfig, site_rows: list[dict[str, str]]) -> dict[str, str]:
    lats = [float(row['lat']) for row in site_rows if row['lat']]
    signed_lons = [_signed_lon(row['lon']) for row in site_rows if row['lon']]
    signed_lons = [lon for lon in signed_lons if lon is not None]

    vgp_rows = [row for row in site_rows if row['vgp_lat'] and row['vgp_lon']]
    if vgp_rows:
        pole_mean = ipmag.fisher_mean(
            dec=[float(row['vgp_lon']) for row in vgp_rows],
            inc=[float(row['vgp_lat']) for row in vgp_rows],
        )
        pole_lat = f"{pole_mean['inc']:.3f}"
        pole_lon = f"{pole_mean['dec'] % 360:.3f}"
        pole_alpha95 = f"{pole_mean['alpha95']:.3f}"
        pole_k = f"{pole_mean['k']:.3f}"
        pole_n = str(pole_mean['n'])
    else:
        pole_lat = pole_lon = pole_alpha95 = pole_k = pole_n = ''

    return {
        'location': config.location,
        'location_type': config.location_type,
        'result_name': config.result_name,
        'result_type': config.result_type,
        'result_quality': config.result_quality,
        'method_codes': config.method_codes,
        'citations': config.citations,
        'geologic_classes': config.geologic_classes,
        'geologic_types': config.geologic_types,
        'lithologies': config.lithologies,
        'lat_s': f'{min(lats):.3f}' if lats else '',
        'lat_n': f'{max(lats):.3f}' if lats else '',
        'lon_w': f'{min(signed_lons):.3f}' if signed_lons else '',
        'lon_e': f'{max(signed_lons):.3f}' if signed_lons else '',
        'age': _stringify(config.age),
        'age_sigma': _stringify(config.age_sigma),
        'age_low': _stringify(config.age_low),
        'age_high': _stringify(config.age_high),
        'age_unit': config.age_unit,
        'pole_lat': pole_lat,
        'pole_lon': pole_lon,
        'pole_alpha95': pole_alpha95,
        'pole_k': pole_k,
        'pole_n_sites': pole_n,
        'description': config.description,
    }


def _write_table(path: Path, table_name: str, columns: list[str], rows: list[dict[str, str]]) -> None:
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write(f'tab delimited\t{table_name}\n')
        handle.write('\t'.join(columns) + '\n')
        for row in rows:
            handle.write('\t'.join(row.get(column, '') for column in columns) + '\n')


def convert_csvs_to_magic(config: PoleMagicConfig, script_dir: Path) -> tuple[Path, Path]:
    site_rows = _build_site_rows(config, script_dir)
    if not site_rows:
        raise ValueError(f'No rows were converted for {config.location}')

    sites_path = script_dir / 'sites.txt'
    locations_path = script_dir / 'locations.txt'
    extra_cols = [col for col in site_rows[0].keys() if col not in SITE_COLS]
    _write_table(sites_path, 'sites', SITE_COLS + extra_cols, site_rows)
    _write_table(locations_path, 'locations', LOCATION_COLS, [_build_location_row(config, site_rows)])
    print(f'Wrote {len(site_rows)} site rows to {sites_path}')
    print(f'Wrote 1 location row to {locations_path}')
    return sites_path, locations_path