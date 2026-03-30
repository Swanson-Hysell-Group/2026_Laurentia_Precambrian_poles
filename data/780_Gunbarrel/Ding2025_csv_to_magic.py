"""
Convert Ding2025a.csv to MagIC 3.0 sites.txt and locations.txt, then combine
into a MagIC upload file using ipmag.upload_magic().

Source: Ding, J., Evans, D. A. D., Kilian, T., Mitchell, R. N., Gong, Z.,
    Chamberlain, K., et al. (2025). Local rotations biasing supercontinental
    configurations: Revisiting a Key ca. 780 Ma paleomagnetic pole for
    Laurentia and Rodinia. Journal of Geophysical Research: Solid Earth, 130,
    e2025JB031762. https://doi.org/10.1029/2025JB031762

MagIC data model: v3.0 (https://www2.earthref.org/MagIC/data-models/3.0)

Each site produces two rows in sites.txt:
    - dir_tilt_correction = 0: geographic coordinates (Dg, Ig)
    - dir_tilt_correction = 100: stratigraphic coordinates (Ds, Is) with VGPs

Notes on metadata choices:
    - Wyoming craton sites (BT50-BT72, J20B24, BT-8, MM1) are mafic dikes;
      "fine- to medium-grained diabase consisting of clinopyroxene,
      orthopyroxene, plagioclase, with minor biotite, hornblende, and
      interstitial granophyre" (Ding et al., 2025, Section 2).
    - Slave craton sites (Gun, Mar) are Tsezotene sills from the Mackenzie
      Mountains (Park, Buchan & Gandhi, 1995).
    - Cal_mean is the Calder sheet at Hottah Lake; "gabbro-diabase, composed
      of plagioclase and pyroxene" (Ding et al., 2025, Section 2).
    - New Ding et al. sites used stepwise thermal demagnetization (LP-DIR-T),
      PCA component fitting (DE-BFL), and Fisher mean (DE-FM).
    - Legacy sites (BT-8, MM1 from Harlan et al. 1997; Gun, Mar from Park,
      Buchan & Gandhi 1995) only have DE-FM since lab methods are not detailed.
    - Mar has no original DOI (1995 pre-DOI publication); directions are
      reported in Ding et al. (2025), so that DOI is used.
    - Longitude converted from degrees West to MagIC 0-360 East convention.
    - Citation semicolons converted to MagIC colon-delimited format.
    - Footnote markers (non-breaking space + "a") stripped from site names.
"""

import csv
from pathlib import Path
import pmagpy.ipmag as ipmag

# --- Configuration -----------------------------------------------------------

LOCATION = 'Gunbarrel LIP'

# Default age bounds for undated Gunbarrel LIP sites (Ma)
LIP_AGE_LOW = '775'
LIP_AGE_HIGH = '780'

MAGIC_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_sigma', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95',
    'dir_n_samples', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description'
]

# Per-site metadata from Ding et al. (2025) and references therein.
# Keys are base site names (before any footnote suffixes).
SITE_META = {
    # New Ding et al. (2025) Wyoming craton dikes - Beartooth Mountains
    'BT50':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'BT51':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'BT54':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'BT5658':   {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'BT64':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'BT68':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'BT69':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'BT72':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
    'J20B24':   {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},

    # Harlan et al. (1997) Wyoming craton dikes
    'BT-8':     {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'DE-FM'},
    'MM1':      {'geologic_types': 'Volcanic Dike', 'lithologies': 'Diabase',
                 'method_codes': 'DE-FM'},

    # Slave craton - Tsezotene sills, Mackenzie Mountains
    'Gun':      {'geologic_types': 'Sill', 'lithologies': 'Diabase',
                 'method_codes': 'DE-FM'},
    'Mar':      {'geologic_types': 'Sill', 'lithologies': 'Diabase',
                 'method_codes': 'DE-FM',
                 'citation_override': '10.1029/2025JB031762'},

    # Slave craton - Calder sheet, Hottah Lake (this study + prior)
    'Cal_mean': {'geologic_types': 'Sill', 'lithologies': 'Diabase',
                 'method_codes': 'LP-DIR-T:DE-BFL:DE-FM'},
}

# --- Conversion --------------------------------------------------------------


def clean_site_name(raw_name):
    """Strip footnote markers (e.g. non-breaking space + 'a') from site IDs."""
    return raw_name.replace('\xa0', ' ').split()[0]


def convert_lon_west_to_east(lon_w):
    """Convert longitude in degrees West to MagIC 0-360 East convention."""
    return round(360 - float(lon_w), 3)


def format_citations(raw_citations):
    """Convert semicolon-delimited citations to MagIC colon-delimited format."""
    if not raw_citations:
        return ''
    return ':'.join(c.strip() for c in raw_citations.split(';'))


def main():
    script_dir = Path(__file__).parent
    csv_path = script_dir / 'Ding2025a.csv'
    out_path = script_dir / 'sites.txt'

    rows_out = []

    with open(csv_path, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            site = clean_site_name(row['ID'].strip())
            meta = SITE_META.get(site, {})

            lat = row['Slat (°N)'].strip()
            lon = str(convert_lon_west_to_east(row['Slong (°W)'].strip()))

            n = row['n/N'].strip()
            dg = row['Dg (°)'].strip()
            ig = row['Ig (°)'].strip()
            ds = row['Ds (°)'].strip()
            is_ = row['Is (°)'].strip()

            k = row['k'].strip()
            if k in ('\u2013', ''):  # en-dash or empty
                k = ''
            a95 = row['\u03b195(\u00b0)'].strip()  # α95(°)

            plat = row['Plat (°N)'].strip()
            plon = row['Plon (°E)'].strip()
            dp = row['vgp_dp'].strip()
            dm = row['vgp_dm'].strip()

            age = row['age'].strip()
            age_sigma = row['age_sigma'].strip()
            age_unit = row['age_unit'].strip() if row['age_unit'].strip() else 'Ma'

            # Sites without specific ages get LIP-wide age bounds
            if age and age_sigma:
                age_low = ''
                age_high = ''
            else:
                age_low = LIP_AGE_LOW
                age_high = LIP_AGE_HIGH

            if 'citation_override' in meta:
                citations = meta['citation_override']
            else:
                citations = format_citations(row['citations'].strip())

            comments = row['comments'].strip()

            common = {
                'site': site,
                'location': LOCATION,
                'result_type': 'i',
                'result_quality': 'g',
                'method_codes': meta.get('method_codes', 'DE-FM'),
                'citations': citations,
                'geologic_classes': 'Igneous',
                'geologic_types': meta.get('geologic_types', 'Volcanic Dike'),
                'lithologies': meta.get('lithologies', 'Diabase'),
                'lat': lat,
                'lon': lon,
                'age': age,
                'age_sigma': age_sigma,
                'age_low': age_low,
                'age_high': age_high,
                'age_unit': age_unit,
                'dir_k': k,
                'dir_alpha95': a95,
                'dir_n_samples': n,
                'description': comments,
            }

            # Geographic coordinates row (tilt_correction = 0)
            geo = dict(common)
            geo['dir_tilt_correction'] = '0'
            geo['dir_dec'] = dg
            geo['dir_inc'] = ig
            geo['vgp_lat'] = ''
            geo['vgp_lon'] = ''
            geo['vgp_dp'] = ''
            geo['vgp_dm'] = ''
            rows_out.append(geo)

            # Stratigraphic coordinates row (tilt_correction = 100)
            strat = dict(common)
            strat['dir_tilt_correction'] = '100'
            strat['dir_dec'] = ds
            strat['dir_inc'] = is_
            strat['vgp_lat'] = plat
            strat['vgp_lon'] = plon
            strat['vgp_dp'] = dp
            strat['vgp_dm'] = dm
            rows_out.append(strat)

    with open(out_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(MAGIC_COLS) + '\n')
        for r in rows_out:
            f.write('\t'.join(r[col] for col in MAGIC_COLS) + '\n')

    print(f'Wrote {len(rows_out)} rows ({len(rows_out) // 2} sites '
          f'x 2 tilt corrections) to {out_path.name}')

    # --- Generate locations.txt from site-level data -------------------------
    write_locations(script_dir, rows_out)


def write_locations(script_dir, site_rows):
    """Calculate mean pole from site VGPs and write MagIC locations.txt.

    Args:
        script_dir: Path to the output directory.
        site_rows: List of site row dicts (both geographic and stratigraphic).
    """
    # Filter to stratigraphic rows only
    strat_rows = [r for r in site_rows if r['dir_tilt_correction'] == '100']

    # Fisher mean of VGPs
    vgp_lons = [float(r['vgp_lon']) for r in strat_rows if r['vgp_lon']]
    vgp_lats = [float(r['vgp_lat']) for r in strat_rows if r['vgp_lat']]
    pole_mean = ipmag.fisher_mean(dec=vgp_lons, inc=vgp_lats)

    # Collect metadata across all sites
    all_sites = [r['site'] for r in strat_rows]
    unique_citations = list(dict.fromkeys(
        r['citations'] for r in strat_rows if r['citations']))
    all_citations = ':'.join(unique_citations)
    unique_lithologies = list(dict.fromkeys(
        r['lithologies'] for r in strat_rows if r['lithologies']))
    all_lithologies = ':'.join(unique_lithologies)
    all_method_codes = ':'.join(sorted(set(
        code for r in strat_rows
        for code in r['method_codes'].split(':') if r['method_codes'])))

    # Geographic bounds
    lats = [float(r['lat']) for r in strat_rows]
    lons = [float(r['lon']) for r in strat_rows]

    # Nominal pole age and bounds (Ma)
    age_mid = 780
    age_low = 778
    age_high = 782

    loc_cols = [
        'location', 'location_type', 'result_name', 'result_type',
        'result_quality', 'method_codes', 'citations', 'geologic_classes',
        'lithologies', 'lat_s', 'lat_n', 'lon_w', 'lon_e',
        'age', 'age_low', 'age_high', 'age_unit',
        'dir_tilt_correction',
        'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k', 'pole_n_sites',
        'sites', 'description'
    ]

    loc_row = {
        'location': LOCATION,
        'location_type': 'Region',
        'result_name': 'Gunbarrel LIP ca. 780 Ma pole',
        'result_type': 'a',
        'result_quality': 'g',
        'method_codes': all_method_codes + ':DE-VGP',
        'citations': all_citations,
        'geologic_classes': 'Igneous',
        'lithologies': all_lithologies,
        'lat_s': f'{min(lats):.3f}',
        'lat_n': f'{max(lats):.3f}',
        'lon_w': f'{min(lons):.3f}',
        'lon_e': f'{max(lons):.3f}',
        'age': f'{age_mid:.1f}',
        'age_low': f'{age_low:.1f}',
        'age_high': f'{age_high:.1f}',
        'age_unit': 'Ma',
        'dir_tilt_correction': '100',
        'pole_lat': f'{pole_mean["inc"]:.1f}',
        'pole_lon': f'{pole_mean["dec"]:.1f}',
        'pole_alpha95': f'{pole_mean["alpha95"]:.1f}',
        'pole_k': f'{pole_mean["k"]:.1f}',
        'pole_n_sites': str(pole_mean['n']),
        'sites': ':'.join(all_sites),
        'description': 'Mean pole from 14 Gunbarrel LIP intrusions '
                        '(Ding et al. 2025)',
    }

    loc_path = script_dir / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(loc_cols) + '\n')
        f.write('\t'.join(loc_row[col] for col in loc_cols) + '\n')

    print(f'Wrote {loc_path.name}: pole at {loc_row["pole_lat"]}°N, '
          f'{loc_row["pole_lon"]}°E, A95={loc_row["pole_alpha95"]}°, '
          f'N={loc_row["pole_n_sites"]}')


def combine_and_validate(script_dir):
    """Combine MagIC tables into upload.txt and validate.

    Uses ipmag.upload_magic() to find all MagIC tables (sites.txt,
    locations.txt) in the directory, combine them into a single upload.txt,
    and validate against the MagIC data model.

    Args:
        script_dir: Path to the directory containing MagIC table files.
    """
    result = ipmag.upload_magic(dir_path=str(script_dir),
                                input_dir_path=str(script_dir))
    if result[0] is False:
        print(f'Validation issue: {result[1]}')
    else:
        print(f'Created {result[0]}')


if __name__ == '__main__':
    main()
    combine_and_validate(Path(__file__).parent)
