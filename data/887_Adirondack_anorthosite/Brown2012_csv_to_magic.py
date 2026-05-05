"""
Convert Brown2012_sites.csv to MagIC 3.0 sites.txt and locations.txt, then
combine into a MagIC upload file using ipmag.upload_magic().

Source: Brown, L. L., & McEnroe, S. A. (2012). Paleomagnetism and magnetic
    mineralogy of Grenville metamorphic and igneous rocks, Adirondack
    Highlands, USA. Precambrian Research, 212-213, 57-74.
    https://doi.org/10.1016/j.precamres.2012.04.012

Three poles are produced, matching Brown & McEnroe (2012) Table 3:
    1. Adirondack microcline gneiss pole (14 GMS sites)
    2. Adirondack metamorphic anorthosite pole (11 anorthosite sites + 3
       associated metamorphic rock sites = 14 sites total, per Brown 2012
       Table 3 grouping)
    3. Adirondack fayalite granite pole (8 sites)

Each pole is one row in locations.txt; each site is one row in sites.txt with
its 'location' field pointing to the matching pole.

Age source for the 887 +/- 23 Ma metamorphic anorthosite pole: Zhang, Y.,
Hodgin, E. B., Schmitz, M. D., Schwartz, D., Mohr, M. T., Crowley, J. L., &
Swanson-Hysell, N. L. (in press). Protracted high-grade metamorphism and slow
cooling in the Adirondacks recalibrate Neoproterozoic global paleogeography.
Tectonics. The 887 +/- 23 Ma estimate is for the magnetite-carried remanence
acquisition age in the Marcy massif anorthosite, derived from probabilistic
thermal history inversion of U-Pb zircon, garnet, titanite, rutile, and apatite
data. The 3 'associated metamorphic rocks' sites are folded into the
anorthosite pole following Brown 2012 Table 3, since they share
magnetite-dominated remanence carriers and similar metamorphic history; the
Zhang et al. age constraint applies to magnetite-blocked remanence and is
adopted at the site level for those 3 sites as well. The microcline gneiss
pole (ilmeno-hematite carrier) and fayalite granite pole have no comparable
modern age constraint and are reported with a broad 900 +/- 100 Ma age.

MagIC data model: v3.0 (https://www2.earthref.org/MagIC/data-models/3.0)

Notes on metadata choices:
    - One row per site at dir_tilt_correction = 0 (geographic). These basement
      Adirondack rocks have no bedding tilt to correct. VGPs are populated on
      that single geographic row since each VGP in the source CSV was computed
      from the geographic-coordinate site mean direction.
    - Longitude in the source CSV is signed degrees (negative = degrees West).
      Converted to MagIC 0-360 East convention.
    - The 'dir_comp_name' field in the source CSV ('AR' = Adirondack Reverse,
      'AN' = Adirondack Normal) is not a MagIC controlled-vocabulary field and
      is dropped. Polarity is used internally to antipode normal-polarity VGPs
      to common reverse polarity before Fisher averaging at the pole level.
    - All lab procedures used progressive AF and thermal demagnetization with
      PCA component fitting and Fisher site means (Brown & McEnroe 2012,
      Section 4.1), encoded as LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM.
    - 13 of the 36 sites are from R. Hargraves (per. comm., 2002) reported in
      Brown & McEnroe (2012, Table 1); same DOI applies as the data source.
"""

import csv
from pathlib import Path
import pmagpy.ipmag as ipmag

# --- Configuration -----------------------------------------------------------

CITATION = '10.1016/j.precamres.2012.04.012'

METHOD_CODES = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM'

LOC_GMS = 'Adirondack microcline gneiss'
LOC_ANORTH = 'Adirondack metamorphic anorthosite'
LOC_GRANITE = 'Adirondack fayalite granite'

SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_sigma', 'age_unit',
    'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95',
    'dir_n_samples', 'vgp_lat', 'vgp_lon', 'description'
]

LOC_COLS = [
    'location', 'location_type', 'result_name', 'result_type',
    'result_quality', 'method_codes', 'citations', 'geologic_classes',
    'lithologies', 'lat_s', 'lat_n', 'lon_w', 'lon_e',
    'age', 'age_sigma', 'age_unit',
    'dir_tilt_correction',
    'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k', 'pole_n_sites',
    'sites', 'description'
]

# Map source-CSV lithology strings to per-site MagIC metadata. The 'location'
# field assigns each site to the pole grouping defined in POLES below.
LITHOLOGY_META = {
    'microcline gneiss': {
        'location': LOC_GMS,
        'geologic_classes': 'Metamorphic',
        'geologic_types': 'Not Specified',
        'lithologies': 'Gneiss',
        'age': '900', 'age_sigma': '100',
    },
    'metamorphic anorthosites': {
        'location': LOC_ANORTH,
        'geologic_classes': 'Metamorphic',
        'geologic_types': 'Pluton',
        'lithologies': 'Meta Anorthosite',
        'age': '887', 'age_sigma': '23',
    },
    'associated metamorphic rocks': {
        # Grouped with the anorthosite pole per Brown 2012 Table 3
        'location': LOC_ANORTH,
        'geologic_classes': 'Metamorphic',
        'geologic_types': 'Not Specified',
        'lithologies': 'Gneiss',
        'age': '887', 'age_sigma': '23',
    },
    'post-metamorphic fayalite granites': {
        'location': LOC_GRANITE,
        'geologic_classes': 'Igneous',
        'geologic_types': 'Pluton',
        'lithologies': 'Granite',
        'age': '900', 'age_sigma': '100',
    },
}

# Per-pole metadata used to fill the locations.txt rows. The Fisher-mean pole
# fields (pole_lat/lon/alpha95/k/n_sites) are computed from contributing site
# VGPs at runtime and merged in.
POLES = {
    LOC_GMS: {
        'result_name': 'Adirondack microcline gneiss pole',
        'geologic_classes': 'Metamorphic',
        'lithologies': 'Gneiss',
        'age': '900', 'age_sigma': '100',
        'description': (
            'Mean pole from 14 microcline gneiss sites of the Adirondack '
            'Highlands reported in Brown and McEnroe (2012). Remanence is '
            'carried by ilmeno-hematite. Age 900 +/- 100 Ma reflects the '
            'broad uncertainty in the timing of remanence acquisition during '
            'post-Grenvillian cooling; Brown and McEnroe (2012) estimate ca. '
            '960 Ma based on closure of lamellar magnetism in ilmeno-hematite '
            'using cooling curves of Mezger et al. (1991, 1992). '
            'Normal-polarity site VGPs were antipoded before Fisher averaging.'
        ),
    },
    LOC_ANORTH: {
        'result_name': 'Adirondack metamorphic anorthosite ca. 887 Ma pole',
        'geologic_classes': 'Metamorphic',
        'lithologies': 'Meta Anorthosite:Gneiss',
        'age': '887', 'age_sigma': '23',
        'description': (
            'Mean pole following the Brown and McEnroe (2012) Table 3 '
            'grouping that combines 11 metamorphic anorthosite sites of the '
            'Marcy Massif with 3 associated metamorphic rock sites '
            '(metasyenite, metadiorite, granitic gneiss) for a total of 14 '
            'sites. All 14 sites carry a magnetite-dominated remanence. The '
            '887 +/- 23 Ma age is the magnetite-carried remanence acquisition '
            'age estimated by Zhang et al. (in press, Tectonics) from '
            'probabilistic thermal history inversion of U-Pb zircon, garnet, '
            'titanite, rutile, and apatite data from the Adirondack '
            'Highlands. Normal-polarity site VGPs were antipoded before '
            'Fisher averaging.'
        ),
    },
    LOC_GRANITE: {
        'result_name': 'Adirondack fayalite granite pole',
        'geologic_classes': 'Igneous',
        'lithologies': 'Granite',
        'age': '900', 'age_sigma': '100',
        'description': (
            'Mean pole from 8 post-metamorphic fayalite granite sites of the '
            'Wanakena and Ausable Forks plutons reported in Brown and McEnroe '
            '(2012). Remanence is carried by magnetite. Age 900 +/- 100 Ma '
            'reflects the broad uncertainty in the timing of remanence '
            'acquisition; Brown and McEnroe (2012) estimate ca. 990 Ma using '
            'cooling curves and the magnetite Curie temperature, with the '
            'host plutons crystallizing at ca. 1047 Ma. Normal-polarity site '
            'VGPs were antipoded before Fisher averaging.'
        ),
    },
}


def lon_signed_to_east(lon_signed):
    """Convert signed longitude (negative = West) to MagIC 0-360 East."""
    return float(lon_signed) % 360


def antipode(vgp_lat, vgp_lon):
    """Return the antipodal VGP for inverting normal-polarity directions."""
    return -vgp_lat, (vgp_lon + 180) % 360


def main():
    script_dir = Path(__file__).parent
    csv_path = script_dir / 'Brown2012_sites.csv'

    site_rows = []
    with open(csv_path, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            lith_key = row['lithology'].strip()
            meta = LITHOLOGY_META[lith_key]

            site = row['site'].strip()
            lat = f"{float(row['lat']):.3f}"
            lon = f"{lon_signed_to_east(row['lon']):.3f}"

            out = {
                'site': site,
                'location': meta['location'],
                'result_type': 'i',
                'result_quality': 'g',
                'method_codes': METHOD_CODES,
                'citations': CITATION,
                'geologic_classes': meta['geologic_classes'],
                'geologic_types': meta['geologic_types'],
                'lithologies': meta['lithologies'],
                'lat': lat,
                'lon': lon,
                'age': meta['age'],
                'age_sigma': meta['age_sigma'],
                'age_unit': 'Ma',
                'dir_tilt_correction': '0',
                'dir_dec': row['dir_dec'].strip(),
                'dir_inc': row['dir_inc'].strip(),
                'dir_k': row['k'].strip(),
                'dir_alpha95': row['a95'].strip(),
                'dir_n_samples': row['N'].strip(),
                'vgp_lat': row['vgp_lat'].strip(),
                'vgp_lon': row['vgp_lon'].strip(),
                'description': lith_key,
            }
            # Bookkeeping fields used by write_locations and stripped before write
            out['_polarity'] = row['dir_comp_name'].strip()
            site_rows.append(out)

    # --- Write sites.txt ----------------------------------------------------
    sites_path = script_dir / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in site_rows:
            f.write('\t'.join(r[col] for col in SITE_COLS) + '\n')
    print(f'Wrote {sites_path.name}: {len(site_rows)} site rows')

    # --- Write locations.txt -----------------------------------------------
    write_locations(script_dir, site_rows)


def write_locations(script_dir, site_rows):
    """Compute the three Brown 2012 Table 3 poles and write locations.txt.

    Args:
        script_dir: Path to the output directory.
        site_rows: List of site row dicts with bookkeeping key '_polarity'
            ('AR' or 'AN') and a 'location' field assigning each site to a
            pole grouping.
    """
    loc_rows = []
    for loc_name in [LOC_GMS, LOC_ANORTH, LOC_GRANITE]:
        meta = POLES[loc_name]
        contrib = [r for r in site_rows if r['location'] == loc_name]

        # Bring all VGPs to common (reverse) polarity before Fisher mean
        common_pol_lats, common_pol_lons = [], []
        for r in contrib:
            vlat, vlon = float(r['vgp_lat']), float(r['vgp_lon'])
            if r['_polarity'] == 'AN':
                vlat, vlon = antipode(vlat, vlon)
            common_pol_lats.append(vlat)
            common_pol_lons.append(vlon)

        pole_mean = ipmag.fisher_mean(dec=common_pol_lons,
                                      inc=common_pol_lats)

        lats = [float(r['lat']) for r in contrib]
        lons = [float(r['lon']) for r in contrib]
        site_names = [r['site'] for r in contrib]

        loc_rows.append({
            'location': loc_name,
            'location_type': 'Region',
            'result_name': meta['result_name'],
            'result_type': 'a',
            'result_quality': 'g',
            'method_codes': METHOD_CODES + ':DE-VGP',
            'citations': CITATION,
            'geologic_classes': meta['geologic_classes'],
            'lithologies': meta['lithologies'],
            'lat_s': f'{min(lats):.3f}',
            'lat_n': f'{max(lats):.3f}',
            'lon_w': f'{min(lons):.3f}',
            'lon_e': f'{max(lons):.3f}',
            'age': meta['age'],
            'age_sigma': meta['age_sigma'],
            'age_unit': 'Ma',
            'dir_tilt_correction': '0',
            'pole_lat': f'{pole_mean["inc"]:.1f}',
            'pole_lon': f'{pole_mean["dec"]:.1f}',
            'pole_alpha95': f'{pole_mean["alpha95"]:.1f}',
            'pole_k': f'{pole_mean["k"]:.1f}',
            'pole_n_sites': str(pole_mean['n']),
            'sites': ':'.join(site_names),
            'description': meta['description'],
        })

    loc_path = script_dir / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(LOC_COLS) + '\n')
        for row in loc_rows:
            f.write('\t'.join(row[col] for col in LOC_COLS) + '\n')

    print(f'Wrote {loc_path.name}: {len(loc_rows)} pole rows')
    for row in loc_rows:
        print(f'  {row["location"]}: '
              f'{row["pole_lat"]}°N, {row["pole_lon"]}°E, '
              f'A95={row["pole_alpha95"]}°, N={row["pole_n_sites"]}')


def combine_and_validate(script_dir):
    """Combine MagIC tables into upload.txt and validate.

    Uses ipmag.upload_magic() to find all MagIC tables (sites.txt,
    locations.txt) in the directory, combine them into a single upload.txt,
    and validate against the MagIC data model.
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
