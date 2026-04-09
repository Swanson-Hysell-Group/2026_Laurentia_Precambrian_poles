# Getting Site-Level Paleomagnetic Data into MagIC Format

This guide walks through how to take published site-level paleomagnetic data from a paper and format it for the [MagIC database](https://www2.earthref.org/MagIC) using the [MagIC data model v3.0](https://www2.earthref.org/MagIC/data-models/3.0). The examples use data from the ca. 780 Ma Gunbarrel LIP compiled by [https://doi.org/10.1029/2025JB031762](https://doi.org/10.1029/2025JB031762).

:::{note}
This guide focuses on the two MagIC tables most relevant to paleomagnetic pole compilation: **sites** (individual site directions and VGPs) and **locations** (the mean pole). It does not cover the specimens, samples, or measurements tables, which would be populated when contributing full demagnetization data. While the ideal MagIC contribution includes measurement-level data, having site-level data in a standardized format is valuable in its own right. Site-level directions and VGPs enable advances such as building apparent polar wander paths directly from site-level data with propagation of uncertainty (e.g., [https://doi.org/10.1029/2023GL103436](https://doi.org/10.1029/2023GL103436)).
:::

## The MagIC Data Model Hierarchy

MagIC organizes paleomagnetic data in a hierarchical structure:

```
contribution
  └── locations    ← collection of sites; where the mean pole lives
        └── sites        ← units with a common age and magnetization
              └── samples      ← oriented samples from a site
                    └── specimens    ← sub-samples being measured
                          └── measurements  ← individual measurement steps
```

Two levels of this hierarchy are central to pole compilation:

- A **site** in MagIC is a unit with a common age and magnetization — typically a single lava flow, a dike, a sill, or a sedimentary horizon. The site-level result is a mean direction (and its associated VGP) averaged from the samples collected at that site. Each site is an independent spot reading of the ancient field.

- A **location** is a collection of sites grouped together to produce an averaged result. For paleomagnetic pole compilation, the location-level result is the mean pole — a Fisher average of the site-level VGPs. The location row records this pole along with its statistical parameters (A95, *κ*, *N*) and the list of contributing sites.

When compiling published site-level data for pole calculation, you typically need only the **sites** and **locations** tables. Each table is stored as a tab-delimited `.txt` file with a one-line header identifying the table type.

## What You Need from the Paper

Before starting, extract the following from the publication's data table:

- **Site name** — unique identifier for each site
- **Site coordinates** — latitude (°N) and longitude (°E, 0–360)
- **Directions** — declination and inclination, ideally in both geographic and tilt-corrected coordinates
- **Statistics** — Fisher's *k* (or *κ*) and *α*₉₅, number of samples (*n*)
- **VGP coordinates** — pole latitude (°N) and pole longitude (°E) computed from the tilt-corrected direction
- **Age information** — either a radiometric age with uncertainty or bounding age constraints
- **Lithology and geologic type** — values must be chosen that match the [EarthRef controlled vocabularies](https://www2.earthref.org/vocabularies)
- **Method codes** — what demagnetization and analysis methods were used with them specified using values from the [MagIC method codes list](https://www2.earthref.org/MagIC/method-codes)
- **Citation DOIs** — for every data source

## The `sites.txt` Table

See the Gunbarrel example: {download}`sites.txt <../data/780_Gunbarrel/sites.txt>`

### File Structure

A MagIC `sites.txt` file starts with a header row identifying the table type, followed by the column names, then data rows:

```
tab delimited	sites
site	location	result_type	result_quality	method_codes	...
BT50	Gunbarrel LIP	i	g	LP-DIR-T:DE-BFL:DE-FM	...
```

### Key Columns

The columns used in the Gunbarrel example, with their MagIC data model definitions. Columns marked **required** must be populated for the file to pass validation; **conditional** columns are required depending on which other columns are filled.

| Column | Required | Description | Example Value |
|--------|----------|-------------|---------------|
| `site` | required | Unique site name | `BT50` |
| `location` | required | Location name (must match `locations.txt`) | `Gunbarrel LIP` |
| `result_type` | | `i` = individual site result | `i` |
| `result_quality` | recommended | `g` = good, `b` = bad | `g` |
| `method_codes` | required | Colon-delimited list of method codes | `LP-DIR-T:DE-BFL:DE-FM` |
| `citations` | required | Colon-delimited DOIs | `10.1029/2025JB031762` |
| `geologic_classes` | required | Controlled vocabulary | `Igneous` |
| `geologic_types` | required | Controlled vocabulary | `Volcanic Dike` |
| `lithologies` | required | Controlled vocabulary | `Diabase` |
| `lat` | required | Site latitude (°N, -90 to 90) | `45.028` |
| `lon` | required | Site longitude (°E, 0 to 360) | `250.116` |
| `age` | conditional | Radiometric age (required unless `age_low`/`age_high` given) | `779.5` |
| `age_sigma` | | Age uncertainty, 1σ | `1.5` |
| `age_low` | conditional | Lower age bound (required unless `age` given) | `775` |
| `age_high` | conditional | Upper age bound (required unless `age` given) | `780` |
| `age_unit` | required | Unit for age values | `Ma` |
| `dir_tilt_correction` | | `0` = geographic, `100` = tilt-corrected | `100` |
| `dir_dec` | | Declination (°) | `267.4` |
| `dir_inc` | | Inclination (°) | `-15.1` |
| `dir_k` | | Fisher precision parameter *κ* | `38.1` |
| `dir_alpha95` | | 95% confidence cone (°) | `9.9` |
| `dir_n_samples` | | Number of samples in the mean | `8` |
| `vgp_lat` | | VGP latitude (°N) | `-7.3` |
| `vgp_lon` | | VGP longitude (°E) | `156.5` |
| `vgp_dp` | | VGP semi-axis parallel to site meridian (°) | `5.2` |
| `vgp_dm` | | VGP semi-axis perpendicular to site meridian (°) | `10.2` |
| `bed_dip_direction` | | Dip direction of bedding (°, 0–360) | `135` |
| `bed_dip` | | Dip of bedding measured in the dip direction (°) | `20` |
| `description` | | Free-text notes | |

:::{tip}
Including `bed_dip_direction` and `bed_dip` is recommended when bedding orientation data are available. These fields enable tilt corrections to be recalculated and facilitate fold tests and baked contact tests — field tests that constrain the age of magnetization (R4 in the Meert et al., 2020 reliability criteria).
:::

### Tilt Correction and Multiple Rows Per Site

Each row in `sites.txt` contains one result — a single mean declination and inclination in a specified coordinate system, identified by the `dir_tilt_correction` value. When both geographic and tilt-corrected data are provided, a site will have two rows:

1. **Geographic coordinates** (`dir_tilt_correction` = `0`): Declination and inclination before tilt correction. VGP columns are left empty since VGPs are conventionally calculated from tilt-corrected directions.

2. **Tilt-corrected coordinates** (`dir_tilt_correction` = `100`): Declination and inclination after full tilt correction, with VGP latitude, longitude, dp, and dm populated.

This convention allows both coordinate systems to be preserved in a single table.

:::{admonition} When no tilt correction is applied
For intrusions with no bedding tilt correction (e.g., the Slave craton sills Gun, Mar, and Cal_mean in Gunbarrel LIP example dataset), there is no distinction between geographic and tilt-corrected coordinates. In this case, a single row per site with `dir_tilt_correction` = `0` is sufficient. The Gunbarrel example includes duplicate rows for these sites so that the tilt-corrected query returns all sites uniformly, but this is a convenience choice rather than a requirement.
:::

### Common Conventions and Pitfalls

**Longitude convention:** MagIC uses 0–360°E. If the paper reports longitude in degrees West, convert: `lon_E = 360 - lon_W`. For example, the Ding et al. (2025) source table reports 109.884°W, which becomes 250.116°E.

**Citations as DOIs:** Use the DOI (without the `https://doi.org/` prefix), not a formatted reference string. Multiple citations are colon-delimited: `10.3133/pp1580:10.1130/g19944.1`.

**Method codes:** Colon-delimited codes describing the laboratory and analytical methods. The most common for site-level directional data:
- `LP-DIR-T` — stepwise thermal demagnetization
- `LP-DIR-AF` — stepwise alternating field demagnetization
- `DE-BFL` — best-fit line (PCA) component analysis
- `DE-FM` — Fisher mean
- `DE-VGP` — VGP calculation from mean direction (used at the location level)

The full list of method codes is available at [https://www2.earthref.org/MagIC/method-codes](https://www2.earthref.org/MagIC/method-codes).

**Controlled vocabularies:** The `geologic_classes`, `geologic_types`, and `lithologies` columns use controlled vocabularies from the MagIC data model. Use exact matches from the [EarthRef controlled vocabularies](https://www2.earthref.org/vocabularies). Common values include:
- `geologic_classes`: `Igneous`, `Sedimentary`, `Metamorphic`
- `geologic_types`: `Volcanic Dike`, `Sill`, `Lava Flow`, `Sedimentary Layer`
- `lithologies`: `Diabase`, `Basalt`, `Granite`, `Sandstone`, `Limestone`

**Age handling:** If a site has a radiometric age, populate `age` and `age_sigma` and leave `age_low`/`age_high` empty. If only bracketing constraints are available, populate `age_low` and `age_high` instead. The `age_unit` column is always required.

## The `locations.txt` Table

See the Gunbarrel example: {download}`locations.txt <../data/780_Gunbarrel/locations.txt>`

The `locations.txt` table holds the **mean paleomagnetic pole** calculated from the site-level VGPs. This is the result that gets used in apparent polar wander path construction.

### File Structure

```
tab delimited	locations
location	location_type	result_name	result_type	...
Gunbarrel LIP	Region	Gunbarrel LIP ca. 780 Ma pole	a	...
```

### Key Columns

| Column | Required | Description | Example Value |
|--------|----------|-------------|---------------|
| `location` | required | Location name (matches `sites.txt`) | `Gunbarrel LIP` |
| `location_type` | required | Type of location | `Region` |
| `result_name` | | Descriptive name for the pole result | `Gunbarrel LIP ca. 780 Ma pole` |
| `result_type` | | `a` = averaged result | `a` |
| `result_quality` | recommended | `g` = good | `g` |
| `method_codes` | required | All methods used, plus `DE-VGP` | `DE-BFL:DE-FM:LP-DIR-T:DE-VGP` |
| `citations` | required | All DOIs contributing to the pole | |
| `geologic_classes` | required | Shared geologic class | `Igneous` |
| `lithologies` | required | All lithologies represented | `Diabase` |
| `lat_s` / `lat_n` | required | Geographic bounding box, south/north | `43.818` / `65.667` |
| `lon_w` / `lon_e` | required | Geographic bounding box, west/east | `241.537` / `250.612` |
| `age` | conditional | Nominal pole age (required unless `age_low`/`age_high` given) | `780.0` |
| `age_low` / `age_high` | conditional | Age range for the pole (required unless `age` given) | `778.0` / `782.0` |
| `age_unit` | required | Age unit | `Ma` |
| `dir_tilt_correction` | | Tilt correction applied to site data | `100` |
| `pole_lat` | | Mean pole latitude (°N) | `3.2` |
| `pole_lon` | | Mean pole longitude (°E) | `151.5` |
| `pole_alpha95` | | Pole A95 (°) | `8.0` |
| `pole_k` | | Pole Fisher *κ* | `25.8` |
| `pole_n_sites` | | Number of sites in the pole | `14` |
| `sites` | | Colon-delimited list of site names | `BT50:BT51:BT54:...` |
| `description` | | Free-text description | |

The `pole_lat` and `pole_lon` come from a Fisher mean of the individual site VGPs. The `sites` column lists all site names contributing to the mean, colon-delimited. The `DE-VGP` method code is added at the location level to indicate that the pole was calculated as an average of site-level VGPs.

## Two Ways to Create These Files

### Option 1: Manual Entry in a Spreadsheet

For a small number of sites, the most straightforward approach is to type the data directly into a spreadsheet application (Excel, Google Sheets, LibreOffice Calc):

1. Create a new spreadsheet with the column headers from the tables above.
2. Enter site data from the paper, one row per site per tilt correction.
3. Be careful with the conventions (longitude 0–360°E, colon-delimited method codes, DOIs for citations).
4. Export as tab-delimited text.
5. Add the MagIC header line (`tab delimited	sites`) as the first line of the file.

:::{tip}
The MagIC database provides a [template spreadsheet](https://www2.earthref.org/MagIC/data-models/3.0) that includes all columns with built-in validation for controlled vocabularies. Using the template helps avoid typos in column names and controlled vocabulary values.
:::

### Option 2: Programmatic Conversion with Python

When working with data already in a structured format (CSV, Excel), a Python script can automate the conversion and reduce transcription errors. The Gunbarrel example uses this approach in {download}`Ding2025_csv_to_magic.py <../data/780_Gunbarrel/Ding2025_csv_to_magic.py>`.

The script:

1. **Reads the source CSV** containing published site data (declination, inclination, VGP coordinates, ages, etc.).
2. **Applies unit conversions** — e.g., converting longitude from degrees West to MagIC's 0–360°E convention.
3. **Maps metadata** — assigns geologic types, lithologies, and method codes to each site based on a lookup dictionary.
4. **Generates two rows per site** — one for geographic coordinates (`dir_tilt_correction` = 0) and one for tilt-corrected coordinates (`dir_tilt_correction` = 100).
5. **Writes `sites.txt`** with the MagIC header and tab-delimited data.
6. **Computes the mean pole** from the site VGPs using `ipmag.fisher_mean()` and writes `locations.txt`.

Key elements of the programmatic approach:

```python
# Define the MagIC column order
MAGIC_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_sigma', 'age_low', 'age_high', 'age_unit',
    'dir_tilt_correction', 'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95',
    'dir_n_samples', 'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description'
]

# Convert longitude from degrees West to MagIC 0-360 East
lon_east = 360 - lon_west

# Write with MagIC header
with open('sites.txt', 'w') as f:
    f.write('tab delimited\tsites\n')
    f.write('\t'.join(MAGIC_COLS) + '\n')
    for row in rows:
        f.write('\t'.join(row[col] for col in MAGIC_COLS) + '\n')
```

The advantage of the scripted approach is reproducibility — the conversion can be re-run if errors are found or if upstream data are updated — and the script itself serves as documentation of the choices made during data formatting.

## Validation and Upload

A MagIC contribution is a single `.txt` file containing all of the individual table files concatenated together, separated by `>>>>>>>>>>` on its own line (see the Gunbarrel example: {download}`Gunbarrel-LIP_30.Mar.2026.txt <../data/780_Gunbarrel/Gunbarrel-LIP_30.Mar.2026.txt>`):

```
tab delimited	locations
location	location_type	result_name	...
Gunbarrel LIP	Region	Gunbarrel LIP ca. 780 Ma pole	...
>>>>>>>>>>
tab delimited	sites
site	location	result_type	...
BT50	Gunbarrel LIP	i	...
...
```

You can create this file manually by copying the contents of `sites.txt` and `locations.txt` into a single file with the delimiter between them, or programmatically with PmagPy:

```python
import pmagpy.ipmag as ipmag
ipmag.upload_magic(dir_path='path/to/data/', input_dir_path='path/to/data/')
```

The resulting file can be validated locally using PmagPy's validation functions or through the [MagIC upload interface](https://www2.earthref.org/MagIC/upload), which performs server-side validation against the data model. Either way, validation checks that column names match the data model, controlled vocabulary values are recognized, and required fields are populated.

## Summary Checklist

When converting published site data to MagIC format, verify:

- [ ] Site names are unique identifiers
- [ ] Longitudes are in 0–360°E (not degrees West)
- [ ] Two rows per site if both geographic and tilt-corrected data are available
- [ ] VGP coordinates populated only on the tilt-corrected row
- [ ] Ages use either `age`/`age_sigma` or `age_low`/`age_high` (not both)
- [ ] `age_unit` is filled for every row
- [ ] Method codes, geologic types, and lithologies match MagIC controlled vocabularies
- [ ] Citations are DOIs, colon-delimited for multiple sources
- [ ] The `location` value in `sites.txt` matches the `location` value in `locations.txt`
- [ ] The `sites` column in `locations.txt` lists all contributing site names
- [ ] `result_type` is `i` in `sites.txt` and `a` in `locations.txt`

## MagIC resources 

- [MagIC database](https://www2.earthref.org/MagIC) — the Magnetics Information Consortium database for archiving and searching paleomagnetic data
- [MagIC data model v3.0](https://www2.earthref.org/MagIC/data-models/3.0) — column definitions, validation rules, and table templates
- [MagIC method codes](https://www2.earthref.org/MagIC/method-codes) — full list of laboratory and analytical method codes
- [EarthRef controlled vocabularies](https://www2.earthref.org/vocabularies) — accepted values for geologic classes, types, lithologies, and other controlled fields
- [MagIC upload](https://www2.earthref.org/MagIC/upload) — interface for submitting data to the MagIC database
