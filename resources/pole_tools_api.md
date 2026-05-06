# `pole_tools` API Reference

Utility functions for Laurentia paleomagnetic pole assessment.

Provides routines for loading and rotating poles into the Laurentia reference
frame, computing mean poles from MagIC site data, evaluating reliability
criteria (Deenen et al., 2011; Meert et al., 2020), and plotting poles in
the context of the Laurentia APWP.

## `Deenen_A_95max`

```python
Deenen_A_95max(N)
```

Calculates the maximum A95 threshold from Deenen et al. (2011).

A95 values above this threshold suggest the data are too dispersed
for a reliable pole.

**Parameters**

- **N** (`int`) — Number of sites (or samples) used in the pole calculation.

**Returns**

- A95_max in degrees.

---

## `Deenen_A_95min`

```python
Deenen_A_95min(N)
```

Calculates the minimum A95 threshold from Deenen et al. (2011).

A95 values below this threshold suggest the data may not adequately
sample paleosecular variation (PSV).

**Parameters**

- **N** (`int`) — Number of sites (or samples) used in the pole calculation.

**Returns**

- A95_min in degrees.

---

## `Deenen_test`

```python
Deenen_test(N, A_95)
```

Evaluates whether A95 falls within the Deenen et al. (2011) envelope.

Tests whether the observed A95 is consistent with adequate sampling of
paleosecular variation by checking against N-dependent A95_min and
A95_max thresholds. Prints a pass/fail message.

**Parameters**

- **N** (`int`) — Number of sites used in the pole calculation.
- **A_95** (`float`) — Observed A95 (95% confidence radius) in degrees.

---

## `R2_test`

```python
R2_test(pole_name, pole_df)
```

Evaluates a paleomagnetic pole against the R2 reliability criteria.

Checks four sub-criteria from Meert et al. (2020) R2: sample number
(N >= 25), site number (B >= 8), Fisher precision parameter
(10 <= K <= 70 for adequate PSV sampling), and the Deenen et al. (2011)
A95 envelope. Prints a pass/fail message for each sub-criterion.

**Parameters**

- **pole_name** (`str`) — Name of the rock unit matching a value in the pole_df 'ROCKNAME' column.
- **pole_df** (`pd.DataFrame`) — Poles with columns ROCKNAME, A95, N, B, and KD.

---

## `compute_mean_direction`

```python
compute_mean_direction(sites_tc, unify_polarity=False, flip=False)
```

Computes the Fisher mean direction from site-level declinations and inclinations.

Sites with NaN in either ``dir_dec`` or ``dir_inc`` are dropped before
averaging. The remaining directions are unified to a single polarity with
``pmag.flip(..., combine=True)``; if ``flip`` is True, that unified set is
then flipped 180° via ``ipmag.do_flip`` before computing the Fisher mean.

**Parameters**

- **sites_tc** (`pd.DataFrame`) — Tilt-corrected site data with columns ``dir_dec`` and ``dir_inc``.
- **unify_polarity** (`bool`) — If True, unifies directions to a single polarity
- **flip** (`bool`) — If True, applies a 180° flip to the polarity-unified directions prior to averaging (e.g. to report the mean in the opposite polarity).

**Returns**

- tuple[list, dict]: ``(dir_block_unified, dir_mean)`` where ``dir_block_unified`` is the list of polarity-unified (and optionally flipped) site directions as ``[dec, inc]`` pairs, and ``dir_mean`` is the Fisher mean from ``ipmag.fisher_mean`` with keys ``dec``, ``inc``, ``n``, ``r``, ``k``, ``alpha95``, and ``csd``.

---

## `compute_mean_direction_from_vgps`

```python
compute_mean_direction_from_vgps(sites_tc, study_lon, study_lat, unify_polarity=False, flip=False)
```

Computes the Fisher mean direction from site VGPs converted to 
directions at a common study location.

Each site VGP (``vgp_lon``, ``vgp_lat``) is converted to a direction
(declination, inclination) at the supplied ``study_lon``/``study_lat`` via
``pmag.vgp_di``. This is appropriate when sites span a small region and a
single representative location is used to express the mean as a direction.
Sites with NaN in either VGP column are dropped before conversion. The
resulting directions are unified to a single polarity with
``pmag.flip(..., combine=True)``; if ``flip`` is True, that unified set is
then flipped 180° via ``ipmag.do_flip`` before computing the Fisher mean.

**Parameters**

- **sites_tc** (`pd.DataFrame`) — Tilt-corrected site data with columns ``vgp_lon`` and ``vgp_lat``.
- **study_lon** (`float`) — Longitude in degrees of the common study site at which directions are computed from the VGPs.
- **study_lat** (`float`) — Latitude in degrees of the common study site.
- **unify_polarity** (`bool`) — If True, unifies directions to a single polarity.
- **flip** (`bool`) — If True, applies a 180° flip to the polarity-unified directions prior to averaging.

**Returns**

- tuple[list, dict]: ``(dir_block_unified, dir_mean)`` where ``dir_block_unified`` is the list of polarity-unified (and optionally flipped) directions at the study site as ``[dec, inc]`` pairs, and ``dir_mean`` is the Fisher mean from ``ipmag.fisher_mean`` with keys ``dec``, ``inc``, ``n``, ``r``, ``k``, ``alpha95``, and ``csd``.

---

## `compute_mean_pole`

```python
compute_mean_pole(sites_tc, unify_polarity=False, flip=False)
```

Computes the Fisher mean VGP pole from site-level VGPs.

Sites with NaN in either ``vgp_lon`` or ``vgp_lat`` are dropped before
averaging. The remaining VGPs are unified to a single polarity with
``pmag.flip(..., combine=True)``; if ``flip`` is True, that unified set is
then flipped 180° via ``ipmag.do_flip`` before computing the Fisher mean.

**Parameters**

- **sites_tc** (`pd.DataFrame`) — Tilt-corrected site data with columns ``vgp_lon`` and ``vgp_lat``.
- **unify_polarity** (`bool`) — If True, unifies VGPs to a single polarity
- **flip** (`bool`) — If True, applies a 180° flip to the polarity-unified VGPs prior to averaging (e.g. to report the mean in the opposite polarity).

**Returns**

- tuple[list, dict]: ``(vgp_block_unified, pole_mean)`` where ``vgp_block`` is the list of site VGPs (optionally unified and/or flipped)  as ``[lon, lat]`` pairs, and ``pole_mean`` is the Fisher mean from ``ipmag.fisher_mean`` with keys ``dec``, ``inc``, ``n``, ``r``, ``k``, ``alpha95``, and ``csd``, where ``dec``/``inc`` correspond to the mean pole longitude/latitude.

---

## `get_Laurentia_poles`

```python
get_Laurentia_poles(file_name='../data/Kringdalen_w_Laurentia.xlsx', sheet_name='Laurentia')
```

Loads Laurentia poles and rotates them into a common reference frame.

Poles from Scotland, Greenland, and Svalbard terranes are rotated into the
Laurentia reference frame using published Euler poles. Poles from Laurentia
and Trans-Hudson orogen are kept in their original coordinates. Unrecognized
terranes receive NaN for rotated coordinates.

**Parameters**

- **file_name** (`str`) — Path to the Excel file containing pole data. Expected columns include PLAT, PLONG, Terrane, ROCKNAME, nominal age, and A95.
- **sheet_name** (`str`) — Name of the sheet to read from the Excel file.

**Returns**

- pd.DataFrame: Original pole data with added PLAT_rotated and PLONG_rotated columns containing poles in the Laurentia reference frame.

---

## `get_Laurentia_stricto_poles`

```python
get_Laurentia_stricto_poles(file_name='../data/Kringdalen_w_Laurentia.xlsx', sheet_name='Laurentia')
```

Returns only poles from the Laurentia terrane (sensu stricto).

Filters the full rotated pole dataset to include only entries where
Terrane == 'Laurentia', excluding Scotland, Greenland, Svalbard, and
Trans-Hudson orogen poles.

**Parameters**

- **file_name** (`str`) — Path to the Excel file containing pole data.
- **sheet_name** (`str`) — Name of the sheet to read from the Excel file.

**Returns**

- pd.DataFrame: Subset of poles with Terrane == 'Laurentia', including rotated coordinates from ``get_Laurentia_poles``.

---

## `load_magic_sites`

```python
load_magic_sites(sites_path)
```

Loads a MagIC sites.txt file and splits by tilt correction.

Reads a tab-delimited MagIC sites table (skipping the header row) and
returns separate DataFrames for geographic (dir_tilt_correction == 0)
and tilt-corrected (dir_tilt_correction == 100) coordinates.

**Parameters**

- **sites_path** (`str`) — Path to a MagIC-format sites.txt file.

**Returns**

- tuple[pd.DataFrame, pd.DataFrame]: (sites_geo, sites_tc) DataFrames for geographic and tilt-corrected coordinates respectively.

---

## `make_nordic_summary`

```python
make_nordic_summary(terrane, rockname, sites, dir_mean, pole_mean, study_lon, study_lat, component_comment='', tests='', f_factor=1, pole_mean_unflattened=None, R1=None, R2=None, R3=None, R4=None, R5=None, R6=None, R7=None, Grade=None, nominal_age=None, lomagage=None, himagage=None, REF_method=None, POLE_AUTHORS=None, YEAR=None, JOURNAL=None, VOLUME=None, VPAGES='', TITLE=None, COMMENT='')
```

---

## `plot_apwp_context`

```python
plot_apwp_context(Laurentia_poles, pole_plat, pole_plon, pole_A95, age_min=540, age_max=1780, central_longitude=160, central_latitude=0, projection='mollweide', excluded_terranes=('Laurentia-Scotland', 'Laurentia-Svalbard'), figsize=(12, 12))
```

Plots a pole in the context of the Laurentia Precambrian APWP.

Shows the Laurentia apparent polar wander path color-coded by age with
the target pole highlighted in green. By default, only includes
Laurentia and Greenland (rotated) poles; Svalbard and Scotland poles
are excluded via ``excluded_terranes``. Uses rotated coordinates
throughout.

**Parameters**

- **Laurentia_poles** (`pd.DataFrame`) — Output of ``get_Laurentia_poles`` with columns PLONG_rotated, PLAT_rotated, A95, nominal age, Terrane, and ROCKNAME.
- **pole_plat** (`float`) — Latitude of the pole to highlight in degrees.
- **pole_plon** (`float`) — Longitude of the pole to highlight in degrees.
- **pole_A95** (`float`) — A95 of the pole to highlight in degrees.
- **age_min** (`float`) — Minimum age for filtering in Ma.
- **age_max** (`float`) — Maximum age for filtering in Ma.
- **central_longitude** (`float`) — Center longitude for the projection.
- **central_latitude** (`float`) — Center latitude for the orthographic projection. Ignored when ``projection='mollweide'``.
- **projection** (`str`) — Map projection to use. Either ``'mollweide'`` (default) or ``'orthographic'``.
- **excluded_terranes** (`tuple[str, ...] or None`) — Terrane labels to exclude from the plotted APWP. Defaults to Scotland and Svalbard. Pass ``None`` or an empty tuple to include all rotated terranes.
- **figsize** (`tuple`) — Figure size as (width, height) in inches.

**Returns**

- matplotlib.axes.Axes: The map axis.

---

## `plot_pole_overlap`

```python
plot_pole_overlap(ROCKNAME, Precambrian_poles, Phanerozoic_poles, pole_plat=None, pole_plon=None, pole_A95=None, pole_age=None)
```

Plots all poles younger than the specified pole in both polarities.

Creates a Mollweide projection map showing Precambrian and Phanerozoic
poles that are younger than the pole identified by ROCKNAME. Both normal
and antipodal polarities are plotted. The target pole is highlighted in
green. This is used for the R7 criterion (Meert et al., 2020) to check
whether the pole resembles any younger pole.

Pole coordinates default to the values in the Precambrian_poles DataFrame
but can be overridden with the optional arguments (e.g. when the pole has
been recalculated from MagIC site data).

**Parameters**

- **ROCKNAME** (`str`) — Name of the rock unit to use as the age cutoff. Must match a value in the Precambrian_poles 'ROCKNAME' column.
- **Precambrian_poles** (`pd.DataFrame`) — Precambrian poles with columns ROCKNAME, nominal age, PLONG_rotated, PLAT_rotated, PLONG, PLAT, and A95.
- **Phanerozoic_poles** (`pd.DataFrame`) — Phanerozoic reference poles with columns Lon, Lat, a95, and Age (e.g. Torsvik et al., 2012).
- **pole_plat** (`float or None`) — Override pole latitude in degrees.
- **pole_plon** (`float or None`) — Override pole longitude in degrees.
- **pole_A95** (`float or None`) — Override pole A95 in degrees.
- **pole_age** (`float or None`) — Override pole age in Ma for filtering.

---

## `plot_site_map`

```python
plot_site_map(sites, zoom_start=4, tiles='OpenStreetMap', color='firebrick', radius=5)
```

Builds an interactive folium map of paleomagnetic site locations.

Longitudes in MagIC sites tables are stored in 0–360° convention; this
function shifts them to the −180/180° convention expected by folium.
Duplicate site rows (e.g., geographic and tilt-corrected entries for
the same site) are collapsed by site name.

**Parameters**

- **sites** (`pd.DataFrame`) — Site data with columns ``site``, ``lat``, and ``lon`` (longitude in 0–360°).
- **zoom_start** (`int`) — Initial zoom level for the folium map.
- **tiles** (`str`) — Folium tile layer name (e.g., 'OpenStreetMap', 'CartoDB positron').
- **color** (`str`) — Outline color of the site markers.
- **radius** (`float`) — Marker radius in pixels.

**Returns**

- folium.Map: Interactive map with a CircleMarker per site, labeled with the site name on hover and a popup showing coordinates.

---

## `plot_vgps_and_pole`

```python
plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=150, central_latitude=0, figsize=(8, 8))
```

Plots individual site VGPs and the mean pole on an orthographic map.

Each VGP is labeled with its site name. The mean pole is shown in red
with its A95 confidence circle.

**Parameters**

- **vgp_block** (`list`) — List of VGPs as [lon, lat] pairs.
- **pole_mean** (`dict`) — Mean pole dictionary from ``ipmag.fisher_mean`` with keys dec, inc, n, alpha95.
- **central_longitude** (`float`) — Center longitude for the orthographic projection.
- **central_latitude** (`float`) — Center latitude for the orthographic projection.
- **figsize** (`tuple`) — Figure size as (width, height) in inches.

**Returns**

- matplotlib.axes.Axes: The orthographic map axis.
