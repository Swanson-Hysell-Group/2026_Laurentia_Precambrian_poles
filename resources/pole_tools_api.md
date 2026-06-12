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

## `_int_or_blank`

```python
_int_or_blank(x)
```

Coerce to an integer; blank for missing values.

---

## `_round_or_blank`

```python
_round_or_blank(x, ndigits=1)
```

Round a numeric value to ``ndigits`` decimals; blank for missing values.

---

## `assess_R2`

```python
assess_R2(sites_tc, pole_mean, verbose=True)
```

Evaluate the Meert et al. (2020) R2 criteria from recreated site data.

Reports the three R2 sub-criteria of Meert et al. (2020) (see
``resources/Meert2020_R_criteria.md``) for the site-level data:

- **(a) Demagnetization** (advisory): at least two stepwise methods
  (AF and thermal), inferred from the MagIC ``method_codes`` (``LP-DIR-AF``
  and ``LP-DIR-T``).
- **(b) Component analysis** (advisory): PCA best-fit lines (``DE-BFL``) or
  great circles (``DE-BFP``), inferred from ``method_codes``.
- **(c) Adequate PSV sampling:** the Deenen et al. (2011) A95 envelope on
  the VGP distribution (``12*N^-0.40 <= A95 <= 82*N^-0.63``) together with
  the statistical thresholds N >= 25 samples, 10 <= K <= 70, B >= 8 sites
  with at least 3 samples per site.

Sub-criteria (a) and (b) are reported for information only and are *not*
used to set the R2 score: MagIC ``method_codes`` are inconsistently applied
across studies, and the two-demag-method requirement is not universally
expected (for example AF is inappropriate for hematite-bearing rocks).
These sub-criteria should be evaluated against the source publication. The
returned R2 score reflects the quantitative PSV-sampling sub-criterion (c),
which is the substantive, reproducible part of R2.

**Parameters**

- **sites_tc** (`pd.DataFrame`) — Tilt-corrected site data with at least ``dir_n_samples`` and (optionally) ``method_codes``.
- **pole_mean** (`dict`) — Fisher VGP mean from ``compute_mean_pole`` with keys ``n`` (number of site VGPs = B), ``k`` (K), and ``alpha95`` (the VGP-distribution A95).
- **verbose** (`bool`) — If True, print a per-item report.

**Returns**

- ``{'a': bool|None, 'b': bool|None, 'c': bool, 'R2': int, 'B': int, 'N_samples': int, 'K': float, 'A95': float, 'min_samples_per_site': int, 'deenen_pass': bool}``. ``R2`` equals the score of sub-criterion (c); ``a``/``b`` are advisory (None if ``method_codes`` is absent).

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

## `corrected_pole_note`

```python
corrected_pole_note(method, plat, plon, f=None, f_range=None, zeta95=None, eta95=None, a95=None, monte_carlo='', extra='')
```

Format a note describing a study's inclination-corrected pole for ``Comment``.

The Nordic f-columns carry a standardized blanket flattening factor (0.6 for
sediments, 1 for crystalline rocks), so a study's own inclination-shallowing
result — which may not reduce to a single f and may carry a fully propagated
(Kent-ellipse) uncertainty — is recorded in the ``Comment`` field instead.
Build the string higher in the notebook and include it (with any other
context) in the ``COMMENT`` passed to ``make_nordic_summary``.

**Parameters**

- **method** (`str`) — correction method, e.g. ``'E/I (Tauxe & Kent, 2004)'``.
- **plat, plon** (`float`) — corrected pole latitude and longitude in degrees.
- **f** (`float or None`) — single flattening factor, if applicable.
- **f_range** (`tuple or None`) — ``(low, high)`` 95% flattening range.
- **zeta95, eta95** (`float or None`) — Kent-ellipse semi-axes in degrees.
- **a95** (`float or None`) — circular A95 in degrees (if not an ellipse).
- **monte_carlo** (`str`) — uncertainty-propagation note, e.g. ``'Pierce et al. (2022) Monte Carlo'``.
- **extra** (`str`) — any trailing text.

**Returns**

- a one-sentence note.

---

## `fetch_magic_contribution`

```python
fetch_magic_contribution(contribution_id, dir_path, verbose=True)
```

Fetch a published MagIC contribution by ID, with a local fallback.

Tries to download contribution ``contribution_id`` from earthref.org/MagIC
with ``ipmag.download_magic_from_id`` and unpack it with
``ipmag.download_magic`` into per-table files (``sites.txt``,
``locations.txt``, ...) in ``dir_path``. The combined contribution file is
saved as ``magic_contribution_<id>.txt`` and thereby cached locally, so
after the first successful run the data remain available offline: if a later
run cannot reach MagIC, the cached combined file is unpacked instead, and if
even that is absent any pre-existing local table files are left in place.

**Parameters**

- **contribution_id** (`str or int`) — MagIC contribution ID (e.g. ``20696``).
- **dir_path** (`str`) — Directory to download into and read the cache from.
- **verbose** (`bool`) — If True, print a status message.

**Returns**

- ``'magic'`` if freshly downloaded from MagIC, ``'cache'`` if the locally cached contribution file was unpacked, or ``'local'`` if it fell back to pre-existing local table files.

---

## `fishqq_vgps`

```python
fishqq_vgps(sites_tc, unify_polarity=True)
```

Fisher Q-Q test on site VGPs to assess the shape of the VGP distribution.

Uses ``ipmag.fishqq`` to test whether the set of site virtual geomagnetic
poles (VGPs) is consistent with a Fisher distribution. The test compares
the longitude (declination) component against a uniform distribution and
the latitude (inclination) component against an exponential distribution
(Fisher et al., 1987), returning the Mu (uniform) and Me (exponential)
statistics with their critical values and a pass/fail message. A VGP
distribution that adequately samples paleosecular variation is commonly
elongate and may be reported as non-Fisherian, so this test is interpreted
alongside the SVEI test (``svei_test_vgps``) rather than as a strict
reliability cutoff.

Sites with NaN in either ``vgp_lon`` or ``vgp_lat`` are dropped. The VGPs
are brought to a common polarity with ``pmag.flip(..., combine=True)`` when
``unify_polarity`` is True so that a dual-polarity unit is treated as one
mode.

**Parameters**

- **sites_tc** (`pd.DataFrame`) — Tilt-corrected site data with columns ``vgp_lon`` and ``vgp_lat``.
- **unify_polarity** (`bool`) — If True, bring VGPs to a common polarity before the test.

**Returns**

- dict or tuple[dict, dict]: The dictionary (or, for a two-mode dataset, pair of dictionaries) returned by ``ipmag.fishqq``, with keys including ``N``, ``Mu``, ``Mu_critical``, ``Me``, ``Me_critical``, and ``Test_result``.

---

## `get_Laurentia_poles`

```python
get_Laurentia_poles(file_name='../data/Laurentia_poles.csv', sheet_name='Laurentia')
```

Loads Laurentia poles and rotates them into a common reference frame.

Poles from Scotland, Greenland, and Svalbard terranes are rotated into the
Laurentia reference frame using published Euler poles. Poles from Laurentia
and Trans-Hudson orogen are kept in their original coordinates. Unrecognized
terranes receive NaN for rotated coordinates.

The default source is ``data/Laurentia_poles.csv`` (exported from the
``Laurentia`` sheet of ``Kringdalen_w_Laurentia.xlsx``). A ``.csv`` path is
read with ``read_csv``; any other extension is read with ``read_excel``
using ``sheet_name``.

The compilation mixes two conventions: some poles report ``A95`` and
``nominal age`` directly, while others (Kringdalen-native) report ``DP``/
``DM`` and ``lomagage``/``himagage`` instead. So that the whole path is
available downstream regardless of convention, the ``A95`` and
``nominal age`` columns are filled, **in memory only**, from the midpoints
of ``DP``/``DM`` and ``lomagage``/``himagage`` wherever they are not given
explicitly; the source columns are left unchanged.

**Parameters**

- **file_name** (`str`) — Path to the pole-data file (CSV by default; an Excel workbook is also accepted). Expected columns include PLAT, PLONG, Terrane, ROCKNAME, and either nominal age / A95 or lomagage / himagage / DP / DM.
- **sheet_name** (`str`) — Sheet to read when ``file_name`` is an Excel workbook.

**Returns**

- pd.DataFrame: Pole data with ``A95`` / ``nominal age`` filled from the DP/DM and lomagage/himagage midpoints where absent, plus added PLAT_rotated and PLONG_rotated columns in the Laurentia reference frame.

---

## `get_Laurentia_stricto_poles`

```python
get_Laurentia_stricto_poles(file_name='../data/Laurentia_poles.csv', sheet_name='Laurentia')
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

## `kent_a95_approx`

```python
kent_a95_approx(zeta95, eta95)
```

Approximate a circular A95 from Kent-ellipse 95% semi-axes.

Inclination-shallowing-corrected sedimentary poles are reported with a Kent
confidence ellipse (semi-axes ``zeta95`` and ``eta95``) that propagates the
flattening-factor uncertainty (e.g. via the Monte Carlo approach of Pierce et
al., 2022). The Nordic Workshop compilation format stores a single circular
A95, so this returns the radius of the circle with the same area as the
ellipse — the geometric mean of the two semi-axes, ``sqrt(zeta95 * eta95)`` —
as the A95 to record there. The full Kent ellipse should still be reported in
the notebook.

**Parameters**

- **zeta95** (`float`) — Major-axis semi-angle of the 95% Kent ellipse (degrees).
- **eta95** (`float`) — Minor-axis semi-angle of the 95% Kent ellipse (degrees).

**Returns**

- Equal-area circular A95 approximation in degrees.

---

## `load_magic_sites`

```python
load_magic_sites(sites_path, drop_bad=True)
```

Loads a MagIC sites.txt file and splits by tilt correction.

Reads a tab-delimited MagIC sites table (skipping the header row) and
returns separate DataFrames for geographic (dir_tilt_correction == 0)
and tilt-corrected (dir_tilt_correction == 100) coordinates.

Sites flagged ``result_quality == 'b'`` (e.g. a duplicate re-measurement of
a flow that is retained in the contribution for completeness but excluded
from the pole) are dropped by default so that pole, direction, and count
calculations use only the accepted sites. Pass ``drop_bad=False`` to keep
them (e.g. to inspect the full contribution).

**Parameters**

- **sites_path** (`str`) — Path to a MagIC-format sites.txt file.
- **drop_bad** (`bool`) — If True (default), exclude rows with ``result_quality == 'b'``.

**Returns**

- tuple[pd.DataFrame, pd.DataFrame]: (sites_geo, sites_tc) DataFrames for geographic and tilt-corrected coordinates respectively.

---

## `make_nordic_summary`

```python
make_nordic_summary(terrane, rockname, sites, dir_mean, pole_mean, study_lon, study_lat, lithology='crystalline', f_factor=None, pole_mean_unflattened=None, component_comment='', tests='', gpmdb_number='', magic_id='', percent_reversed='', demag_code='', R1=None, R2=None, R3=None, R4='', R5=None, R6=None, R7=None, Grade=None, nominal_age=None, lomagage=None, himagage=None, REF_method=None, POLE_AUTHORS=None, YEAR=None, JOURNAL=None, VOLUME=None, VPAGES='', TITLE=None, COMMENT='')
```

Build a one-row Nordic-compilation summary for a pole.

Returns a single-row ``pandas.DataFrame`` whose columns are exactly
``NORDIC_COLUMNS`` in order, so it can be pasted straight into the Nordic
Workshop spreadsheet. Pole position is PLAT = pole latitude
(``pole_mean['inc']``), PLONG = pole longitude (``pole_mean['dec']``); these
are VGP-Fisher-mean poles with a circular confidence (a single A95), so the
oval semi-axes **DP and DM are left blank** rather than set equal to A95.

Flattening (``f`` columns): a standardized blanket flattening factor is used —
**f = 0.6 for ``lithology='sedimentary'``, f = 1 for ``lithology='crystalline'``**
(igneous/metamorphic). The corrected inclination ``INCf = unsquish(INC, f)``
is written to both INCf columns, and the flattening-corrected pole
(PLATf/PLONf, A95f) is computed from the corrected mean direction at the
study locality and written to both f-blocks (identical). For crystalline
rocks (f = 1) the f-block equals the (circular) main block, so DPf/DMf are
likewise blank; a sedimentary flattening correction yields a genuine
confidence oval, so DPf/DMf are reported and A95f = sqrt(DPf*DMf). A study's
own inclination-shallowing determination (which may not be a single f and may
have a propagated ellipse, e.g. E/I in Jacobsville) is recorded in
``COMMENT`` — build it with ``corrected_pole_note``.

``R4`` holds the field-test letter code(s) from ``resources/field_test_codes.md``
(e.g. ``'C'`` baked contact, ``'c'`` inverse baked contact, ``'g'``/``'G'``
conglomerate, ``'f'``/``'F'`` fold, ``'U'`` unconformity), not a number. A
populated R4 (a positive field test) contributes 1 to the ``R`` total; an
empty R4 (or ``'0'``) contributes 0. R1–R3 and R5–R7 are 0/1 integers.

The Van der Voo Q-criteria columns (40, 24, 10, 16, Q2–Q7, Q(7)) are left
blank (this project scores the modern Meert et al. (2020) R-criteria instead).
Numbers round to the tenth, except SLAT/SLONG and f (hundredth) and the ages
(nominal age / lomagage / himagage, rounded to the nearest integer). RESULT# is
built from ``gpmdb_number`` and ``magic_id`` as ``"GPMDB:<n> MagIC:<id>"``.
ROCKNAME appears twice (a deliberate duplicate in the Nordic format).

---

## `plot_Deenen_test`

```python
plot_Deenen_test(mean_pole, figsize=(7, 4), ax=None)
```

Plots the Deenen et al. (2011) A95 envelope with the observed pole.

Shades the acceptable A95 range as a function of N between A95_min and
A95_max curves and overlays the observed pole. The N axis extends to 80
by default but expands automatically if ``mean_pole['n']`` exceeds 80.

**Parameters**

- **mean_pole** (`dict`) — Mean pole dictionary (e.g. from ``ipmag.fisher_mean``) with keys ``n`` and ``alpha95``.
- **figsize** (`tuple`) — Figure size when a new figure is created.
- **ax** (`matplotlib.axes.Axes or None`) — Axis to plot on. If None, a new figure and axis are created.

**Returns**

- matplotlib.axes.Axes: The axis containing the plot.

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
plot_pole_overlap(ROCKNAME, Precambrian_poles, Phanerozoic_poles, pole_plat=None, pole_plon=None, pole_A95=None, pole_age=None, show=True)
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
- **show** (`bool`) — If True (default), call ``plt.show()`` so the figure is rendered reliably in notebooks.

**Returns**

- matplotlib.axes.Axes: The map axis.

---

## `plot_site_directions`

```python
plot_site_directions(sites, color='blue', marker='o', markersize=20, title=None, ax=None, show=True)
```

Plots site mean directions on an equal-area net without unifying polarity.

Builds a direction block from the ``dir_dec``/``dir_inc`` columns (dropping
rows with NaN in either) and plots it with ``ipmag.plot_di``. No polarity
unification or flipping is applied, so mixed-polarity data plot in their
measured directions — downward (positive inclination) directions as filled
symbols and upward (negative inclination) directions as open symbols. This is
useful for inspecting the data before ``compute_mean_direction`` or
``compute_mean_pole`` unify polarity to compute the Fisher mean.

**Parameters**

- **sites** (`pd.DataFrame`) — Site data with columns ``dir_dec`` and ``dir_inc``.
- **color** (`str`) — Symbol color passed to ``ipmag.plot_di``.
- **marker** (`str`) — Symbol marker passed to ``ipmag.plot_di``.
- **markersize** (`int`) — Symbol size passed to ``ipmag.plot_di``.
- **title** (`str`) — Title for the plot. If None, no title is set.
- **ax** (`matplotlib.axes.Axes`) — Existing equal-area axis to draw on. If None, a new net is created with ``ipmag.plot_net``.
- **show** (`bool`) — If True, calls ``plt.show()`` after plotting.

**Returns**

- The direction block as a list of ``[dec, inc]`` pairs (the non-unified data that were plotted).

---

## `plot_site_map`

```python
plot_site_map(sites, zoom_start=4, tiles='OpenStreetMap', color='firebrick', radius=5)
```

Builds an interactive folium map of paleomagnetic site locations.

Longitudes in MagIC sites tables are stored in 0–360° convention; this
function shifts them to the −180/180° convention expected by folium.
Duplicate site rows (e.g., geographic and tilt-corrected entries for
the same site) are collapsed by site name; where both coordinate frames
are present the tilt-corrected row (``dir_tilt_correction == 100``) is
kept so the displayed direction is the tilt-corrected one.

When the site table carries them, the site-mean direction in
tilt-corrected coordinates (declination, inclination, and α95) is shown
both in each marker's hover tooltip and in its click popup (the popup
also reports the site coordinates).

**Parameters**

- **sites** (`pd.DataFrame`) — Site data with columns ``site``, ``lat``, and ``lon`` (longitude in 0–360°). If present, ``dir_dec``, ``dir_inc``, and ``dir_alpha95`` are shown in each popup, and ``dir_tilt_correction`` is used to prefer tilt-corrected rows.
- **zoom_start** (`int`) — Initial zoom level for the folium map.
- **tiles** (`str`) — Folium tile layer name (e.g., 'OpenStreetMap', 'CartoDB positron').
- **color** (`str`) — Outline color of the site markers.
- **radius** (`float`) — Marker radius in pixels.

**Returns**

- folium.Map: Interactive map with a CircleMarker per site, labeled with the site name on hover and a popup showing the coordinates and (where available) the tilt-corrected mean direction.

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

---

## `reversal_test`

```python
reversal_test(sites, plot=True, random_seed=None)
```

Runs the McFadden & McElhinny (1990) and bootstrap reversal tests.

Builds a direction block from the ``dir_dec``/``dir_inc`` columns (dropping
rows with NaN in either) and passes the full mixed-polarity set to both
PmagPy reversal tests, which internally separate the normal and reversed
modes. The McFadden & McElhinny (1990) test reports a Watson V common-mean
statistic with an A/B/C classification (or a negative/indeterminate result);
the bootstrap test compares the two modes' Cartesian components and passes
only if their bootstrapped confidence bounds overlap in x, y, and z.

**Parameters**

- **sites** (`pd.DataFrame`) — Site data with columns ``dir_dec`` and ``dir_inc``.
- **plot** (`bool`) — If True, draws the stereonet (MM1990) and the bootstrap cumulative-distribution plots.
- **random_seed** (`int`) — Seed for the Monte Carlo / bootstrap resampling; pass a fixed value for reproducible notebook output.

**Returns**

- The McFadden & McElhinny (1990) result from ``ipmag.reversal_test_MM1990``, ``(classification, angle, critical_angle, label)``.

---

## `save_nordic_summary`

```python
save_nordic_summary(summary, filename, output_dir='../data/nordic_summaries')
```

Save a Nordic summary (single-row DataFrame) to a CSV file.

Writes ``summary`` as a one-row CSV with the exact ``NORDIC_COLUMNS`` header
(including the intentional duplicate labels) into ``output_dir`` (created if
needed).

**Parameters**

- **summary** (`pd.DataFrame`) — One-row summary from ``make_nordic_summary``.
- **filename** (`str`) — Output CSV filename, typically the notebook name (e.g. '780_Gunbarrel' or '780_Gunbarrel.csv').
- **output_dir** (`str`) — Directory in which to write the CSV.

**Returns**

- The full path to the written CSV file.

---

## `svei_test_vgps`

```python
svei_test_vgps(sites_tc, study_lon, study_lat, model='TK03_GAD', kappa=-1, num_sims=1000, plot=True)
```

SVEI test of the VGP scatter shape against a statistical PSV field model.

Each site VGP is converted to a direction at a common locality
(``study_lon``, ``study_lat``) and the resulting directional distribution
is evaluated with the SVEI test implemented in ``pmagpy.svei``. The test
compares the elongation (E = tau2/tau3) and the azimuth of the minor
eigenvector (V2dec) of the distribution against Monte Carlo realizations of
a giant-Gaussian-process paleosecular-variation field model (default
TK03.GAD; Tauxe & Kent, 2004) at the paleolatitude implied by the data. It
reports whether the observed scatter shape is consistent with
adequately-sampled PSV; for an undeformed unit at a given paleolatitude the
field model predicts a N-S elongation whose magnitude depends on latitude.

**Parameters**

- **sites_tc** (`pd.DataFrame`) — Tilt-corrected site data with columns ``vgp_lon`` and ``vgp_lat``.
- **study_lon** (`float`) — Longitude of the common locality used to convert VGPs to directions, in degrees.
- **study_lat** (`float`) — Latitude of the common locality, in degrees.
- **model** (`str`) — Name of the GGP field model passed to ``svei`` (e.g. ``'TK03_GAD'``, ``'BCE19_GAD'``).
- **kappa** (`float`) — Within-site Fisher precision used in the simulations; -1 specifies infinite kappa (no within-site uncertainty).
- **num_sims** (`int`) — Number of Monte Carlo simulations for the E and V2dec confidence bounds.
- **plot** (`bool`) — Whether ``svei`` makes its diagnostic plots.

**Returns**

- The ``svei.svei_test`` result dictionary, with keys including ``lat``, ``E``, ``Esim_min``, ``Esim_max``, ``E_result``, ``V2dec``, ``V2sim_min``, ``V2sim_max``, ``V2_result``, ``A2I``, ``A2D``, and ``H`` (0 if the field-model null cannot be rejected, 1 if rejected).
