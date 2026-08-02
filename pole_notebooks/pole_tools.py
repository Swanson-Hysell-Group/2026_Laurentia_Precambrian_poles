"""Utility functions for Laurentia paleomagnetic pole assessment.

Provides routines for loading and rotating poles into the Laurentia reference
frame, computing mean poles from MagIC site data, evaluating reliability
criteria (Deenen et al., 2011; Meert et al., 2020), and plotting poles in
the context of the Laurentia APWP.
"""

import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pmagpy.svei as svei
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import os

Torsvik2012_poles = pd.read_excel('../data/Torsvik2012.xlsx')
Torsvik2012_Laurentia = Torsvik2012_poles[4:187]

# Euler poles [pole_lat, pole_lon, CCW rotation angle] (degrees) used to rotate
# poles from terranes that have separated from Laurentia back into the Laurentia
# reference frame. Greenland: Cenozoic opening of Baffin Bay and the Labrador Sea
# (Roest & Srivastava, 1989); Scotland: opening of the Atlantic (Torsvik et al.,
# 2017); Svalbard: Maloof et al. (2006).
TERRANE_EULER_POLES = {
    'Laurentia-Greenland':      [67.5, -118.5, -13.8],
    'Laurentia-Greenland-Nain': [67.5, -118.5, -13.8],
    'Laurentia-Scotland':       [78.6, 161.9, -31.0],
    'Laurentia-Svalbard':       [-81.0, 125.0, 68.0],
}

def get_Laurentia_poles(file_name='../data/Laurentia_poles.csv', sheet_name='Laurentia',
                        recent_file_name='../data/nordic_summaries/nordic_summaries_combined.csv',
                        recent_age_max=1779):
    """Loads Laurentia poles and rotates them into a common reference frame.

    Poles from Scotland, Greenland, and Svalbard terranes are rotated into the
    Laurentia reference frame using published Euler poles. Poles from Laurentia
    and Trans-Hudson orogen are kept in their original coordinates. Unrecognized
    terranes receive NaN for rotated coordinates.

    By default the path is a **hybrid**: poles with ``nominal age`` at or below
    ``recent_age_max`` (1779 Ma) are taken from ``recent_file_name`` — the
    project's own recompiled Nordic summaries (``nordic_summaries_combined.csv``,
    the per-notebook poles recreated from MagIC site data) — while older poles
    are taken from ``file_name`` (the legacy ``data/Laurentia_poles.csv``
    compilation, exported from the ``Laurentia`` sheet of
    ``Kringdalen_w_Laurentia.xlsx``). The two sources share the core pole columns
    (Terrane, ROCKNAME, PLAT, PLONG, A95, nominal age, ...); any columns unique
    to one source are NaN-filled on rows from the other. Set ``recent_file_name``
    to ``None`` to use ``file_name`` alone for all ages (the legacy behavior). A
    ``.csv`` path is read with ``read_csv``; any other extension is read with
    ``read_excel`` using ``sheet_name``.

    The compilation mixes two conventions: some poles report ``A95`` and
    ``nominal age`` directly, while others (Kringdalen-native) report ``DP``/
    ``DM`` and ``lomagage``/``himagage`` instead. So that the whole path is
    available downstream regardless of convention, the ``A95`` and
    ``nominal age`` columns are filled, **in memory only**, from the midpoints
    of ``DP``/``DM`` and ``lomagage``/``himagage`` wherever they are not given
    explicitly; the source columns are left unchanged. This fill happens before
    the age split so that both sources are partitioned on a populated
    ``nominal age``.

    Args:
        file_name (str): Path to the legacy pole-data file used for poles older
            than ``recent_age_max`` (CSV by default; an Excel workbook is also
            accepted). Expected columns include PLAT, PLONG, Terrane, ROCKNAME,
            and either nominal age / A95 or lomagage / himagage / DP / DM.
        sheet_name (str): Sheet to read when a file is an Excel workbook.
        recent_file_name (str or None): Path to the recompiled Nordic summaries
            used for poles at or below ``recent_age_max``. ``None`` (or a missing
            file) falls back to ``file_name`` for all ages.
        recent_age_max (float): Age cutoff in Ma. Poles with ``nominal age`` <=
            this value come from ``recent_file_name``; older poles come from
            ``file_name``.

    Returns:
        pd.DataFrame: Pole data with ``A95`` / ``nominal age`` filled from the
        DP/DM and lomagage/himagage midpoints where absent, plus added
        PLAT_rotated and PLONG_rotated columns in the Laurentia reference frame.
    """
    def _load(path):
        if str(path).lower().endswith('.csv'):
            df = pd.read_csv(path)
        else:
            df = pd.read_excel(path, sheet_name=sheet_name)
        cols = set(df.columns)
        if 'nominal age' in cols and {'lomagage', 'himagage'} <= cols:
            df['nominal age'] = df['nominal age'].fillna(
                (df['lomagage'] + df['himagage']) / 2)
        if 'A95' in cols and {'DP', 'DM'} <= cols:
            df['A95'] = df['A95'].fillna((df['DP'] + df['DM']) / 2)
        return df

    Laurentia_poles = _load(file_name)

    if recent_file_name is not None and os.path.exists(recent_file_name):
        recent = _load(recent_file_name)
        # Recompiled Nordic poles for ages <= cutoff; legacy compilation for
        # older poles. A unit whose recompiled age is <= the cutoff but whose
        # legacy age is just above it (e.g. the ECMB dykes, 1779 vs 1780 Ma)
        # would otherwise appear twice, so legacy poles sharing a ROCKNAME with a
        # recompiled pole are dropped — the recompiled version supersedes them.
        recent = recent[recent['nominal age'] <= recent_age_max]
        recent_names = set(recent['ROCKNAME'].dropna())
        Laurentia_poles = Laurentia_poles[
            (Laurentia_poles['nominal age'] > recent_age_max) &
            (~Laurentia_poles['ROCKNAME'].isin(recent_names))]
        Laurentia_poles = pd.concat([recent, Laurentia_poles], ignore_index=True)
        Laurentia_poles = Laurentia_poles.sort_values(
            'nominal age').reset_index(drop=True)
    elif recent_file_name is not None:
        warnings.warn(
            f"recent pole file not found at {recent_file_name}; "
            f"using {file_name} for all ages.")

    plat_rot = []
    plon_rot = []
    for _, row in Laurentia_poles.iterrows():
        terrane = row['Terrane']
        if terrane in TERRANE_EULER_POLES:
            plat, plon = pmag.pt_rot(TERRANE_EULER_POLES[terrane],
                                     [row['PLAT']], [row['PLONG']])
            plat_rot.append(plat[0])
            plon_rot.append(plon[0])
        elif terrane in ('Laurentia', 'Laurentia-Trans-Hudson orogen'):
            plat_rot.append(row['PLAT'])
            plon_rot.append(row['PLONG'])
        else:
            plat_rot.append(float('nan'))
            plon_rot.append(float('nan'))

    Laurentia_poles['PLAT_rotated'] = plat_rot
    Laurentia_poles['PLONG_rotated'] = plon_rot

    return Laurentia_poles

def get_Laurentia_stricto_poles(file_name='../data/Laurentia_poles.csv', sheet_name='Laurentia',
                                recent_file_name='../data/nordic_summaries/nordic_summaries_combined.csv',
                                recent_age_max=1779):
    """Returns only poles from the Laurentia terrane (sensu stricto).

    Filters the full rotated pole dataset to include only entries where
    Terrane == 'Laurentia', excluding Scotland, Greenland, Svalbard, and
    Trans-Hudson orogen poles. The underlying path is the same hybrid
    (recompiled Nordic poles <= ``recent_age_max``, legacy for older poles)
    built by ``get_Laurentia_poles``.

    Args:
        file_name (str): Legacy pole-data file for poles older than
            ``recent_age_max``.
        sheet_name (str): Name of the sheet to read from an Excel workbook.
        recent_file_name (str or None): Recompiled Nordic summaries for poles
            at or below ``recent_age_max``; ``None`` uses ``file_name`` alone.
        recent_age_max (float): Age cutoff in Ma (see ``get_Laurentia_poles``).

    Returns:
        pd.DataFrame: Subset of poles with Terrane == 'Laurentia', including
        rotated coordinates from ``get_Laurentia_poles``.
    """
    Laurentia_poles = get_Laurentia_poles(file_name=file_name, sheet_name=sheet_name,
                                          recent_file_name=recent_file_name,
                                          recent_age_max=recent_age_max)
    Laurentia_stricto_poles = Laurentia_poles[(Laurentia_poles['Terrane']=='Laurentia')]
    return Laurentia_stricto_poles

def plot_pole_overlap(ROCKNAME, Precambrian_poles, Phanerozoic_poles,
                      pole_plat=None, pole_plon=None, pole_A95=None,
                      pole_age=None, show=True):
    """Plots all poles younger than the specified pole in both polarities.

    Creates a Mollweide projection map showing Precambrian and Phanerozoic
    poles that are younger than the pole identified by ROCKNAME. Both normal
    and antipodal polarities are plotted. The target pole is highlighted in
    green. This is used for the R7 criterion (Meert et al., 2020) to check
    whether the pole resembles any younger pole.

    Pole coordinates default to the values in the Precambrian_poles DataFrame
    but can be overridden with the optional arguments (e.g. when the pole has
    been recalculated from MagIC site data).

    Args:
        ROCKNAME (str): Name of the rock unit to use as the age cutoff. Must
            match a value in the Precambrian_poles 'ROCKNAME' column.
        Precambrian_poles (pd.DataFrame): Precambrian poles with columns
            ROCKNAME, nominal age, PLONG_rotated, PLAT_rotated, PLONG,
            PLAT, and A95.
        Phanerozoic_poles (pd.DataFrame): Phanerozoic reference poles with
            columns Lon, Lat, a95, and Age (e.g. Torsvik et al., 2012).
        pole_plat (float or None): Override pole latitude in degrees.
        pole_plon (float or None): Override pole longitude in degrees.
        pole_A95 (float or None): Override pole A95 in degrees.
        pole_age (float or None): Override pole age in Ma for filtering.
        show (bool): If True (default), call ``plt.show()`` so the figure is
            rendered reliably in notebooks.

    Returns:
        matplotlib.axes.Axes: The map axis.
    """

    pole_index = Precambrian_poles.loc[
        Precambrian_poles['ROCKNAME'] == ROCKNAME
    ].index
    has_match = len(pole_index) > 0

    this_row = None
    if has_match:
        this_row = Precambrian_poles.loc[pole_index].iloc[0]

    def _normalized_terrane(value):
        if pd.isna(value):
            return None
        return str(value).strip()

    def _is_laurentian_subblock(terrane):
        terrane = _normalized_terrane(terrane)
        if terrane is None:
            return False
        return terrane.startswith('Laurentia-') and terrane not in (
            'Laurentia-Scotland',
            'Laurentia-Svalbard',
            'Laurentia-Greenland',
        )

    def _row_value_with_fallback(primary, fallback=None):
        if this_row is None:
            return None
        value = this_row.get(primary)
        if pd.notna(value):
            return value
        if fallback is not None:
            return this_row.get(fallback)
        return value

    if pole_age is None:
        if has_match:
            pole_age = this_row.get('nominal age')
            if pd.isna(pole_age) and {'lomagage', 'himagage'} <= set(Precambrian_poles.columns):
                lo = this_row.get('lomagage')
                hi = this_row.get('himagage')
                if pd.notna(lo) and pd.notna(hi):
                    pole_age = 0.5 * (float(lo) + float(hi))
        else:
            raise ValueError(
                f"ROCKNAME {ROCKNAME!r} not found in Precambrian_poles; "
                "provide pole_age when plotting a custom pole."
            )
    if pole_plon is None:
        if has_match:
            pole_plon = _row_value_with_fallback('PLONG_rotated', 'PLONG')
        else:
            raise ValueError(
                f"ROCKNAME {ROCKNAME!r} not found in Precambrian_poles; "
                "provide pole_plon when plotting a custom pole."
            )
    if pole_plat is None:
        if has_match:
            pole_plat = _row_value_with_fallback('PLAT_rotated', 'PLAT')
        else:
            raise ValueError(
                f"ROCKNAME {ROCKNAME!r} not found in Precambrian_poles; "
                "provide pole_plat when plotting a custom pole."
            )
    if pole_A95 is None:
        if has_match:
            pole_A95 = this_row.get('A95')
            if pd.isna(pole_A95) and {'DP', 'DM'} <= set(Precambrian_poles.columns):
                dp = this_row.get('DP')
                dm = this_row.get('DM')
                if pd.notna(dp) and pd.notna(dm):
                    pole_A95 = np.sqrt(float(dp) * float(dm))
        else:
            raise ValueError(
                f"ROCKNAME {ROCKNAME!r} not found in Precambrian_poles; "
                "provide pole_A95 when plotting a custom pole."
            )

    pole_plon = pd.to_numeric(pole_plon, errors='coerce')
    pole_plat = pd.to_numeric(pole_plat, errors='coerce')
    pole_A95 = pd.to_numeric(pole_A95, errors='coerce')

    ax = ipmag.make_mollweide_map(add_land=False, central_longitude=140, figsize=(20,20))

    age_min = 0
    age_max = pole_age

    pole_terrane = _normalized_terrane(this_row.get('Terrane')) if this_row is not None else None
    Precambrian_poles_filtered = Precambrian_poles[
        Precambrian_poles['nominal age'] <= age_max
    ].copy()

    # For pre-amalgamation Laurentian blocks, compare only within the same block
    # until Laurentia is treated as amalgamated in the terrane field.
    if _is_laurentian_subblock(pole_terrane):
        filtered_terranes = Precambrian_poles_filtered['Terrane'].apply(_normalized_terrane)
        keep_mask = (
            filtered_terranes.eq(pole_terrane)
            | filtered_terranes.eq('Laurentia')
            | filtered_terranes.eq('Laurentia-Trans-Hudson orogen')
        )
        Precambrian_poles_filtered = Precambrian_poles_filtered[keep_mask].copy()

    ipmag.plot_poles_colorbar(ax, Phanerozoic_poles['Lon'].tolist(), Phanerozoic_poles['Lat'].tolist(), Phanerozoic_poles['a95'].tolist(), 
                              Phanerozoic_poles['Age'].tolist(),age_min,age_max,colormap='coolwarm',colorbar=False)

    Torsvik2012_Lon_reversed = Phanerozoic_poles['Lon']+180
    Torsvik2012_Lat_reversed = -Phanerozoic_poles['Lat']
    ipmag.plot_poles_colorbar(ax, Torsvik2012_Lon_reversed.tolist(), Torsvik2012_Lat_reversed.tolist(), Phanerozoic_poles['a95'].tolist(), 
                              Phanerozoic_poles['Age'].tolist(),age_min,age_max,marker='s',colormap='coolwarm',colorbar=False)

    for n in Phanerozoic_poles.index:
        ax.text(Phanerozoic_poles['Lon'][n], Phanerozoic_poles['Lat'][n],
                str(int(Phanerozoic_poles['Age'][n])),transform=ccrs.PlateCarree(),fontsize=6)
        ax.text(Torsvik2012_Lon_reversed[n], Torsvik2012_Lat_reversed[n],
                str(int(Phanerozoic_poles['Age'][n])),transform=ccrs.PlateCarree(),fontsize=6)

    ipmag.plot_poles_colorbar(ax, Precambrian_poles_filtered['PLONG'].tolist(), Precambrian_poles_filtered['PLAT'].tolist(), Precambrian_poles_filtered['A95'].tolist(), 
                              Precambrian_poles_filtered['nominal age'].tolist(),age_min,age_max,colormap='coolwarm',colorbar=False)

    Precambrian_poles_filtered_Lon_reversed = Precambrian_poles_filtered['PLONG']+180
    Precambrian_poles_filtered_Lat_reversed = -Precambrian_poles_filtered['PLAT']
    ipmag.plot_poles_colorbar(ax, Precambrian_poles_filtered_Lon_reversed.tolist(), Precambrian_poles_filtered_Lat_reversed.tolist(), 
                              Precambrian_poles_filtered['A95'].tolist(), 
                              Precambrian_poles_filtered['nominal age'].tolist(),age_min,age_max,colormap='coolwarm')

    for n in Precambrian_poles_filtered.index:
        age_label = str(int(Precambrian_poles_filtered['nominal age'][n]))
        ax.text(Precambrian_poles_filtered['PLONG'][n],Precambrian_poles_filtered['PLAT'][n],
                age_label,transform=ccrs.PlateCarree(),fontsize=6)
        ax.text(Precambrian_poles_filtered_Lon_reversed[n],Precambrian_poles_filtered_Lat_reversed[n],
                age_label,transform=ccrs.PlateCarree(),fontsize=6)

    # Draw the A95 confidence circle (best-effort; may fail near geographic poles).
    if np.isfinite(pole_A95):
        try:
            ipmag.plot_pole(ax, pole_plon, pole_plat,
                            pole_A95, filled_pole=True, fill_color='green',
                            fill_alpha=0.5, zorder=998)
            ipmag.plot_pole(ax, 180 + pole_plon, -pole_plat,
                            pole_A95, filled_pole=True, fill_color='green',
                            fill_alpha=0.5, zorder=998)
        except Exception:
            pass  # circle rendering failed; the star below still guarantees visibility

    # Always draw a bright star marker on top of everything — zorder=1000 ensures
    # it is above all other map artists regardless of how many poles are plotted.
    for _plon, _plat in [(pole_plon, pole_plat), (180 + pole_plon, -pole_plat)]:
        ax.plot(_plon, _plat, marker='*', color='lime', markersize=22,
                markeredgecolor='darkgreen', markeredgewidth=1.5,
                transform=ccrs.PlateCarree(), zorder=1000)

    if show:
        plt.show()
    return ax

def plot_apwp_context(Laurentia_poles, pole_plat, pole_plon, pole_A95,
                      age_min=540, age_max=1780, central_longitude=160,
                      central_latitude=0, projection='mollweide',
                      excluded_terranes=('Laurentia-Scotland',
                                         'Laurentia-Svalbard'),
                      figsize=(12, 12)):
    """Plots a pole in the context of the Laurentia Precambrian APWP.

    Shows the Laurentia apparent polar wander path color-coded by age with
    the target pole highlighted in green. By default, only includes
    Laurentia and Greenland (rotated) poles; Svalbard and Scotland poles
    are excluded via ``excluded_terranes``. Uses rotated coordinates
    throughout.

    Args:
        Laurentia_poles (pd.DataFrame): Output of ``get_Laurentia_poles``
            with columns PLONG_rotated, PLAT_rotated, A95, nominal age,
            Terrane, and ROCKNAME.
        pole_plat (float): Latitude of the pole to highlight in degrees.
        pole_plon (float): Longitude of the pole to highlight in degrees.
        pole_A95 (float): A95 of the pole to highlight in degrees.
        age_min (float): Minimum age for filtering in Ma.
        age_max (float): Maximum age for filtering in Ma.
        central_longitude (float): Center longitude for the projection.
        central_latitude (float): Center latitude for the orthographic
            projection. Ignored when ``projection='mollweide'``.
        projection (str): Map projection to use. Either ``'mollweide'``
            (default) or ``'orthographic'``.
        excluded_terranes (tuple[str, ...] or None): Terrane labels to
            exclude from the plotted APWP. Defaults to Scotland and
            Svalbard. Pass ``None`` or an empty tuple to include all
            rotated terranes.
        figsize (tuple): Figure size as (width, height) in inches.

    Returns:
        matplotlib.axes.Axes: The map axis.
    """
    if projection == 'mollweide':
        ax = ipmag.make_mollweide_map(central_longitude=central_longitude,
                                       figsize=figsize)
    elif projection == 'orthographic':
        ax = ipmag.make_orthographic_map(central_longitude=central_longitude,
                                          central_latitude=central_latitude,
                                          figsize=figsize)
    else:
        raise ValueError(
            f"projection must be 'mollweide' or 'orthographic', got {projection!r}"
        )

    if excluded_terranes is None:
        excluded_terranes = ()

    path_poles = Laurentia_poles[
        (Laurentia_poles['nominal age'] >= age_min) &
        (Laurentia_poles['nominal age'] <= age_max) &
        (Laurentia_poles['PLAT_rotated'].notna()) &
        (~Laurentia_poles['Terrane'].isin(excluded_terranes))
    ]

    ipmag.plot_poles_colorbar(ax,
                              path_poles['PLONG_rotated'].tolist(),
                              path_poles['PLAT_rotated'].tolist(),
                              path_poles['A95'].tolist(),
                              path_poles['nominal age'].tolist(),
                              age_min, age_max,
                              colormap='viridis')

    for n in path_poles.index:
        ax.text(path_poles['PLONG_rotated'][n] + 2,
                path_poles['PLAT_rotated'][n] + 2,
                str(int(path_poles['nominal age'][n])),
                transform=ccrs.PlateCarree(), fontsize=6, color='gray')

    ipmag.plot_pole(ax, pole_plon, pole_plat, pole_A95,
                    color='green', markersize=60, filled_pole=True,
                    fill_color='green', fill_alpha=0.4)

    ax.set_title(f'Laurentia APWP ({age_min}–{age_max} Ma) with pole at '
                 f'{pole_plat:.1f}°N, {pole_plon:.1f}°E')
    return ax

def Deenen_A_95min(N):
    """Calculates the minimum A95 threshold from Deenen et al. (2011).

    A95 values below this threshold suggest the data may not adequately
    sample paleosecular variation (PSV).

    Args:
        N (int): Number of sites (or samples) used in the pole calculation.

    Returns:
        float: A95_min in degrees.
    """
    A_95=12*N**(-0.4)
    return A_95
    
def Deenen_A_95max(N):
    """Calculates the maximum A95 threshold from Deenen et al. (2011).

    A95 values above this threshold suggest the data are too dispersed
    for a reliable pole.

    Args:
        N (int): Number of sites (or samples) used in the pole calculation.

    Returns:
        float: A95_max in degrees.
    """
    A_95=82*N**(-0.63)
    return A_95

def Deenen_test(N,A_95):
    """Evaluates whether A95 falls within the Deenen et al. (2011) envelope.

    Tests whether the observed A95 is consistent with adequate sampling of
    paleosecular variation by checking against N-dependent A95_min and
    A95_max thresholds. Prints a pass/fail message.

    Args:
        N (int): Number of sites used in the pole calculation.
        A_95 (float): Observed A95 (95% confidence radius) in degrees.
    """
    Deenen_min = Deenen_A_95min(N)
    Deenen_max = Deenen_A_95max(N)
    
    if A_95 < Deenen_min:
        print('A_95 of ' + str(round(A_95,1)) + ' is too small for Deenen criteria of ' +
              str(round(Deenen_min,1)) + ' for this number of sites')
    elif A_95 > Deenen_max:
        print('A_95 of ' + str(round(A_95,1)) + ' is too large for Deenen criteria of ' +
              str(round(Deenen_max,1)) + ' for this number of sites')
    else:
        print('A_95 of ' + str(round(A_95,1)) + ' passes Deenen et al. (2011) criteria of being between ' +
              str(round(Deenen_min,1)) + ' and ' + str(round(Deenen_max,1)) + ' for this number of sites')
        
def R2_test(pole_name,pole_df):
    """Evaluates a paleomagnetic pole against the R2 reliability criteria.

    Checks four sub-criteria from Meert et al. (2020) R2: sample number
    (N >= 25), site number (B >= 8), Fisher precision parameter
    (10 <= K <= 70 for adequate PSV sampling), and the Deenen et al. (2011)
    A95 envelope. Prints a pass/fail message for each sub-criterion.

    Args:
        pole_name (str): Name of the rock unit matching a value in the
            pole_df 'ROCKNAME' column.
        pole_df (pd.DataFrame): Poles with columns ROCKNAME, A95, N, B,
            and KD.
    """
    this_pole = pole_df[pole_df['ROCKNAME'] == pole_name]
    this_pole.reset_index(inplace=True)
    
    A95 = pd.to_numeric(this_pole['A95'][0], errors='coerce')
    N   = pd.to_numeric(this_pole['N'][0],   errors='coerce')
    B   = pd.to_numeric(this_pole['B'][0],   errors='coerce')
    KD  = pd.to_numeric(this_pole['KD'][0],  errors='coerce')
    
    if N >= 25:
        print('N = ' + str(round(N)) + ' (N ≥ 25; sufficient sample number);')
    else:
        print('N = ' + str(round(N)) + ' (N < 25; insufficient sample number);')  
        
    if pd.isna(B):
        print('B = n/a (published grand mean of unit means; site number not reported);')
    elif B >= 8:
        print('B = ' + str(round(B)) + ' (B ≥ 8; sufficient site number);')
    else:
        print('B = ' + str(round(B)) + ' (B < 8; insufficient site number);')

    if KD >= 70:
        print('K = ' + str(round(KD)) + ' (K ≥ 70; concern about underrepresenting PSV);')
    elif KD >= 10:
        print('K = ' + str(round(KD)) + ' (70 ≥ K ≥ 10; meets PSV criteria);')
    else:
        print('K = ' + str(round(KD)) + ' (10 ≥ K; low K, too dispersed);')

    if not pd.isna(B):
        Deenen_test(B,A95)

def assess_R2(sites_tc, pole_mean, verbose=True):
    """Evaluate the Meert et al. (2020) R2 criteria from recreated site data.

    Reports the three R2 sub-criteria of Meert et al. (2020) (see
    ``resources/Meert2020_R_criteria.md``) for the site-level data:

    - **(a) Demagnetization** (advisory): at least two stepwise methods
      (AF and thermal), inferred from the MagIC ``method_codes`` (``LP-DIR-AF``
      and ``LP-DIR-T``).
    - **(b) Component analysis** (advisory): PCA best-fit lines (``DE-BFL``) or
      great circles (``DE-BFP``), inferred from ``method_codes``.
    - **(c) Adequate PSV sampling:** the Deenen et al. (2011) A95 envelope on
      the VGP distribution (``12*N^-0.40 <= A95 <= 82*N^-0.63``) together with
      the statistical thresholds N >= 25 samples, 10 <= K <= 70, and B >= 8
      sites. The traditional >= 3-samples-per-site guideline is reported for
      information only and is NOT applied as a pass/fail gate: Sapienza et al.
      (2023, doi:10.1029/2023JB027211) show that paleopole accuracy is
      dominated by the number of independent sites (B), with diminishing
      returns from additional in-site samples, so a few sites with < 3 samples
      do not compromise a well-sampled multi-site pole.

    Sub-criteria (a) and (b) are reported for information only and are *not*
    used to set the R2 score: MagIC ``method_codes`` are inconsistently applied
    across studies, and the two-demag-method requirement is not universally
    expected (for example AF is inappropriate for hematite-bearing rocks).
    These sub-criteria should be evaluated against the source publication. The
    returned R2 score reflects the quantitative PSV-sampling sub-criterion (c),
    which is the substantive, reproducible part of R2.

    Args:
        sites_tc (pd.DataFrame): Tilt-corrected site data with at least
            ``dir_n_samples`` and (optionally) ``method_codes``.
        pole_mean (dict): Fisher VGP mean from ``compute_mean_pole`` with keys
            ``n`` (number of site VGPs = B), ``k`` (K), and ``alpha95`` (the
            VGP-distribution A95).
        verbose (bool): If True, print a per-item report.

    Returns:
        dict: ``{'a': bool|None, 'b': bool|None, 'c': bool, 'R2': int,
        'B': int, 'N_samples': int, 'K': float, 'A95': float,
        'min_samples_per_site': int, 'deenen_pass': bool}``. ``R2`` equals the
        score of sub-criterion (c); ``a``/``b`` are advisory (None if
        ``method_codes`` is absent).
    """
    B = int(pole_mean['n'])
    K = pole_mean['k']
    A95 = pole_mean['alpha95']
    N_samples = int(sites_tc['dir_n_samples'].dropna().sum())
    min_per_site = int(sites_tc['dir_n_samples'].dropna().min())

    deenen_min = Deenen_A_95min(B)
    deenen_max = Deenen_A_95max(B)
    deenen_pass = deenen_min <= A95 <= deenen_max

    c_n = N_samples >= 25
    c_k = 10 <= K <= 70
    c_b = B >= 8
    # The >= 3-samples-per-site guideline is reported for information only and
    # is NOT a pass/fail gate (Sapienza et al., 2023): paleopole accuracy is
    # dominated by the number of sites (B), not the number of samples per site.
    c = bool(deenen_pass and c_n and c_k and c_b)

    a = b = None
    if 'method_codes' in sites_tc.columns:
        codes = ':'.join(sites_tc['method_codes'].dropna().astype(str)).upper()
        a = ('LP-DIR-AF' in codes) and ('LP-DIR-T' in codes)
        b = ('DE-BFL' in codes) or ('DE-BFP' in codes)

    # R2 score is driven by the quantitative PSV-sampling sub-criterion (c);
    # (a) and (b) are advisory only (see docstring).
    R2 = int(c)

    if verbose:
        def mark(x):
            return '?' if x is None else ('yes' if x else 'no')
        print('R2 assessment (Meert et al., 2020):')
        print('  (a) two demag methods (AF + thermal)  [advisory, from '
              f'method_codes]: {mark(a)}')
        print('  (b) PCA / great-circle component analysis  [advisory, from '
              f'method_codes]: {mark(b)}')
        print('      note: method codes are inconsistently applied; confirm '
              '(a)/(b) against the source publication.')
        print('  (c) adequate PSV sampling [scored]:')
        print(f'        N = {N_samples} samples (>= 25): '
              f'{"PASS" if c_n else "FAIL"}')
        print(f'        K = {K:.1f} (10-70): {"PASS" if c_k else "FAIL"}')
        print(f'        B = {B} sites (>= 8): {"PASS" if c_b else "FAIL"}')
        print(f'        min {min_per_site} samples/site '
              '[advisory only, not a gate; Sapienza et al., 2023]')
        print(f'        A95 = {A95:.1f} in Deenen envelope '
              f'[{deenen_min:.1f}, {deenen_max:.1f}]: '
              f'{"PASS" if deenen_pass else "FAIL"}')
        print(f'  => R2 (from sub-criterion c) = {R2}')

    return {'a': a, 'b': b, 'c': c, 'R2': R2, 'B': B, 'N_samples': N_samples,
            'K': K, 'A95': A95, 'min_samples_per_site': min_per_site,
            'deenen_pass': deenen_pass}


def plot_Deenen_test(mean_pole, figsize=(7, 4), ax=None):
    """Plots the Deenen et al. (2011) A95 envelope with the observed pole.

    Shades the acceptable A95 range as a function of N between A95_min and
    A95_max curves and overlays the observed pole. The N axis extends to 80
    by default but expands automatically if ``mean_pole['n']`` exceeds 80.

    Args:
        mean_pole (dict): Mean pole dictionary (e.g. from
            ``ipmag.fisher_mean``) with keys ``n`` and ``alpha95``.
        figsize (tuple): Figure size when a new figure is created.
        ax (matplotlib.axes.Axes or None): Axis to plot on. If None, a new
            figure and axis are created.

    Returns:
        matplotlib.axes.Axes: The axis containing the plot.
    """
    n_obs = mean_pole['n']
    n_max = max(80, int(np.ceil(n_obs)) + 5)
    N_range = np.arange(5, n_max + 1)
    A95_min_curve = Deenen_A_95min(N_range)
    A95_max_curve = Deenen_A_95max(N_range)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ax.fill_between(N_range, A95_min_curve, A95_max_curve, alpha=0.2,
                    color='green', label='Acceptable range')
    ax.plot(N_range, A95_min_curve, 'g--', linewidth=1, label=r'$A_{95,min}$')
    ax.plot(N_range, A95_max_curve, 'g--', linewidth=1, label=r'$A_{95,max}$')
    ax.plot(n_obs, mean_pole['alpha95'], 'r*', markersize=15, zorder=5,
            label='pole')
    ax.set_xlabel('Number of sites (N)', fontsize=12)
    ax.set_ylabel(r'$A_{95}$ (°)', fontsize=12)
    ax.set_title('Deenen et al. (2011) test', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(5, n_max)
    ax.set_ylim(0, 25)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    return ax

def kent_a95_approx(zeta95, eta95):
    """Approximate a circular A95 from Kent-ellipse 95% semi-axes.

    Inclination-shallowing-corrected sedimentary poles are reported with a Kent
    confidence ellipse (semi-axes ``zeta95`` and ``eta95``) that propagates the
    flattening-factor uncertainty (e.g. via the Monte Carlo approach of Pierce et
    al., 2022). The Nordic Workshop compilation format stores a single circular
    A95, so this returns the radius of the circle with the same area as the
    ellipse — the geometric mean of the two semi-axes, ``sqrt(zeta95 * eta95)`` —
    as the A95 to record there. The full Kent ellipse should still be reported in
    the notebook.

    Args:
        zeta95 (float): Major-axis semi-angle of the 95% Kent ellipse (degrees).
        eta95 (float): Minor-axis semi-angle of the 95% Kent ellipse (degrees).

    Returns:
        float: Equal-area circular A95 approximation in degrees.
    """
    return float(np.sqrt(zeta95 * eta95))

def fishqq_vgps(sites_tc, unify_polarity=True):
    """Fisher Q-Q test on site VGPs to assess the shape of the VGP distribution.

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

    Args:
        sites_tc (pd.DataFrame): Tilt-corrected site data with columns
            ``vgp_lon`` and ``vgp_lat``.
        unify_polarity (bool): If True, bring VGPs to a common polarity before
            the test.

    Returns:
        dict or tuple[dict, dict]: The dictionary (or, for a two-mode dataset,
        pair of dictionaries) returned by ``ipmag.fishqq``, with keys including
        ``N``, ``Mu``, ``Mu_critical``, ``Me``, ``Me_critical``, and
        ``Test_result``.
    """
    vgp_sites = sites_tc.dropna(subset=['vgp_lon', 'vgp_lat'])
    vgp_block = ipmag.make_di_block(vgp_sites['vgp_lon'].tolist(),
                                    vgp_sites['vgp_lat'].tolist())
    if unify_polarity:
        vgp_block = pmag.flip(vgp_block, combine=True)
    try:
        return ipmag.fishqq(di_block=vgp_block, data_type='poles')
    except TypeError:
        # Backward compatible with pmagpy versions that do not accept data_type.
        return ipmag.fishqq(di_block=vgp_block)


def svei_test_vgps(sites_tc, study_lon, study_lat, model='TK03_GAD',
                   kappa=-1, num_sims=1000, plot=True):
    """SVEI test of the VGP scatter shape against a statistical PSV field model.

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

    Args:
        sites_tc (pd.DataFrame): Tilt-corrected site data with columns
            ``vgp_lon`` and ``vgp_lat``.
        study_lon (float): Longitude of the common locality used to convert
            VGPs to directions, in degrees.
        study_lat (float): Latitude of the common locality, in degrees.
        model (str): Name of the GGP field model passed to ``svei`` (e.g.
            ``'TK03_GAD'``, ``'BCE19_GAD'``).
        kappa (float): Within-site Fisher precision used in the simulations;
            -1 specifies infinite kappa (no within-site uncertainty).
        num_sims (int): Number of Monte Carlo simulations for the E and V2dec
            confidence bounds.
        plot (bool): Whether ``svei`` makes its diagnostic plots.

    Returns:
        dict: The ``svei.svei_test`` result dictionary, with keys including
        ``lat``, ``E``, ``Esim_min``, ``Esim_max``, ``E_result``, ``V2dec``,
        ``V2sim_min``, ``V2sim_max``, ``V2_result``, ``A2I``, ``A2D``, and
        ``H`` (0 if the field-model null cannot be rejected, 1 if rejected).
    """
    dir_block, _ = compute_mean_direction_from_vgps(
        sites_tc, study_lon, study_lat, unify_polarity=False)
    return svei.svei_test(np.array(dir_block), model_name=model, kappa=kappa,
                          num_sims=num_sims, plot=plot)


def fetch_magic_contribution(contribution_id, dir_path, verbose=True):
    """Fetch a published MagIC contribution by ID, with a local fallback.

    Tries to download contribution ``contribution_id`` from earthref.org/MagIC
    with ``ipmag.download_magic_from_id`` and unpack it with
    ``ipmag.download_magic`` into per-table files (``sites.txt``,
    ``locations.txt``, ...) in ``dir_path``. The combined contribution file is
    saved as ``magic_contribution_<id>.txt`` and thereby cached locally, so
    after the first successful run the data remain available offline: if a later
    run cannot reach MagIC, the cached combined file is unpacked instead, and if
    even that is absent any pre-existing local table files are left in place.

    Args:
        contribution_id (str or int): MagIC contribution ID (e.g. ``20696``).
        dir_path (str): Directory to download into and read the cache from.
        verbose (bool): If True, print a status message.

    Returns:
        str: ``'magic'`` if freshly downloaded from MagIC, ``'cache'`` if the
        locally cached contribution file was unpacked, or ``'local'`` if it fell
        back to pre-existing local table files.
    """
    combined = f'magic_contribution_{contribution_id}.txt'
    combined_path = os.path.join(dir_path, combined)

    def _unpack():
        ipmag.download_magic(infile=combined, dir_path=dir_path,
                             input_dir_path=dir_path, print_progress=False)

    try:
        ok, msg = ipmag.download_magic_from_id(contribution_id, directory=dir_path)
        if not ok:
            raise RuntimeError(msg)
        _unpack()
        if verbose:
            print(f'-I- Using MagIC contribution {contribution_id} '
                  f'downloaded from earthref.org')
        return 'magic'
    except Exception as ex:
        if os.path.exists(combined_path):
            if verbose:
                print(f'-W- Could not fetch from MagIC ({ex}); '
                      f'unpacking local cache {combined}')
            _unpack()
            return 'cache'
        if verbose:
            print(f'-W- Could not fetch from MagIC ({ex}) and no local cache; '
                  f'using existing local table files in {dir_path}')
        return 'local'


def load_magic_sites(sites_path, drop_bad=True):
    """Loads a MagIC sites.txt file and splits by tilt correction.

    Reads a tab-delimited MagIC sites table (skipping the header row) and
    returns separate DataFrames for geographic (dir_tilt_correction == 0)
    and tilt-corrected (dir_tilt_correction == 100) coordinates.

    Sites flagged ``result_quality == 'b'`` (e.g. a duplicate re-measurement of
    a flow that is retained in the contribution for completeness but excluded
    from the pole) are dropped by default so that pole, direction, and count
    calculations use only the accepted sites. Pass ``drop_bad=False`` to keep
    them (e.g. to inspect the full contribution).

    Args:
        sites_path (str): Path to a MagIC-format sites.txt file.
        drop_bad (bool): If True (default), exclude rows with
            ``result_quality == 'b'``.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: (sites_geo, sites_tc) DataFrames
        for geographic and tilt-corrected coordinates respectively.
    """
    sites = pd.read_csv(sites_path, sep='\t', skiprows=1)
    if drop_bad and 'result_quality' in sites.columns:
        sites = sites[sites['result_quality'] != 'b'].reset_index(drop=True)
    sites_geo = sites[sites['dir_tilt_correction'] == 0].reset_index(drop=True)
    sites_tc = sites[sites['dir_tilt_correction'] == 100].reset_index(drop=True)
    return sites_geo, sites_tc

def compute_mean_pole(sites_tc, unify_polarity=False, flip=False):
    """Computes the Fisher mean VGP pole from site-level VGPs.

    Sites with NaN in either ``vgp_lon`` or ``vgp_lat`` are dropped before
    averaging. The remaining VGPs are unified to a single polarity with
    ``pmag.flip(..., combine=True)``; if ``flip`` is True, that unified set is
    then flipped 180° via ``ipmag.do_flip`` before computing the Fisher mean.

    Args:
        sites_tc (pd.DataFrame): Tilt-corrected site data with columns
            ``vgp_lon`` and ``vgp_lat``.
        unify_polarity (bool): If True, unifies VGPs to a single polarity
        flip (bool): If True, applies a 180° flip to the polarity-unified VGPs
            prior to averaging (e.g. to report the mean in the opposite
            polarity).

    Returns:
        tuple[list, dict]: ``(vgp_block_unified, pole_mean)`` where
        ``vgp_block`` is the list of site VGPs (optionally unified and/or flipped) 
        as ``[lon, lat]`` pairs, and ``pole_mean`` is the
        Fisher mean from ``ipmag.fisher_mean`` with keys ``dec``, ``inc``,
        ``n``, ``r``, ``k``, ``alpha95``, and ``csd``, where ``dec``/``inc``
        correspond to the mean pole longitude/latitude.
    """
    vgp_sites = sites_tc.dropna(subset=['vgp_lon', 'vgp_lat'])
    vgp_lons = vgp_sites['vgp_lon'].tolist()
    vgp_lats = vgp_sites['vgp_lat'].tolist()
    vgp_block = ipmag.make_di_block(vgp_lons, vgp_lats)
    if unify_polarity:
        vgp_block = pmag.flip(vgp_block, combine=True)
    if flip:
        vgp_block = ipmag.do_flip(di_block=vgp_block)
    pole_mean = ipmag.fisher_mean(di_block=vgp_block)
    return vgp_block, pole_mean

def compute_mean_direction(sites_tc, unify_polarity=False, flip=False):
    """Computes the Fisher mean direction from site-level declinations and inclinations.

    Sites with NaN in either ``dir_dec`` or ``dir_inc`` are dropped before
    averaging. The remaining directions are unified to a single polarity with
    ``pmag.flip(..., combine=True)``; if ``flip`` is True, that unified set is
    then flipped 180° via ``ipmag.do_flip`` before computing the Fisher mean.

    Args:
        sites_tc (pd.DataFrame): Tilt-corrected site data with columns
            ``dir_dec`` and ``dir_inc``.
        unify_polarity (bool): If True, unifies directions to a single polarity
        flip (bool): If True, applies a 180° flip to the polarity-unified
            directions prior to averaging (e.g. to report the mean in the
            opposite polarity).

    Returns:
        tuple[list, dict]: ``(dir_block_unified, dir_mean)`` where
        ``dir_block_unified`` is the list of polarity-unified (and optionally
        flipped) site directions as ``[dec, inc]`` pairs, and ``dir_mean`` is
        the Fisher mean from ``ipmag.fisher_mean`` with keys ``dec``, ``inc``,
        ``n``, ``r``, ``k``, ``alpha95``, and ``csd``.
    """
    dir_sites = sites_tc.dropna(subset=['dir_dec', 'dir_inc'])
    dir_decs = dir_sites['dir_dec'].tolist()
    dir_incs = dir_sites['dir_inc'].tolist()
    dir_block = ipmag.make_di_block(dir_decs, dir_incs)
    if unify_polarity:
        dir_block = pmag.flip(dir_block, combine=True)
    if flip:
        dir_block = ipmag.do_flip(di_block=dir_block)
    dir_mean = ipmag.fisher_mean(di_block=dir_block)
    return dir_block, dir_mean

def compute_mean_direction_from_vgps(sites_tc, study_lon, study_lat, 
                                     unify_polarity=False, flip=False):
    """Computes the Fisher mean direction from site VGPs converted to 
    directions at a common study location.

    Each site VGP (``vgp_lon``, ``vgp_lat``) is converted to a direction
    (declination, inclination) at the supplied ``study_lon``/``study_lat`` via
    ``pmag.vgp_di``. This is appropriate when sites span a small region and a
    single representative location is used to express the mean as a direction.
    Sites with NaN in either VGP column are dropped before conversion. The
    resulting directions are unified to a single polarity with
    ``pmag.flip(..., combine=True)``; if ``flip`` is True, that unified set is
    then flipped 180° via ``ipmag.do_flip`` before computing the Fisher mean.

    Args:
        sites_tc (pd.DataFrame): Tilt-corrected site data with columns
            ``vgp_lon`` and ``vgp_lat``.
        study_lon (float): Longitude in degrees of the common study site at
            which directions are computed from the VGPs.
        study_lat (float): Latitude in degrees of the common study site.
        unify_polarity (bool): If True, unifies directions to a single polarity.
        flip (bool): If True, applies a 180° flip to the polarity-unified
            directions prior to averaging.

    Returns:
        tuple[list, dict]: ``(dir_block_unified, dir_mean)`` where
        ``dir_block_unified`` is the list of polarity-unified (and optionally
        flipped) directions at the study site as ``[dec, inc]`` pairs, and
        ``dir_mean`` is the Fisher mean from ``ipmag.fisher_mean`` with keys
        ``dec``, ``inc``, ``n``, ``r``, ``k``, ``alpha95``, and ``csd``.
    """
    vgp_sites = sites_tc.dropna(subset=['vgp_lon', 'vgp_lat'])
    decs = []
    incs = []
    for vgp_lon, vgp_lat in zip(vgp_sites['vgp_lon'], vgp_sites['vgp_lat']):
        dec, inc = pmag.vgp_di(vgp_lat, vgp_lon, study_lat, study_lon)
        decs.append(dec)
        incs.append(inc)
    dir_block = ipmag.make_di_block(decs, incs)
    if unify_polarity:
        dir_block = pmag.flip(dir_block, combine=True)
    if flip:
        dir_block = ipmag.do_flip(di_block=dir_block)
    dir_mean = ipmag.fisher_mean(di_block=dir_block)
    return dir_block, dir_mean

def plot_site_directions(sites, color='blue', marker='o', markersize=20,
                         title=None, ax=None, show=True):
    """Plots site mean directions on an equal-area net without unifying polarity.

    Builds a direction block from the ``dir_dec``/``dir_inc`` columns (dropping
    rows with NaN in either) and plots it with ``ipmag.plot_di``. No polarity
    unification or flipping is applied, so mixed-polarity data plot in their
    measured directions — downward (positive inclination) directions as filled
    symbols and upward (negative inclination) directions as open symbols. This is
    useful for inspecting the data before ``compute_mean_direction`` or
    ``compute_mean_pole`` unify polarity to compute the Fisher mean.

    Args:
        sites (pd.DataFrame): Site data with columns ``dir_dec`` and ``dir_inc``.
        color (str): Symbol color passed to ``ipmag.plot_di``.
        marker (str): Symbol marker passed to ``ipmag.plot_di``.
        markersize (int): Symbol size passed to ``ipmag.plot_di``.
        title (str, optional): Title for the plot. If None, no title is set.
        ax (matplotlib.axes.Axes, optional): Existing equal-area axis to draw on.
            If None, a new net is created with ``ipmag.plot_net``.
        show (bool): If True, calls ``plt.show()`` after plotting.

    Returns:
        list: The direction block as a list of ``[dec, inc]`` pairs (the
        non-unified data that were plotted).
    """
    dir_sites = sites.dropna(subset=['dir_dec', 'dir_inc'])
    dir_block = ipmag.make_di_block(dir_sites['dir_dec'].tolist(),
                                    dir_sites['dir_inc'].tolist())
    if ax is None:
        ipmag.plot_net()
    ipmag.plot_di(di_block=dir_block, color=color, marker=marker,
                  markersize=markersize)
    if title is not None:
        plt.title(title)
    if show:
        plt.show()
    return dir_block

def reversal_test(sites, plot=True, random_seed=None):
    """Runs the McFadden & McElhinny (1990) and bootstrap reversal tests.

    Builds a direction block from the ``dir_dec``/``dir_inc`` columns (dropping
    rows with NaN in either) and passes the full mixed-polarity set to both
    PmagPy reversal tests, which internally separate the normal and reversed
    modes. The McFadden & McElhinny (1990) test reports a Watson V common-mean
    statistic with an A/B/C classification (or a negative/indeterminate result);
    the bootstrap test compares the two modes' Cartesian components and passes
    only if their bootstrapped confidence bounds overlap in x, y, and z.

    Args:
        sites (pd.DataFrame): Site data with columns ``dir_dec`` and ``dir_inc``.
        plot (bool): If True, draws the stereonet (MM1990) and the bootstrap
            cumulative-distribution plots.
        random_seed (int, optional): Seed for the Monte Carlo / bootstrap
            resampling; pass a fixed value for reproducible notebook output.

    Returns:
        tuple: The McFadden & McElhinny (1990) result from
        ``ipmag.reversal_test_MM1990``, ``(classification, angle,
        critical_angle, label)``.
    """
    dir_sites = sites.dropna(subset=['dir_dec', 'dir_inc'])
    decs = dir_sites['dir_dec'].tolist()
    incs = dir_sites['dir_inc'].tolist()
    dir_block = ipmag.make_di_block(decs, incs)
    # split into two modes the same way the tests do: by the principal
    # eigenvector (pmag.flip -> doprinc), not by inclination sign
    pol1, pol2 = pmag.flip(dir_block)
    print(f'{len(dir_block)} directions (polarity 1: {len(pol1)}, polarity 2: {len(pol2)})\n')

    print('McFadden & McElhinny (1990) reversal test')
    print('-' * 42)
    mm = ipmag.reversal_test_MM1990(dec=decs, inc=incs, plot_stereo=plot,
                                    random_seed=random_seed)
    if plot:
        plt.show()

    print('\nBootstrap reversal test')
    print('-' * 42)
    ipmag.reversal_test_bootstrap(dec=decs, inc=incs, plot_stereo=plot,
                                  random_seed=random_seed)
    if plot:
        plt.show()
    return mm

def plot_vgps_and_pole(vgp_block, pole_mean, central_longitude=150,
                       central_latitude=0, figsize=(8, 8)):
    """Plots individual site VGPs and the mean pole on an orthographic map.

    Each VGP is labeled with its site name. The mean pole is shown in red
    with its A95 confidence circle.

    Args:
        vgp_block (list): List of VGPs as [lon, lat] pairs.
        pole_mean (dict): Mean pole dictionary from ``ipmag.fisher_mean``
            with keys dec, inc, n, alpha95.
        central_longitude (float): Center longitude for the orthographic
            projection.
        central_latitude (float): Center latitude for the orthographic
            projection.
        figsize (tuple): Figure size as (width, height) in inches.

    Returns:
        matplotlib.axes.Axes: The orthographic map axis.
    """
    ax = ipmag.make_orthographic_map(central_longitude=central_longitude,
                                     central_latitude=central_latitude,
                                     figsize=figsize, land_edge_color='none')
    ipmag.plot_vgp(ax, di_block=vgp_block, color='blue', markersize=30, alpha=0.5)
    ipmag.plot_pole(ax, pole_mean['dec'], pole_mean['inc'], pole_mean['alpha95'],
                    color='red', markersize=60, filled_pole=True,
                    fill_color='red', fill_alpha=0.3)
    ax.set_title(f'Mean pole: {pole_mean["inc"]:.1f}°N, '
                 f'{pole_mean["dec"]:.1f}°E, A95={pole_mean["alpha95"]:.1f}°, '
                 f'N={pole_mean["n"]}')
    return ax

def plot_site_map(sites, zoom_start=4, tiles='OpenStreetMap',
                  color='firebrick', radius=5):
    """Builds an interactive folium map of paleomagnetic site locations.

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

    Args:
        sites (pd.DataFrame): Site data with columns ``site``, ``lat``,
            and ``lon`` (longitude in 0–360°). If present, ``dir_dec``,
            ``dir_inc``, and ``dir_alpha95`` are shown in each popup, and
            ``dir_tilt_correction`` is used to prefer tilt-corrected rows.
        zoom_start (int): Initial zoom level for the folium map.
        tiles (str): Folium tile layer name (e.g., 'OpenStreetMap',
            'CartoDB positron').
        color (str): Outline color of the site markers.
        radius (float): Marker radius in pixels.

    Returns:
        folium.Map: Interactive map with a CircleMarker per site, labeled
        with the site name on hover and a popup showing the coordinates and
        (where available) the tilt-corrected mean direction.
    """
    import folium

    dir_cols = [c for c in ('dir_dec', 'dir_inc', 'dir_alpha95')
                if c in sites.columns]
    keep = ['site', 'lat', 'lon'] + dir_cols

    if 'dir_tilt_correction' in sites.columns:
        # Sort so tilt-corrected (100) rows come last, then keep one row per
        # site preferring the tilt-corrected entry.
        ordered = sites.sort_values('dir_tilt_correction')
        site_locs = ordered[keep].drop_duplicates(
            subset='site', keep='last').copy()
    else:
        site_locs = sites[keep].drop_duplicates(subset='site').copy()
    site_locs['lon'] = ((site_locs['lon'] + 180) % 360) - 180

    m = folium.Map(
        location=[site_locs['lat'].mean(), site_locs['lon'].mean()],
        zoom_start=zoom_start,
        tiles=tiles,
    )

    def _fmt(value):
        """Format a numeric direction value, or '–' if missing."""
        try:
            x = float(value)
        except (TypeError, ValueError):
            return '–'
        return '–' if np.isnan(x) else f'{x:.1f}'

    for _, row in site_locs.iterrows():
        popup_html = (f"<b>{row['site']}</b><br>"
                      f"{row['lat']:.3f}°N, {row['lon']:.3f}°E")
        tooltip_html = f"<b>{row['site']}</b>"
        if {'dir_dec', 'dir_inc', 'dir_alpha95'} <= set(dir_cols):
            dir_html = (f"Dec {_fmt(row['dir_dec'])}°, "
                        f"Inc {_fmt(row['dir_inc'])}°, "
                        f"α95 {_fmt(row['dir_alpha95'])}°")
            popup_html += f"<br>{dir_html}"
            tooltip_html += f"<br>{dir_html}"
        folium.CircleMarker(
            location=[row['lat'], row['lon']],
            radius=radius,
            color=color,
            fill=True,
            fill_opacity=0.85,
            popup=folium.Popup(popup_html, max_width=220),
            tooltip=tooltip_html,
        ).add_to(m)

    return m


# Exact column names and order of the Nordic Workshop compilation spreadsheet
# (matches the published per-craton CSVs, e.g. the Congo+Kalahari compilation).
# Note the intentional duplicate labels (a second f / INCf / PLATf / ... block
# and a repeated ROCKNAME) — these are preserved so a summary CSV can be pasted
# directly into the Nordic format.
NORDIC_COLUMNS = [
    'Terrane', 'ROCKNAME', 'GPMDB RESULT#', 'MAGIC #', 'COMPONENT', 'TESTS',
    'TILT', 'SLAT', 'SLONG', 'B', 'N', 'DEC', 'INC', 'KD', 'ED95',
    'PLAT', 'PLONG', 'DP', 'DM', 'A95',
    'f (1.0 vs 0.6)', 'INCf', 'PLATf', 'PLONf', 'DPf', 'DMf', 'A95f',
    'f (1.0/0.8/0.6)', 'INCf', 'PLATf', 'PLONf', 'DPf', 'DMf', 'A95f',
    '%REV', 'DEMAGCODE', 'Q1', '24', '10', '16',
    'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q(7)',
    'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'Rsum', 'Grade', 'Grade E+21',
    'nominal age', 'age diff', 'lomagage', 'himagage', 'REF/method', 'ROCKNAME',
    'POLE AUTHORS', 'YEAR', 'JOURNAL', 'VOLUME', 'VPAGES', 'TITLE', 'Comment',
    'Lithology',
]


def _round_or_blank(x, ndigits=1):
    """Round a numeric value to ``ndigits`` decimals; blank for missing values."""
    if x is None or x == '':
        return ''
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return x
    if np.isnan(xf):
        return ''
    return round(xf, ndigits)


def _int_or_blank(x):
    """Coerce to an integer; blank for missing values."""
    if x is None or x == '':
        return ''
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return x
    if np.isnan(xf):
        return ''
    return int(round(xf))


def corrected_pole_note(method, plat, plon, *, f=None, f_range=None,
                        zeta95=None, eta95=None, a95=None, monte_carlo='', extra=''):
    """Format a note describing a study's inclination-corrected pole for ``Comment``.

    The Nordic f-columns carry a standardized blanket flattening factor (0.6 for
    sediments, 1 for crystalline rocks), so a study's own inclination-shallowing
    result — which may not reduce to a single f and may carry a fully propagated
    (Kent-ellipse) uncertainty — is recorded in the ``Comment`` field instead.
    Build the string higher in the notebook and include it (with any other
    context) in the ``COMMENT`` passed to ``make_nordic_summary``.

    Args:
        method (str): correction method, e.g. ``'E/I (Tauxe & Kent, 2004)'``.
        plat, plon (float): corrected pole latitude and longitude in degrees.
        f (float or None): single flattening factor, if applicable.
        f_range (tuple or None): ``(low, high)`` 95% flattening range.
        zeta95, eta95 (float or None): Kent-ellipse semi-axes in degrees.
        a95 (float or None): circular A95 in degrees (if not an ellipse).
        monte_carlo (str): uncertainty-propagation note, e.g.
            ``'Pierce et al. (2022) Monte Carlo'``.
        extra (str): any trailing text.

    Returns:
        str: a one-sentence note.
    """
    hemi = 'N' if plat >= 0 else 'S'
    pos = f"{abs(plat):.1f}°{hemi}, {plon:.1f}°E"
    fbits = []
    if f is not None:
        fbits.append(f"f = {f:.2f}")
    if f_range is not None:
        fbits.append(f"95% range {f_range[0]:.2f}-{f_range[1]:.2f}")
    fstr = (', ' + ', '.join(fbits)) if fbits else ''
    if zeta95 is not None and eta95 is not None:
        unc = f"Kent ellipse zeta95 = {zeta95:.1f}° / eta95 = {eta95:.1f}°"
    elif a95 is not None:
        unc = f"A95 = {a95:.1f}°"
    else:
        unc = ''
    mc = f", uncertainty propagated via {monte_carlo}" if monte_carlo else ''
    note = f"Study inclination-corrected pole - {method}{fstr}: {pos}"
    if unc:
        note += f", {unc}{mc}"
    note += '.'
    if extra:
        note += ' ' + extra
    return note


def _magic_lithology(sites):
    """Build the ``Lithology`` value from the MagIC site table.

    Concatenates the unique ``lithologies`` terms recorded for the pole's sites,
    splitting MagIC's ``':'`` compound delimiter (e.g. ``'Sandstone:Mudstone'`` ->
    ``'Sandstone; Mudstone'``, ``'Basalt'`` -> ``'Basalt'``). Order of first
    appearance is preserved. Returns ``''`` if the ``lithologies`` column is
    absent or empty; pass ``lithology_note`` to ``make_nordic_summary`` to
    override the value directly.

    Args:
        sites (pd.DataFrame): MagIC ``sites`` table for the pole.

    Returns:
        str: Semicolon-joined unique lithology terms.
    """
    if 'lithologies' not in getattr(sites, 'columns', []):
        return ''
    terms = []
    for val in sites['lithologies'].dropna():
        for part in str(val).split(':'):
            part = part.strip()
            if part and part not in terms:
                terms.append(part)
    return '; '.join(terms)


def make_nordic_summary(terrane,
                        rockname,
                        sites,
                        dir_mean,
                        pole_mean,
                        study_lon,
                        study_lat,
                        lithology='crystalline',
                        f_factor=None,
                        f_factor_2=None,
                        pole_mean_unflattened=None,
                        component_comment='',
                        tests='',
                        gpmdb_number='',
                        magic_id='',
                        percent_reversed='',
                        demag_code='',
                        R1=None,
                        R2=None,
                        R3=None,
                        R4='',
                        R5=None,
                        R6=None,
                        R7=None,
                        Grade=None,
                        Grade_E21=None,
                        nominal_age=None,
                        lomagage=None,
                        himagage=None,
                        REF_method=None,
                        POLE_AUTHORS=None,
                        YEAR=None,
                        JOURNAL=None,
                        VOLUME=None,
                        VPAGES='',
                        TITLE=None,
                        COMMENT='',
                        lithology_note=None):
    """Build a one-row Nordic-compilation summary for a pole.

    Returns a single-row ``pandas.DataFrame`` whose columns are exactly
    ``NORDIC_COLUMNS`` in order, so it can be pasted straight into the Nordic
    Workshop spreadsheet. Pole position is PLAT = pole latitude
    (``pole_mean['inc']``), PLONG = pole longitude (``pole_mean['dec']``); these
    are VGP-Fisher-mean poles with a circular confidence (a single A95), so the
    oval semi-axes **DP and DM are left blank** rather than set equal to A95.

    Flattening: there are two f-blocks. The **first** (``f (1.0 vs 0.6)``) uses a
    standardized blanket flattening factor — **f = 0.6 for ``lithology='sedimentary'``,
    f = 1 for ``lithology='crystalline'``** (igneous/metamorphic). The **second**
    (``f (1.0/0.8/0.6)``) is identical to the first unless a distinct ``f_factor_2``
    is supplied, allowing an intermediate factor (e.g. f = 0.8 for a partly
    compacted unit). In each block the corrected inclination ``INCf = unsquish(INC, f)``
    is written, and the flattening-corrected pole (PLATf/PLONf, A95f) is computed
    from the corrected mean direction at the study locality. For f = 1 the f-block
    equals the (circular) main block, so DPf/DMf are blank; a flattening correction
    (f < 1) yields a genuine confidence oval, so DPf/DMf are reported and
    A95f = sqrt(DPf*DMf). A study's own inclination-shallowing determination (which
    may not be a single f and may have a propagated ellipse, e.g. E/I in
    Jacobsville) is recorded in ``COMMENT`` — build it with ``corrected_pole_note``;
    pass it as ``pole_mean_unflattened`` to place it in the first f-block instead.

    ``R4`` holds the field-test letter code(s) from ``resources/field_test_codes.md``
    (e.g. ``'C'`` baked contact, ``'c'`` inverse baked contact, ``'g'``/``'G'``
    conglomerate, ``'f'``/``'F'`` fold, ``'U'`` unconformity), not a number. A
    populated R4 (a positive field test) contributes 1 to the ``R`` total; an
    empty R4 (or ``'0'``) contributes 0. R1–R3 and R5–R7 are 0/1 integers.

    The Van der Voo Q-criteria columns (Q1, 24, 10, 16, Q2–Q7, Q(7)) are left
    blank (this project scores the modern Meert et al. (2020) R-criteria instead;
    their sum is ``Rsum``). ``Grade`` is this project's reliability grade;
    ``Grade_E21`` is the published Evans et al. (2021) grade (a letter, or ``'--'``
    for units not graded there) and is required so it is populated per unit.
    ``GPMDB RESULT#`` holds ``gpmdb_number`` (integer) and ``MAGIC #`` holds
    ``magic_id`` in their own columns. ``age diff`` = himagage − lomagage. The
    ``Lithology`` column is pulled from the MagIC ``sites`` table (``lithologies``)
    via ``_magic_lithology``, or set directly with ``lithology_note``. Numbers round
    to the tenth, except SLAT/SLONG and f (hundredth) and the ages (nominal age /
    lomagage / himagage / age diff, rounded to the nearest integer). ROCKNAME
    appears twice (a deliberate duplicate in the Nordic format).
    """
    if str(lithology).lower() not in ('crystalline', 'sedimentary'):
        raise ValueError("lithology must be 'crystalline' or 'sedimentary'")
    if sites['dir_tilt_correction'].nunique() != 1:
        warnings.warn(
            "Multiple tilt correction values found in sites data; "
            "cannot determine single tilt correction for summary, choosing first value."
        )
    tilt = sites['dir_tilt_correction'].iloc[0]
    site_n = pole_mean['n']
    sample_n = sites['dir_n_samples'].sum()

    # Backward compatibility: existing notebooks may still pass f_factor and
    # pole_mean_unflattened. If f_factor is provided, it overrides lithology.
    if f_factor is None:
        # blanket flattening factor: 0.6 for sediments, 1 for crystalline rocks
        f = 0.6 if str(lithology).lower() == 'sedimentary' else 1.0
    else:
        f = float(f_factor)

    DEC, INC, DAL = dir_mean['dec'], dir_mean['inc'], dir_mean['alpha95']
    PLAT, PLON, A95 = pole_mean['inc'], pole_mean['dec'], pole_mean['alpha95']
    # These are circular VGP-Fisher-mean poles: the confidence is a single A95,
    # so the oval semi-axes DP/DM are left blank rather than set equal to A95.

    r = _round_or_blank

    def _fblock(f_val):
        """Seven Nordic f-columns [f, INCf, PLATf, PLONf, DPf, DMf, A95f] for f_val."""
        inc_corr = ipmag.unsquish([INC], f_val)[0]
        if f_val == 1.0:
            # no shallowing correction: the f-block equals the (circular) main
            # pole, so DPf/DMf are left blank
            return [r(f_val, 2), r(inc_corr), r(PLAT), r(PLON), '', '', r(A95)]
        # a flattening correction (f < 1) yields a genuine confidence oval: the
        # corrected pole and its DPf/DMf come from the corrected mean direction
        # at the locality, so DPf/DMf are reported
        plon_c, plat_c, dpf, dmf = pmag.dia_vgp(DEC, inc_corr, DAL,
                                                study_lat, study_lon)
        return [r(f_val, 2), r(inc_corr), r(plat_c), r(plon_c),
                r(dpf), r(dmf), r((dpf * dmf) ** 0.5)]

    # First flattening block (f (1.0 vs 0.6)): blanket factor, or a study-specific
    # inclination-corrected pole supplied via pole_mean_unflattened.
    if pole_mean_unflattened is not None:
        inc_corr = ipmag.unsquish([INC], f)[0]
        fblock1 = [r(f, 2), r(inc_corr),
                   r(pole_mean_unflattened.get('inc', PLAT)),
                   r(pole_mean_unflattened.get('dec', PLON)),
                   '', '',                       # circular corrected pole; DPf/DMf blank
                   r(pole_mean_unflattened.get('alpha95', A95))]
    else:
        fblock1 = _fblock(f)

    # Second flattening block (f (1.0/0.8/0.6)): identical to the first unless a
    # distinct f_factor_2 is supplied (e.g. f = 0.8 for a partly compacted unit).
    fblock2 = fblock1 if f_factor_2 is None else _fblock(float(f_factor_2))

    required = {'R1': R1, 'R2': R2, 'R3': R3, 'R5': R5, 'R6': R6, 'R7': R7,
                'Grade': Grade, 'Grade_E21': Grade_E21, 'nominal_age': nominal_age,
                'lomagage': lomagage, 'himagage': himagage, 'REF_method': REF_method,
                'POLE_AUTHORS': POLE_AUTHORS, 'YEAR': YEAR, 'JOURNAL': JOURNAL,
                'VOLUME': VOLUME, 'TITLE': TITLE}
    for name, val in required.items():
        if val is None:
            raise ValueError(f"Required parameter '{name}' was not provided.")

    # R4 = field-test letter code(s); a populated R4 contributes 1 to the R total
    r4_str = '' if R4 in (None, '', '0', 0) else str(R4)
    r4_score = 1 if r4_str else 0
    R_total = (int(R1) + int(R2) + int(R3) + r4_score
               + int(R5) + int(R6) + int(R7))

    # GPMDB RESULT# and MAGIC # occupy their own columns in the new format
    gpmdb_col = _int_or_blank(gpmdb_number)
    magic_col = '' if magic_id in ('', None) else magic_id

    # age diff = himagage - lomagage (both required, so both are present)
    age_diff = _int_or_blank(float(himagage) - float(lomagage))

    # Lithology: pull from the MagIC sites table unless given explicitly
    lithology_col = _magic_lithology(sites) if lithology_note is None else lithology_note

    values = [
        terrane, rockname, gpmdb_col, magic_col, component_comment, tests,
        _int_or_blank(tilt), r(study_lat, 2), r(study_lon, 2),
        _int_or_blank(site_n), _int_or_blank(sample_n),
        r(DEC), r(INC), r(dir_mean['k']), r(DAL),
        r(PLAT), r(PLON), '', '', r(A95),                 # DP/DM blank (circular A95)
        *fblock1,                                         # f (1.0 vs 0.6) block
        *fblock2,                                         # f (1.0/0.8/0.6) block
        _int_or_blank(percent_reversed), demag_code,
        '', '', '', '',                                   # Q1, 24, 10, 16 (blank)
        '', '', '', '', '', '', '',                       # Q2..Q(7) (blank)
        _int_or_blank(R1), _int_or_blank(R2), _int_or_blank(R3),
        r4_str, _int_or_blank(R5), _int_or_blank(R6),
        _int_or_blank(R7), R_total, Grade, Grade_E21,
        _int_or_blank(nominal_age), age_diff,
        _int_or_blank(lomagage), _int_or_blank(himagage), REF_method, rockname,
        POLE_AUTHORS, _int_or_blank(YEAR), JOURNAL, VOLUME, VPAGES, TITLE, COMMENT,
        lithology_col,
    ]
    return pd.DataFrame([values], columns=NORDIC_COLUMNS)


def save_nordic_summary(summary, filename, output_dir='../data/nordic_summaries'):
    """Save a Nordic summary (single-row DataFrame) to a CSV file.

    Writes ``summary`` as a one-row CSV with the exact ``NORDIC_COLUMNS`` header
    (including the intentional duplicate labels) into ``output_dir`` (created if
    needed).

    Args:
        summary (pd.DataFrame): One-row summary from ``make_nordic_summary``.
        filename (str): Output CSV filename, typically the notebook name
            (e.g. '780_Gunbarrel' or '780_Gunbarrel.csv').
        output_dir (str): Directory in which to write the CSV.

    Returns:
        str: The full path to the written CSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    if not filename.endswith('.csv'):
        filename += '.csv'
    output_path = os.path.join(output_dir, filename)
    if not isinstance(summary, pd.DataFrame):
        summary = pd.DataFrame([summary])
    # utf-8-sig writes a BOM so Excel reads the UTF-8 special characters
    # (°, ±, en/em-dashes, etc.) correctly instead of as mojibake.
    summary.to_csv(output_path, index=False, encoding='utf-8-sig')
    return output_path


# --- Kent mean poles for sedimentary units ---------------------------------
# The Nordic compilation records the pole as measured (the f = 1, uncorrected
# position) with a circular A95. Detrital remanence in sediments is shallowed
# during compaction, so for a sedimentary unit that position is a minimum
# paleolatitude. The companion Kent table records, for every sedimentary pole,
# an inclination-shallowing-corrected mean pole summarized as a Kent (1982)
# distribution whose 95% confidence ellipse propagates the uncertainty in the
# flattening factor f following Pierce et al. (2022). Two routes populate it:
#
#   quantified   the study determined f from the directional distribution itself
#                (E/I, Tauxe & Kent, 2004; or SVEI, Tauxe et al., 2024) and
#                propagated it with ``ipmag.find_ei_kent`` / ``find_svei_kent``.
#   compilation  no f could be determined from the data, so f is resampled from
#                the compiled distribution of measured flattening factors of
#                Pierce et al. (2022) with ``ipmag.find_compilation_kent``.

KENT_COLUMNS = [
    'Terrane', 'ROCKNAME', 'Lithology', 'SLAT', 'SLONG',
    'nominal age', 'lomagage', 'himagage', 'Grade',
    'f method', 'f source', 'f', 'f low', 'f high',
    'PLAT uncorrected', 'PLONG uncorrected', 'A95 uncorrected',
    'PLAT', 'PLONG', 'Zdec', 'Zinc', 'Zeta', 'Edec', 'Einc', 'Eta', 'A95',
    'n resampled', 'Comment',
]

# Values accepted for ``f_method``; the first two are study determinations of f
# from the data, the third resamples the Pierce et al. (2022) compilation.
KENT_METHODS = {
    'E/I': 'E/I (Tauxe & Kent, 2004)',
    'SVEI': 'SVEI (Tauxe et al., 2024)',
    'compilation': 'compilation f (Pierce et al., 2022)',
}


def compilation_kent_pole(pole_mean, study_lon, study_lat, random_seed=1,
                          map_central_longitude=None, map_central_latitude=None,
                          **kwargs):
    """Kent mean pole from the Pierce et al. (2022) compilation of f factors.

    Thin wrapper around ``ipmag.find_compilation_kent`` for a sedimentary pole
    whose flattening factor could not be determined from its own directional
    distribution (too few directions for E/I or SVEI, or only a published mean
    pole is available). Flattening factors are resampled from the compilation of
    measured f values of Pierce et al. (2022); each resample unsquishes the mean
    direction implied by the pole at the sampling locality, and the resulting
    swarm of mean poles is summarized as a Kent distribution. The uncorrected
    pole and the Kent ellipse are plotted on an orthographic map centered on the
    pole unless a different center is given.

    Args:
        pole_mean (dict): Uncorrected mean pole, as returned by
            ``compute_mean_pole`` — ``dec`` (pole longitude), ``inc`` (pole
            latitude), and ``alpha95``.
        study_lon (float): Sampling-locality longitude in degrees east.
        study_lat (float): Sampling-locality latitude in degrees.
        random_seed (int): Seed for the resampling, so a notebook rerun and the
            committed Kent summary agree. Pass ``None`` for an unseeded run.
        map_central_longitude (float or None): Map center; defaults to the
            uncorrected pole longitude.
        map_central_latitude (float or None): Map center; defaults to the
            uncorrected pole latitude.
        **kwargs: Passed through to ``ipmag.find_compilation_kent``.

    Returns:
        dict: The Kent statistics dictionary (``dec``, ``inc``, ``Zdec``,
        ``Zinc``, ``Zeta``, ``Edec``, ``Einc``, ``Eta``, ``n``).

    Note:
        This needs a PmagPy carrying the fix to ``find_compilation_kent``: it
        previously iterated over ``len(f_from_compilation)`` rather than ``n``,
        so the f distribution was sampled 70 times regardless of ``n`` and the
        Kent mean moved by up to a degree between seeds. With the fix the
        default ``n = 10000`` gives the documented ``n * n_fish`` = 1,000,000
        resampled mean poles (about a minute per pole), and repeat runs agree to
        a few hundredths of a degree. ``random_seed`` is fixed regardless so the
        committed table is exactly reproducible.
    """
    plon, plat = pole_mean['dec'], pole_mean['inc']
    if map_central_longitude is None:
        map_central_longitude = plon
    if map_central_latitude is None:
        map_central_latitude = plat
    kent = ipmag.find_compilation_kent(
        plon, plat, pole_mean['alpha95'], study_lon % 360, study_lat,
        map_central_longitude=map_central_longitude,
        map_central_latitude=map_central_latitude,
        random_seed=random_seed, **kwargs)
    _rasterize_dense_scatter(plt.gca())
    return kent


# A scatter larger than this is rasterized rather than left as vector art.
_RASTERIZE_ABOVE = 10000


def _rasterize_dense_scatter(ax, threshold=_RASTERIZE_ABOVE):
    """Rasterize any very dense scatter on ``ax``, leaving the rest vector.

    The resampled-mean-pole swarm drawn by ``find_compilation_kent`` is
    ``n * n_fish`` points (1,000,000 at the defaults). Left as vector art that
    is tens of megabytes per figure in a notebook saved with the SVG inline
    backend, which the pole notebooks use. The swarm is a density cloud at
    alpha = 0.002, so nothing is lost by rasterizing it; the map, the ellipse,
    and the pole markers stay vector.
    """
    for collection in ax.collections:
        try:
            n_points = len(collection.get_offsets())
        except (AttributeError, TypeError):
            continue
        if n_points > threshold:
            collection.set_rasterized(True)


def plot_kent_ellipse(map_axis, kent, **kwargs):
    """Draw a Kent mean pole and its 95% confidence ellipse on the right hemisphere.

    ``ipmag.plot_pole_ellipse`` passes ``lower`` through to
    ``pmagplotlib.plot_ell``, whose default (``lower=True``) returns the
    lower-hemisphere projection of the ellipse. For a pole at a **southern**
    latitude that projection is the antipode, so the ellipse is drawn about
    170 deg from the pole it belongs to — and silently, because the pole marker
    is plotted separately and still lands in the right place. Selecting
    ``lower`` from the sign of the pole latitude keeps the ellipse around its
    own pole. Use this rather than calling ``plot_pole_ellipse`` directly when a
    pole may be at a southern latitude; ``scripts/build_apwp_figure.py`` carries
    the same wrapper for the compilation figures, and
    ``ipmag.find_compilation_kent`` now applies the same selection to the
    diagnostic map it draws.

    Args:
        map_axis: Cartopy map axis to draw on.
        kent (dict): Kent statistics with ``dec``, ``inc``, ``Zdec``, ``Zinc``,
            ``Zeta``, ``Edec``, ``Einc``, ``Eta``.
        **kwargs: Passed through to ``ipmag.plot_pole_ellipse``.

    Returns:
        The map axis.
    """
    return ipmag.plot_pole_ellipse(map_axis, kent,
                                   lower=kent['inc'] >= 0, **kwargs)


def compilation_f_summary(f_from_compilation=None, percentiles=(2.5, 97.5)):
    """Median and 95% range of the compiled flattening factors.

    Summarizes the distribution of measured f values that
    ``ipmag.find_compilation_kent`` resamples (Pierce et al., 2022, Table S1,
    exposed as ``ipmag.PIERCE2022_F_COMPILATION``). For a pole whose own data
    cannot constrain f, this is the f that is being assumed: the median is the
    central value of the compiled population and the percentiles bound the
    middle 95% of it. It is a property of the compilation, identical for every
    such pole, and is **not** a determination for the unit — unlike the E/I or
    SVEI values, which are.

    Args:
        f_from_compilation (list or None): Flattening factors to summarize;
            defaults to the Pierce et al. (2022) compilation.
        percentiles (tuple): Lower and upper percentiles (default 2.5, 97.5).

    Returns:
        tuple[float, tuple[float, float]]: ``(median, (low, high))``.
    """
    if f_from_compilation is None:
        f_from_compilation = getattr(ipmag, 'PIERCE2022_F_COMPILATION', None)
        if f_from_compilation is None:
            raise AttributeError(
                'this PmagPy does not expose ipmag.PIERCE2022_F_COMPILATION; '
                'pass f_from_compilation explicitly or update PmagPy')
    values = np.asarray(f_from_compilation, dtype=float)
    low, high = np.percentile(values, percentiles)
    return float(np.median(values)), (float(low), float(high))


def make_kent_summary(nordic_summary, kent, f_method, f_source,
                      f=None, f_range=None, comment=''):
    """Build a one-row Kent-pole summary for a sedimentary pole.

    Pairs the inclination-shallowing-corrected Kent mean pole with the
    uncorrected pole recorded in the unit's Nordic summary, so the two are
    always reported together and the size of the correction is explicit. The
    terrane, unit name, lithology, locality, age, and grade are read from
    ``nordic_summary`` rather than restated, so the two tables cannot drift.

    ``A95`` is the equal-area circular approximation to the Kent ellipse
    (``sqrt(Zeta * Eta)``, see ``kent_a95_approx``), carried so that consumers
    that can only take a circular confidence — the paleolatitude figures and the
    pole map — have a defensible radius to use; the full ellipse is in the
    ``Zdec``/``Zinc``/``Zeta`` and ``Edec``/``Einc``/``Eta`` columns.

    Args:
        nordic_summary (pd.DataFrame): The unit's one-row Nordic summary from
            ``make_nordic_summary``.
        kent (dict): Kent statistics for the corrected pole, from
            ``ipmag.find_ei_kent``, ``ipmag.find_svei_kent``, or
            ``compilation_kent_pole``. A published Kent ellipse can be passed as
            a dictionary with the same keys.
        f_method (str): Key of ``KENT_METHODS`` — ``'E/I'``, ``'SVEI'``, or
            ``'compilation'``.
        f_source (str): Where the correction comes from, e.g.
            ``'this study'`` or ``'Zhang et al. (2024)'``.
        f (float or None): Flattening factor, where the study determined one.
            For ``f_method='compilation'`` this defaults to the **median** of
            the resampled compilation (see ``compilation_f_summary``) — the f
            that is in effect being assumed, not one determined for this unit.
        f_range (tuple or None): ``(low, high)`` 95% bounds on f; for
            ``f_method='compilation'`` it defaults to the 2.5th/97.5th
            percentiles of the compilation.
        comment (str): Any further context.

    Returns:
        pd.DataFrame: One-row summary with columns ``KENT_COLUMNS``.
    """
    if f_method not in KENT_METHODS:
        raise ValueError(f"f_method must be one of {sorted(KENT_METHODS)}")
    if f_method == 'compilation' and f is None:
        # report the compiled distribution that was resampled rather than
        # leaving the f columns blank; it is the same for every such pole
        f_default, f_range_default = compilation_f_summary()
        f = f_default
        if f_range is None:
            f_range = f_range_default
    row = nordic_summary.iloc[0]

    def src(name):
        """First column of ``name`` — NORDIC_COLUMNS repeats some labels."""
        return row.iloc[NORDIC_COLUMNS.index(name)]

    r = _round_or_blank
    values = [
        src('Terrane'), src('ROCKNAME'), src('Lithology'), src('SLAT'), src('SLONG'),
        src('nominal age'), src('lomagage'), src('himagage'), src('Grade'),
        KENT_METHODS[f_method], f_source,
        r(f, 2), r(f_range[0], 2) if f_range else '', r(f_range[1], 2) if f_range else '',
        r(src('PLAT')), r(src('PLONG')), r(src('A95')),
        r(kent['inc']), r(kent['dec']),
        r(kent['Zdec']), r(kent['Zinc']), r(kent['Zeta']),
        r(kent['Edec']), r(kent['Einc']), r(kent['Eta']),
        r(kent_a95_approx(kent['Zeta'], kent['Eta'])),
        _int_or_blank(kent.get('n', '')), comment,
    ]
    return pd.DataFrame([values], columns=KENT_COLUMNS)


def save_kent_summary(summary, filename, output_dir='../data/nordic_summaries/kent'):
    """Save a Kent-pole summary (single-row DataFrame) to a CSV file.

    Companion to ``save_nordic_summary``; ``build_kent_poles.py`` concatenates
    these into ``data/nordic_summaries/kent_poles_combined.csv``.

    Args:
        summary (pd.DataFrame): One-row summary from ``make_kent_summary``.
        filename (str): Output CSV filename, typically the notebook name
            (e.g. '755_Chuar_Group' or '755_Chuar_Group.csv').
        output_dir (str): Directory in which to write the CSV.

    Returns:
        str: The full path to the written CSV file.
    """
    return save_nordic_summary(summary, filename, output_dir=output_dir)
