"""Utility functions for Laurentia paleomagnetic pole assessment.

Provides routines for loading and rotating poles into the Laurentia reference
frame, computing mean poles from MagIC site data, evaluating reliability
criteria (Deenen et al., 2011; Meert et al., 2020), and plotting poles in
the context of the Laurentia APWP.
"""

import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import cartopy.crs as ccrs
import pandas as pd

Torsvik2012_poles = pd.read_excel('../data/Torsvik2012.xlsx')
Torsvik2012_Laurentia = Torsvik2012_poles[4:187]

def get_Laurentia_poles(file_name='../data/Kringdalen_w_Laurentia.xlsx', sheet_name='Laurentia'):
    """Loads Laurentia poles and rotates them into a common reference frame.

    Poles from Scotland, Greenland, and Svalbard terranes are rotated into the
    Laurentia reference frame using published Euler poles. Poles from Laurentia
    and Trans-Hudson orogen are kept in their original coordinates. Unrecognized
    terranes receive NaN for rotated coordinates.

    Args:
        file_name (str): Path to the Excel file containing pole data.
            Expected columns include PLAT, PLONG, Terrane, ROCKNAME,
            nominal age, and A95.
        sheet_name (str): Name of the sheet to read from the Excel file.

    Returns:
        pd.DataFrame: Original pole data with added PLAT_rotated and
        PLONG_rotated columns containing poles in the Laurentia reference
        frame.
    """
    Laurentia_poles = pd.read_excel(file_name, sheet_name=sheet_name)

    Euler_poles = {
        'Laurentia-Greenland':      [67.5, -118.5, -13.8],  # [lat, lon, CCW angle]
        'Laurentia-Greenland-Nain': [67.5, -118.5, -13.8],
        'Laurentia-Scotland':       [78.6, 161.9, -31.0],
        'Laurentia-Svalbard':       [-81.0, 125.0, 68.0],
    }

    plat_rot = []
    plon_rot = []
    for _, row in Laurentia_poles.iterrows():
        terrane = row['Terrane']
        if terrane in Euler_poles:
            plat, plon = pmag.pt_rot(Euler_poles[terrane],
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

def get_Laurentia_stricto_poles(file_name='../data/Kringdalen_w_Laurentia.xlsx', sheet_name='Laurentia'):
    """Returns only poles from the Laurentia terrane (sensu stricto).

    Filters the full rotated pole dataset to include only entries where
    Terrane == 'Laurentia', excluding Scotland, Greenland, Svalbard, and
    Trans-Hudson orogen poles.

    Args:
        file_name (str): Path to the Excel file containing pole data.
        sheet_name (str): Name of the sheet to read from the Excel file.

    Returns:
        pd.DataFrame: Subset of poles with Terrane == 'Laurentia', including
        rotated coordinates from ``get_Laurentia_poles``.
    """
    Laurentia_poles = get_Laurentia_poles(file_name=file_name, sheet_name=sheet_name)
    Laurentia_stricto_poles = Laurentia_poles[(Laurentia_poles['Terrane']=='Laurentia')]
    return Laurentia_stricto_poles

def plot_pole_overlap(ROCKNAME, Precambrian_poles, Phanerozoic_poles,
                      pole_plat=None, pole_plon=None, pole_A95=None,
                      pole_age=None):
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
    """

    pole_index = Precambrian_poles.loc[Precambrian_poles['ROCKNAME'] == ROCKNAME].index
    if pole_age is None:
        pole_age = Precambrian_poles['nominal age'].values[pole_index][0]
    if pole_plon is None:
        pole_plon = Precambrian_poles['PLONG_rotated'].values[pole_index][0]
    if pole_plat is None:
        pole_plat = Precambrian_poles['PLAT_rotated'].values[pole_index][0]
    if pole_A95 is None:
        pole_A95 = Precambrian_poles['A95'].values[pole_index][0]

    ax = ipmag.make_mollweide_map(add_land=False, central_longitude=140, figsize=(20,20))

    age_min = 0
    age_max = pole_age

    Precambrian_poles_filtered = Precambrian_poles[Precambrian_poles['nominal age']<=age_max]

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

    ipmag.plot_pole(ax,pole_plon,pole_plat,
                    pole_A95,filled_pole=True,fill_color='green',fill_alpha=0.5)
    ipmag.plot_pole(ax,180+pole_plon,-pole_plat,
                    pole_A95,filled_pole=True,fill_color='green',fill_alpha=0.5)

def plot_apwp_context(Laurentia_poles, pole_plat, pole_plon, pole_A95,
                      age_min=540, age_max=1780, central_longitude=160,
                      figsize=(12, 12)):
    """Plots a pole in the context of the Laurentia Precambrian APWP.

    Shows the Laurentia apparent polar wander path color-coded by age with
    the target pole highlighted in green. Only includes Laurentia and
    Greenland (rotated) poles; excludes Svalbard and Scotland. Uses rotated
    coordinates throughout.

    Args:
        Laurentia_poles (pd.DataFrame): Output of ``get_Laurentia_poles``
            with columns PLONG_rotated, PLAT_rotated, A95, nominal age,
            Terrane, and ROCKNAME.
        pole_plat (float): Latitude of the pole to highlight in degrees.
        pole_plon (float): Longitude of the pole to highlight in degrees.
        pole_A95 (float): A95 of the pole to highlight in degrees.
        age_min (float): Minimum age for filtering in Ma.
        age_max (float): Maximum age for filtering in Ma.
        central_longitude (float): Center longitude for the Mollweide
            projection.
        figsize (tuple): Figure size as (width, height) in inches.

    Returns:
        matplotlib.axes.Axes: The Mollweide map axis.
    """
    ax = ipmag.make_mollweide_map(central_longitude=central_longitude,
                                   figsize=figsize)

    # Exclude Svalbard and Scotland; keep Laurentia, Greenland (rotated),
    # and Trans-Hudson
    excluded_terranes = ('Laurentia-Scotland', 'Laurentia-Svalbard')
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
    
    A95 = this_pole['A95'][0]
    N = this_pole['N'][0]
    B = this_pole['B'][0]
    KD = this_pole['KD'][0]
    
    if N >= 25:
        print('N = ' + str(round(N)) + ' (N ≥ 25; sufficient sample number);')
    else:
        print('N = ' + str(round(N)) + ' (N < 25; insufficient sample number);')  
        
    if B >= 8:
        print('B = ' + str(round(B)) + ' (B ≥ 8; sufficient site number);')
    else:
        print('B = ' + str(round(B)) + ' (B < 8; insufficient site number);') 
        
    if KD >= 70:
        print('K = ' + str(round(KD)) + ' (K ≥ 70; concern about underrepresenting PSV);')
    elif KD >= 10:
        print('K = ' + str(round(KD)) + ' (70 ≥ K ≥ 10; meets PSV criteria);') 
    else: 
        print('K = ' + str(round(KD)) + ' (10 ≥ K; low K, too dispersed);')
        
    Deenen_test(B,A95)

def load_magic_sites(sites_path):
    """Loads a MagIC sites.txt file and splits by tilt correction.

    Reads a tab-delimited MagIC sites table (skipping the header row) and
    returns separate DataFrames for geographic (dir_tilt_correction == 0)
    and tilt-corrected (dir_tilt_correction == 100) coordinates.

    Args:
        sites_path (str): Path to a MagIC-format sites.txt file.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: (sites_geo, sites_tc) DataFrames
        for geographic and tilt-corrected coordinates respectively.
    """
    sites = pd.read_csv(sites_path, sep='\t', skiprows=1)
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
                                     figsize=figsize)
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
    the same site) are collapsed by site name.

    Args:
        sites (pd.DataFrame): Site data with columns ``site``, ``lat``,
            and ``lon`` (longitude in 0–360°).
        zoom_start (int): Initial zoom level for the folium map.
        tiles (str): Folium tile layer name (e.g., 'OpenStreetMap',
            'CartoDB positron').
        color (str): Outline color of the site markers.
        radius (float): Marker radius in pixels.

    Returns:
        folium.Map: Interactive map with a CircleMarker per site, labeled
        with the site name on hover and a popup showing coordinates.
    """
    import folium

    site_locs = sites[['site', 'lat', 'lon']].drop_duplicates(
        subset='site').copy()
    site_locs['lon'] = ((site_locs['lon'] + 180) % 360) - 180

    m = folium.Map(
        location=[site_locs['lat'].mean(), site_locs['lon'].mean()],
        zoom_start=zoom_start,
        tiles=tiles,
    )

    for _, row in site_locs.iterrows():
        folium.CircleMarker(
            location=[row['lat'], row['lon']],
            radius=radius,
            color=color,
            fill=True,
            fill_opacity=0.85,
            popup=folium.Popup(
                f"<b>{row['site']}</b><br>"
                f"{row['lat']:.3f}°N, {row['lon']:.3f}°E",
                max_width=200,
            ),
            tooltip=row['site'],
        ).add_to(m)

    return m