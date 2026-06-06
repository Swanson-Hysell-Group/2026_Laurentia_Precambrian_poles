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
      the statistical thresholds N >= 25 samples, 10 <= K <= 70, B >= 8 sites
      with at least 3 samples per site.

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
    c_minsamp = min_per_site >= 3
    c = bool(deenen_pass and c_n and c_k and c_b and c_minsamp)

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
        print(f'        min {min_per_site} samples/site (>= 3): '
              f'{"PASS" if c_minsamp else "FAIL"}')
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
            label='Our pole')
    ax.set_xlabel('Number of sites (N)', fontsize=12)
    ax.set_ylabel(r'$A_{95}$ (°)', fontsize=12)
    ax.set_title('Deenen et al. (2011) test', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(5, n_max)
    ax.set_ylim(0, 25)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    return ax

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
    return ipmag.fishqq(di_block=vgp_block, data_type='poles')


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


def make_nordic_summary(terrane, 
                        rockname, 
                        sites,
                        dir_mean,
                        pole_mean,
                        study_lon,
                        study_lat,
                        component_comment='',
                        tests='',
                        f_factor=1,
                        pole_mean_unflattened=None,
                        R1=None,
                        R2=None,
                        R3=None,
                        R4=None,
                        R5=None,
                        R6=None,
                        R7=None,
                        Grade=None,
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
                        COMMENT=''
                        ):
    nordic_dict = {}
    nordic_dict['terrane'] = terrane
    nordic_dict['rockname'] = rockname
    nordic_dict['component_comment'] = component_comment
    nordic_dict['tests'] = tests
    if sites['dir_tilt_correction'].nunique() != 1:
        warnings.warn(
            "Multiple tilt correction values found in sites data; "
            "cannot determine single tilt correction for summary, choosing first value."
        )    
    nordic_dict['tilt'] = sites['dir_tilt_correction'].iloc[0]
    nordic_dict['study_lat'] = study_lat
    nordic_dict['study_lon'] = study_lon
    nordic_dict['site_n'] = pole_mean['n']
    nordic_dict['sample_n'] = sites['dir_n_samples'].sum()
    nordic_dict['dir_dec_mean'] = dir_mean['dec']
    nordic_dict['dir_inc_mean'] = dir_mean['inc']
    nordic_dict['dir_inc_mean_abs'] = np.abs(dir_mean['inc'])
    nordic_dict['dir_k'] = dir_mean['k']
    nordic_dict['dir_alpha_95'] = dir_mean['alpha95']
    nordic_dict['pole_lat'] = pole_mean['dec']
    nordic_dict['pole_lon'] = pole_mean['inc']
    nordic_dict['pole_dp'] = ''
    nordic_dict['pole_dm'] = ''
    nordic_dict['pole_A95'] = pole_mean['alpha95']
    nordic_dict['f_factor'] = f_factor
    nordic_dict['inc_f'] = ipmag.unsquish([dir_mean['inc']], f_factor)[0]

    if pole_mean_unflattened is not None:
        nordic_dict['pole_lat_unflattened'] = pole_mean_unflattened['dec']
        nordic_dict['pole_lon_unflattened'] = pole_mean_unflattened['inc']
        nordic_dict['pole_dp_unflattened'] = ''
        nordic_dict['pole_dm_unflattened'] = ''
        nordic_dict['pole_A95_unflattened'] = pole_mean_unflattened['alpha95']
    else:
        nordic_dict['pole_lat_unflattened'] = pole_mean['dec']
        nordic_dict['pole_lon_unflattened'] = pole_mean['inc']
        nordic_dict['pole_dp_unflattened'] = ''
        nordic_dict['pole_dm_unflattened'] = ''
        nordic_dict['pole_A95_unflattened'] = pole_mean['alpha95']
    
    # check any other optional parameter if None then throw error
    optional_params = {
        'R1': R1,
        'R2': R2,
        'R3': R3,
        'R4': R4,
        'R5': R5,
        'R6': R6,
        'R7': R7,
        'Grade': Grade,
        'nominal_age': nominal_age,
        'lomagage': lomagage,
        'himagage': himagage,
        'REF_method': REF_method,
        'POLE_AUTHORS': POLE_AUTHORS,
        'YEAR': YEAR,
        'JOURNAL': JOURNAL,
        'VOLUME': VOLUME,
        'VPAGES': VPAGES,
        'TITLE': TITLE,
        'COMMENT': COMMENT
    }
    for param_name, param_value in optional_params.items():
        if param_value is None:
            raise ValueError(f"Optional parameter '{param_name}' is required but was not provided.")
        nordic_dict[param_name] = param_value
    return nordic_dict


def save_nordic_summary(summary, filename, output_dir='../data/nordic_summaries'):
    """Save a Nordic summary dictionary to a CSV file.

    The summary is written as a single-row CSV with the dictionary keys as
    column headers, into ``output_dir`` (created if it does not yet exist).

    Args:
        summary (dict): Nordic summary dictionary from make_nordic_summary.
        filename (str): Output CSV filename, typically the notebook name
            (e.g. '780_Gunbarrel' or '780_Gunbarrel.csv').
        output_dir (str): Directory in which to write the CSV. Defaults to
            the data/nordic_summaries folder relative to a pole notebook.

    Returns:
        str: The full path to the written CSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    if not filename.endswith('.csv'):
        filename += '.csv'
    output_path = os.path.join(output_dir, filename)
    pd.DataFrame([summary]).to_csv(output_path, index=False)
    return output_path
