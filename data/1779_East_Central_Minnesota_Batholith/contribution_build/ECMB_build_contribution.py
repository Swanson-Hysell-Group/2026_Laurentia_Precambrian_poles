"""Assemble the East-Central Minnesota Batholith (ECMB) pole contribution.

This script enhances the published Swanson-Hysell et al. (2021) MagIC
contribution **17072** (DOI 10.1029/2021TC006751) by adding the location-level
paleomagnetic pole result for the ca. 1779 Ma northeast-trending diabase dikes
of the ECMB. The published contribution archives the full site / sample /
specimen / measurement data but its ``locations`` table carries no ``pole_*``
result, so "enhancing" here means *adding* the pole (per the conventions in
``plan.md`` and the Lake Shore Traps / Michipicoten builds).

Provenance of the pole (reproduces Swanson-Hysell et al., 2021; see their
``ECMB_pmag_analysis.ipynb``):

- Source = the site table of contribution 17072 (``sites_17072_source.txt``),
  which carries multiple demagnetization components per site
  (lc = low-coercivity overprint, mc = medium-coercivity ChRM, hc =
  high-coercivity Midcontinent-Rift overprint, plus thermal lt/mt/ht fits).
- The pole is built from the **medium-coercivity (mc) component** — the
  characteristic remanence held by low-Ti titanomagnetite (thermally unblocked
  515-565 deg C) — of the **northeast-trending diabase dikes (NED sites)**.
- The one northwest-trending dike (NWD1, a ca. 1096 Ma Midcontinent Rift dike
  used for the baked-contact test) is excluded, and the granite/other sites
  (RFG1, SCG1) carry no mc component.
- Site means with mc ``dir_alpha95 >= 8 deg`` (NED3, NED17, NED19, NED21) are
  excluded following the paper's site-mean uncertainty cut, leaving **23 sites**.

Outputs (written one level up, the runtime data for the notebook):

- ``../sites.txt``     — the full multi-component site table (unchanged from
                         17072), so the contribution stays complete.
- ``../locations.txt`` — the source location row with the ECMB pole result
                         added (``pole_*`` columns, baked-contact ``ST-C`` code,
                         selection documented in ``description``).

The ECMB dikes are near-vertical intrusions reported in geographic coordinates
(no structural/tilt correction; ``dir_tilt_correction = 0``).

Run from this directory with the project interpreter, e.g.:
    python ECMB_build_contribution.py
"""

import os
import pandas as pd
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.dirname(HERE)

# --- mc-component site-mean uncertainty cut and excluded NW dike ---------------
ALPHA95_CUT = 8.0          # paper's site-mean a95 cut on the mc component
EXCLUDE_SITES = {'NWD1'}   # ca. 1096 Ma NW-trending Midcontinent Rift dike


def read_magic_table(path):
    """Read a MagIC v3 table, returning (header_line, DataFrame)."""
    with open(path) as f:
        header_line = f.readline().rstrip('\n')
    df = pd.read_csv(path, sep='\t', skiprows=1, dtype=str)
    return header_line, df


def main():
    src_sites = os.path.join(HERE, 'sites_17072_source.txt')
    src_locs = os.path.join(HERE, 'locations_17072_source.txt')

    sites_header, sites = read_magic_table(src_sites)
    locs_header, locs = read_magic_table(src_locs)

    # ----- write the full site table through unchanged --------------------------
    out_sites = os.path.join(OUT_DIR, 'sites.txt')
    with open(out_sites, 'w') as f:
        f.write(sites_header + '\n')
    sites.to_csv(out_sites, sep='\t', index=False, mode='a')
    print(f'-I- wrote {out_sites} ({len(sites)} site rows, all components)')

    # ----- select the pole sites (mc component of the NE-trending dikes) --------
    a95 = pd.to_numeric(sites['dir_alpha95'], errors='coerce')
    pole_sites = sites[(sites['dir_comp_name'] == 'mc') &
                       (~sites['site'].isin(EXCLUDE_SITES)) &
                       (sites['result_quality'] == 'g') &
                       (a95 < ALPHA95_CUT)].copy()
    pole_sites['vgp_lon'] = pd.to_numeric(pole_sites['vgp_lon'])
    pole_sites['vgp_lat'] = pd.to_numeric(pole_sites['vgp_lat'])
    n_samples = int(pd.to_numeric(pole_sites['dir_n_samples']).sum())

    pole = ipmag.fisher_mean(pole_sites['vgp_lon'].tolist(),
                             pole_sites['vgp_lat'].tolist())
    site_list = ':'.join(pole_sites['site'].tolist())

    print(f'-I- {len(pole_sites)} mc sites in pole '
          f'(a95 < {ALPHA95_CUT} deg, excluding {sorted(EXCLUDE_SITES)})')
    print(f'-I- pole: lon={pole["dec"]:.2f} lat={pole["inc"]:.2f} '
          f'N={pole["n"]} K={pole["k"]:.1f} A95={pole["alpha95"]:.2f} '
          f'(samples={n_samples})')
    print('-I- published (Swanson-Hysell et al., 2021): '
          'lon=265.8 lat=20.4 N=23 K=45.6 A95=4.5')

    # ----- add the pole result to the location row ------------------------------
    description = (
        'Paleomagnetic pole for the ca. 1779 Ma northeast-trending diabase '
        'dikes of the East-Central Minnesota Batholith (Swanson-Hysell et al., '
        '2021). Built from the medium-coercivity (mc) characteristic remanence '
        '(low-Ti titanomagnetite; thermally unblocked 515-565 C) of the '
        'NE-trending dikes, using the 23 site means with mc dir_alpha95 < 8 deg '
        '(NED3, NED17, NED19, NED21 excluded by the cut; the NW-trending ca. '
        '1096 Ma Midcontinent Rift dike NWD1 and the granite sites excluded). '
        'Directions are in geographic coordinates (near-vertical dikes; no tilt '
        'correction). A positive inverse baked-contact test (the ca. 1096 Ma '
        'NW-trending dike NWD1 bakes NE-trending dike NED17; method code ST-C) '
        'shows the dike remanence predates ca. 1096 Ma and supports a primary '
        'magnetization. Age of 1779.1 +/- 2.3 Ma (95% CI) from new U-Pb dates '
        'bracketing the dikes between the St. Cloud Granite (1781.44 +/- 0.51 '
        'Ma) they intrude and the Richmond Granite (1776.76 +/- 0.49 Ma) they '
        'do not.')

    pole_fields = {
        'age': '1779.1',
        'age_sigma': '2.3',
        'age_unit': 'Ma',
        'citations': '10.1029/2021TC006751',
        'description': description,
        'dir_tilt_correction': '0',
        'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP:ST-C',
        'pole_alpha95': f'{pole["alpha95"]:.1f}',
        'pole_k': f'{pole["k"]:.1f}',
        'pole_lat': f'{pole["inc"]:.1f}',
        'pole_lon': f'{pole["dec"]:.1f}',
        'pole_n_sites': str(int(pole['n'])),
        'pole_n_samples': str(n_samples),
        'result_name': 'East-Central Minnesota Batholith ca. 1779 Ma pole',
        'result_quality': 'g',
        'result_type': 'a',
        'sites': site_list,
    }

    loc = locs.iloc[0].to_dict()
    # the source location row carries age 1779; refine to the dated pole age
    for col in ('age_high', 'age_low'):
        loc.pop(col, None)
    loc.update(pole_fields)

    out_locs = os.path.join(OUT_DIR, 'locations.txt')
    loc_df = pd.DataFrame([loc])
    with open(out_locs, 'w') as f:
        f.write('tab \tlocations\n')
    loc_df.to_csv(out_locs, sep='\t', index=False, mode='a')
    print(f'-I- wrote {out_locs} with the ECMB pole result added')


if __name__ == '__main__':
    main()
