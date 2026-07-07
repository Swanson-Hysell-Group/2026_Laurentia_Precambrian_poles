"""Build the updated East-Central Minnesota Batholith (ECMB) MagIC contribution.

The published contribution earthref.org/MagIC/**20213** (Swanson-Hysell et al.,
2021, *Tectonics*, DOI 10.1029/2021TC006751) archives the sites, samples,
specimens, and measurements of the ECMB study but its ``locations`` table carries
**no location-level pole result**. This script produces an updated contribution
that is **identical to the published 20213 except that the paleomagnetic pole for
the ca. 1779 Ma NE-trending diabase dikes is added to the ``locations`` table**.

Design goals of this rebuild:

- **Single authoritative source.** Everything is read from the locally archived
  full contribution ``../previous_MagIC/magic_contribution_20213.txt``. No network
  download and no separate per-table source copies, so the sites / samples /
  specimens / measurements tables are guaranteed to be *exactly* those of the
  prior contribution (byte-for-byte; their raw text blocks are copied through
  unchanged). This fixes an earlier build that had been assembled from a
  different, less complete contribution (17072) and was missing 5,018
  measurement rows and 5 specimens relative to 20213.
- **Only the locations table changes.** The ``ECMB`` location row gains the
  ``pole_*`` result columns, the pole method codes (incl. the field test), a
  descriptive ``result_name``/``description``, and the refined dated pole age.
- **Canonical MagIC table markers** (``tab delimited<TAB><table>``).
- **No contribution table** in the upload file; MagIC assigns the id, version,
  and contributor on upload (matches the Michipicoten reference build).

Pole provenance (reproduces Swanson-Hysell et al., 2021):
  Fisher mean of the medium-coercivity (mc) VGPs of the NE-trending diabase dikes
  (low-Ti titanomagnetite ChRM, thermally unblocked 515-565 C), keeping site means
  with mc ``dir_alpha95 < 8`` deg. Drops NED3, NED17, NED19, NED21 by the cut;
  excludes the ca. 1096 Ma NW-trending Midcontinent Rift dike NWD1 and the granite
  sites. Result: N = 23 sites, 148 samples, 20.5 degN / 265.8 degE, A95 = 4.5,
  K = 45.6. Directions are geographic (near-vertical dikes; no tilt correction).

Field test: a positive **inverse baked-contact test** -- the ca. 1096 Ma
NW-trending dike NWD1 bakes the NE-trending dike NED17 -- shows the dike remanence
predates ca. 1096 Ma and supports a primary magnetization. MagIC method code
``ST-C-I`` (Inverse contact test).

Outputs:
- ``../locations.txt`` -- source location row + pole result (notebook runtime data).
- ``../Swanson-Hysell2021_updated_ECMB_<DD.Mon.YYYY>.txt`` -- combined upload file.

Run from this directory:
    python build_ECMB_contribution.py
"""

import io
import os
from datetime import date

import pandas as pd
import pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
ECMB_DIR = os.path.dirname(HERE)
SOURCE = os.path.join(ECMB_DIR, 'previous_MagIC', 'magic_contribution_20213.txt')

SEP = '>>>>>>>>>>'
LOCATION = 'ECMB'
ALPHA95_CUT = 8.0          # paper's site-mean a95 cut on the mc component
EXCLUDE_SITES = {'NWD1'}   # ca. 1096 Ma NW-trending Midcontinent Rift dike

# Order of columns written to the locations table (published metadata + pole +
# location-level mean direction + structured field-test results).
LOC_COLS = ['location', 'location_type', 'result_name', 'result_type',
            'result_quality', 'method_codes', 'citations', 'geologic_classes',
            'lithologies', 'lat_s', 'lat_n', 'lon_w', 'lon_e', 'age',
            'age_sigma', 'age_unit', 'dir_tilt_correction', 'dir_dec',
            'dir_inc', 'dir_alpha95', 'dir_k', 'dir_n_sites', 'dir_n_samples',
            'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k', 'pole_n_sites',
            'contact_test', 'sites', 'description']


def read_blocks(path):
    """Split a combined MagIC file into ordered (table_name, marker, body) blocks.

    ``body`` is the raw text of the block *after* the marker line (column header
    plus data rows), preserved verbatim so tables can be copied through unchanged.
    """
    with open(path) as f:
        content = f.read()
    blocks = []
    for chunk in content.split(SEP):
        chunk = chunk.strip('\n')
        if not chunk:
            continue
        marker, _, body = chunk.partition('\n')
        table = marker.split('\t')[-1].strip()
        blocks.append((table, marker, body))
    return blocks


def block_to_df(body):
    return pd.read_csv(io.StringIO(body), sep='\t', dtype=str).fillna('')


def compute_pole(sites):
    """Return (pole, meandir, n_samples, site_list) for the mc NE-dike sites.

    ``pole`` is the Fisher mean of the site VGPs; ``meandir`` is the Fisher mean
    of the site mean directions in geographic coordinates (the pole's underlying
    location-level mean direction). The two carry the same N sites.
    """
    a95 = pd.to_numeric(sites['dir_alpha95'], errors='coerce')
    pole_sites = sites[(sites['dir_comp_name'] == 'mc') &
                       (~sites['site'].isin(EXCLUDE_SITES)) &
                       (sites['result_quality'] == 'g') &
                       (a95 < ALPHA95_CUT)].copy()
    pole = ipmag.fisher_mean(pd.to_numeric(pole_sites['vgp_lon']).tolist(),
                             pd.to_numeric(pole_sites['vgp_lat']).tolist())
    meandir = ipmag.fisher_mean(dec=pd.to_numeric(pole_sites['dir_dec']).tolist(),
                                inc=pd.to_numeric(pole_sites['dir_inc']).tolist())
    n_samples = int(pd.to_numeric(pole_sites['dir_n_samples']).sum())
    site_list = ':'.join(pole_sites['site'].tolist())
    return pole, meandir, n_samples, site_list


def build_location_row(src_loc_df, pole, meandir, n_samples, site_list):
    """Return the ECMB location row (dict) with the pole result added."""
    row = src_loc_df.iloc[0].to_dict()

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
        'NW-trending dike NWD1 bakes NE-trending dike NED17; MagIC method code '
        'ST-C-I) shows the dike remanence predates ca. 1096 Ma and supports a '
        'primary magnetization. Age of 1779.1 +/- 2.3 Ma (95% CI) from new U-Pb '
        'dates bracketing the dikes between the St. Cloud Granite (1781.44 +/- '
        '0.51 Ma) they intrude and the Richmond Granite (1776.76 +/- 0.49 Ma) '
        'they do not.')

    row.update({
        'age': '1779.1',              # refined to the dated pole age (95% CI below)
        'age_sigma': '2.3',           # 95% CI half-width (see description)
        'age_unit': 'Ma',
        'citations': '10.1029/2021TC006751',
        'description': description,
        'dir_tilt_correction': '0',
        'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM:DE-VGP:ST-C-I',
        # location-level mean direction (geographic; underlies the pole)
        'dir_dec': f'{meandir["dec"]:.1f}',
        'dir_inc': f'{meandir["inc"]:.1f}',
        'dir_alpha95': f'{meandir["alpha95"]:.1f}',
        'dir_k': f'{meandir["k"]:.1f}',
        'dir_n_sites': str(int(meandir['n'])),
        'dir_n_samples': str(n_samples),
        # pole (Fisher mean of site VGPs)
        'pole_lat': f'{pole["inc"]:.1f}',
        'pole_lon': f'{pole["dec"]:.1f}',
        'pole_alpha95': f'{pole["alpha95"]:.1f}',
        'pole_k': f'{pole["k"]:.1f}',
        'pole_n_sites': str(int(pole['n'])),
        # structured field-test result (controlled vocabulary); only the test
        # actually performed is recorded -- no columns for tests not done
        'contact_test': 'IC+',        # positive inverse contact test (NWD1 bakes NED17)
        'result_name': 'East-Central Minnesota Batholith ca. 1779 Ma pole',
        'result_quality': 'g',
        'result_type': 'a',
        'sites': site_list,
    })
    # drop bracketing-age columns if present; the pole carries a dated age
    for col in ('age_low', 'age_high'):
        row.pop(col, None)
    return row


def write_table(fh, table, df, cols):
    fh.write(f'tab delimited\t{table}\n')
    df.to_csv(fh, sep='\t', index=False, columns=cols)


def main():
    blocks = read_blocks(SOURCE)
    by_name = {t: (m, b) for t, m, b in blocks}

    sites = block_to_df(by_name['sites'][1])
    src_loc = block_to_df(by_name['locations'][1])

    pole, meandir, n_samples, site_list = compute_pole(sites)
    print(f'-I- pole from 20213: N={int(pole["n"])} lon={pole["dec"]:.1f} '
          f'lat={pole["inc"]:.1f} A95={pole["alpha95"]:.1f} K={pole["k"]:.1f} '
          f'samples={n_samples}')
    print('-I- published target : N=23 lon=265.8 lat=20.5 A95=4.5 K=45.6 samples=148')
    print(f'-I- mean direction : dec={meandir["dec"]:.1f} inc={meandir["inc"]:.1f} '
          f'a95={meandir["alpha95"]:.1f} k={meandir["k"]:.1f}')

    loc_row = build_location_row(src_loc, pole, meandir, n_samples, site_list)
    loc_df = pd.DataFrame([loc_row])

    # ---- notebook runtime locations.txt ------------------------------------
    out_loc = os.path.join(ECMB_DIR, 'locations.txt')
    with open(out_loc, 'w') as f:
        write_table(f, 'locations', loc_df, LOC_COLS)
    print(f'-I- wrote {out_loc}')

    # ---- combined upload file ----------------------------------------------
    # locations regenerated; sites/samples/specimens/measurements copied verbatim
    # from 20213; contribution table omitted (MagIC assigns it on upload).
    stamp = date.today().strftime('%d.%b.%Y')
    out_upload = os.path.join(
        ECMB_DIR, f'Swanson-Hysell2021_updated_ECMB_{stamp}.txt')
    passthrough = ['sites', 'samples', 'specimens', 'measurements']
    with open(out_upload, 'w') as f:
        write_table(f, 'locations', loc_df, LOC_COLS)
        for table in passthrough:
            marker, body = by_name[table]
            f.write(SEP + '\n')
            f.write(f'tab delimited\t{table}\n')  # canonical marker
            f.write(body)
            if not body.endswith('\n'):
                f.write('\n')
    print(f'-I- wrote {out_upload}')
    print('-I- tables in upload: locations (with pole) + '
          'sites/samples/specimens/measurements verbatim from 20213')


if __name__ == '__main__':
    main()
