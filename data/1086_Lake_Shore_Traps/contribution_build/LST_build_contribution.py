"""
Build the enhanced MagIC 3.0 contribution for the Lake Shore Traps pole.

This one-time build produced the contribution now published at
earthref.org/MagIC/20696. The pole notebook loads the data directly from that
published contribution (pt.fetch_magic_contribution('20696', ...)); this script
is retained as the documented record of how the contribution was assembled and
the corrections that were applied.

Approach (study-level MagIC contributions, enhanced):
    MagIC contributions are organized by study. The two site-level directional
    studies for the Lake Shore Traps are already in MagIC as legacy
    contributions, downloaded here as the source of truth:
        - magic_contribution_16334_Diehl_Haig_1994.txt  (Diehl & Haig, 1994)
        - magic_contribution_16335_Kulakov_2013.txt      (Kulakov et al., 2013)
    (retrieved with ipmag.download_magic_from_doi / the EarthRef API
     https://api.earthref.org/v1/MagIC/download?doi=<DOI>).

    Kulakov et al. (2013) is a follow-up to Diehl & Haig (1994), so this script
    merges the Diehl & Haig sites INTO the Kulakov 2013 contribution (reference
    10.1139/cjes-2013-0003) as a single enhanced contribution that carries the
    combined Lake Shore Traps mean pole as its location-level result. The pole
    notebook then pulls solely from this one contribution.

    Per-site provenance is preserved through the site-level ``citations`` field
    (Diehl sites cite 10.1139/e94-034; Kulakov sites cite
    10.1139/cjes-2013-0003).

Enhancements applied to the legacy contributions (site + location level only;
the studies have no measurement-level data available):
    - citations set to DOIs (the legacy tables had "This study");
    - method_codes set to LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM, confirmed from both
      papers (AF + thermal demagnetization, principal-component analysis after
      Kirschvink (1980), Fisher means) — replacing the garbled legacy codes;
    - site-level geologic_classes (Igneous), geologic_types (Lava Flow),
      lithologies (Basalt), and age added;
    - vgp_dp / vgp_dm computed with pmag.dia_vgp;
    - longitudes converted to the MagIC 0-360 East convention;
    - a proper location-level pole result (pole_lat/lon/alpha95/k/n_sites and
      the contributing site list) written to locations.txt.

Age source: Fairchild, L. M., Swanson-Hysell, N. L., Ramezani, J., Sprain,
    C. J., & Bowring, S. A. (2017). The end of Midcontinent Rift magmatism and
    the paleogeography of Laurentia. Lithosphere, 9(1), 117-133.
    https://doi.org/10.1130/L580.1  (U-Pb 206Pb/238U weighted-mean date of
    1085.57 +/- 0.5 Ma on an andesitic flow within the Lake Shore Traps).

DUPLICATE FLOW: site LST28 was measured by both studies. The Diehl & Haig
    (1994) measurement is used in the pole (chronological priority); the Kulakov
    et al. (2013) re-measurement is retained in the contribution as a separate
    site ``LST28_K2013`` flagged ``result_quality = 'b'`` and excluded from the
    pole to avoid double-counting the flow. The pole therefore uses 50 unique
    site means (51 site rows total: 50 good + 1 excluded duplicate).

MagIC data model: v3.0 (https://www2.earthref.org/MagIC/data-models/3.0)
"""

from pathlib import Path
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

LOCATION = 'Lake Shore Traps'

# study -> (downloaded MagIC contribution file, publication DOI)
STUDIES = [
    ('Kulakov_2013', 'magic_contribution_16335_Kulakov_2013.txt',
     '10.1139/cjes-2013-0003'),
    ('Diehl_Haig_1994', 'magic_contribution_16334_Diehl_Haig_1994.txt',
     '10.1139/e94-034'),
]
STUDY_LABEL = {
    'Kulakov_2013': 'Kulakov et al. (2013)',
    'Diehl_Haig_1994': 'Diehl & Haig (1994)',
}
MASTER_DOI = '10.1139/cjes-2013-0003'   # Kulakov 2013 is the master contribution
AGE_CITATION = '10.1130/L580.1'         # Fairchild et al. (2017)
CONGLOMERATE_CITATION = '10.1139/e81-053'  # Palmer, Halls & Pesonen (1981), ST-G

METHOD_CODES = 'LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM'

# Radiometric age: use age + age_sigma (not age_low/age_high) per the MagIC
# data model and resources/magic_data_entry_guide.md (use one or the other).
AGE, AGE_SIGMA = '1085.57', '0.5'

# Correction: the downloaded MagIC contribution 16335 has the wrong Fisher k for
# the ten Silver Island (SI) sites — its SI k column was mis-filled with the
# k values of the adjacent LST1a-LST10a rows. These are the correct k values
# from Kulakov et al. (2013), Table 1 (verified against the paper PDF); all other
# SI fields (Dg/Ig/Ds/Is, alpha95, N, VGP) in 16335 are correct.
SI_K_PAPER = {'SI1': 282, 'SI2': 153, 'SI3': 281, 'SI4': 1025, 'SI5': 109,
              'SI6': 738, 'SI7': 221, 'SI8': 320, 'SI9': 574, 'SI10': 225}

# The duplicated flow LST28 (measured by both studies): keep Diehl & Haig's in
# the pole; retain Kulakov's re-measurement as a separate, excluded site.
DUP = ('Kulakov_2013', 'LST28')
DUP_RENAME = 'LST28_K2013'

SITE_COLS = [
    'site', 'location', 'result_type', 'result_quality', 'method_codes',
    'citations', 'geologic_classes', 'geologic_types', 'lithologies',
    'lat', 'lon', 'age', 'age_sigma', 'age_unit',
    'bed_dip', 'bed_dip_direction', 'dir_tilt_correction',
    'dir_dec', 'dir_inc', 'dir_k', 'dir_alpha95', 'dir_n_samples',
    'vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm', 'description',
]

LOC_COLS = [
    'location', 'location_type', 'result_name', 'result_type',
    'result_quality', 'method_codes', 'citations', 'geologic_classes',
    'lithologies', 'lat_s', 'lat_n', 'lon_w', 'lon_e',
    'age', 'age_sigma', 'age_unit',
    'dir_tilt_correction',
    'pole_lat', 'pole_lon', 'pole_alpha95', 'pole_k', 'pole_n_sites',
    'sites', 'description',
]


def parse_sites_table(path):
    """Return the rows of the 'sites' table of a MagIC contribution file."""
    text = Path(path).read_text(encoding='utf-8-sig')
    blocks = text.split('>>>>>>>>>>')
    for block in blocks:
        lines = [ln for ln in block.splitlines() if ln.strip() != '']
        if not lines:
            continue
        if lines[0].lower().startswith('tab') and lines[0].split('\t')[-1].strip() == 'sites':
            header = lines[1].split('\t')
            return [dict(zip(header, ln.split('\t'))) for ln in lines[2:]]
    raise ValueError(f'no sites table found in {path}')


def lon_east(lon):
    return float(lon) % 360


def main():
    d = Path(__file__).parent
    site_rows = []
    for study, fname, doi in STUDIES:
        for r in parse_sites_table(d / fname):
            site = r['site'].strip().replace(' ', '')
            quality = 'g'
            if (study, site) == DUP:
                # Kulakov's re-measurement of flow LST28: keep it in the
                # contribution but rename and flag it excluded from the pole.
                site = DUP_RENAME
                quality = 'b'
            dec, inc = float(r['dir_dec']), float(r['dir_inc'])
            a95 = float(r['dir_alpha95'])
            lat, lon = float(r['lat']), lon_east(r['lon'])
            _, _, dp, dm = pmag.dia_vgp(dec, inc, a95, lat, lon)
            dir_k = str(SI_K_PAPER[site]) if site in SI_K_PAPER else r['dir_k'].strip()
            member = r.get('member', '').strip()
            desc = f'{member}; {STUDY_LABEL[study]}' if member else STUDY_LABEL[study]
            if site == 'LST28':
                desc += (' — flow also measured by Kulakov et al. (2013), '
                         'retained in this contribution as site LST28_K2013 '
                         '(result_quality=b) but excluded from the pole to avoid '
                         'double-counting; this Diehl & Haig (1994) measurement '
                         'is the one used.')
            elif site == DUP_RENAME:
                desc += (' — Kulakov et al. (2013) re-measurement of flow LST28 '
                         '(Diehl & Haig, 1994); flagged result_quality=b and '
                         'excluded from the pole to avoid double-counting the '
                         'flow (the Diehl & Haig measurement is used).')
            site_rows.append({
                'site': site, 'location': LOCATION,
                'result_type': 'i', 'result_quality': quality,
                'method_codes': METHOD_CODES, 'citations': doi,
                'geologic_classes': 'Igneous', 'geologic_types': 'Lava Flow',
                'lithologies': 'Basalt',
                'lat': f'{lat:.3f}', 'lon': f'{lon:.3f}',
                'age': AGE, 'age_sigma': AGE_SIGMA, 'age_unit': 'Ma',
                'bed_dip': r.get('bed_dip', '').strip(),
                'bed_dip_direction': r.get('bed_dip_direction', '').strip(),
                'dir_tilt_correction': '100',
                'dir_dec': f'{dec:.1f}', 'dir_inc': f'{inc:.1f}',
                'dir_k': dir_k, 'dir_alpha95': f'{a95:.1f}',
                'dir_n_samples': r['dir_n_samples'].strip(),
                'vgp_lat': r['vgp_lat'].strip(), 'vgp_lon': r['vgp_lon'].strip(),
                'vgp_dp': f'{dp:.1f}', 'vgp_dm': f'{dm:.1f}',
                'description': desc,
            })

    sites_path = d / 'sites.txt'
    with open(sites_path, 'w') as f:
        f.write('tab delimited\tsites\n')
        f.write('\t'.join(SITE_COLS) + '\n')
        for r in site_rows:
            f.write('\t'.join(str(r[c]) for c in SITE_COLS) + '\n')
    print(f'Wrote {sites_path.name}: {len(site_rows)} site rows')

    write_locations(d, site_rows)


def write_locations(d, site_rows):
    # the pole uses only accepted (result_quality 'g') sites; the excluded
    # duplicate (LST28_K2013, 'b') is in sites.txt but not in the mean
    pole_rows = [r for r in site_rows if r['result_quality'] == 'g']
    block = ipmag.make_di_block([float(r['vgp_lon']) for r in pole_rows],
                                [float(r['vgp_lat']) for r in pole_rows])
    block = pmag.flip(block, combine=True)
    pole = ipmag.fisher_mean(di_block=block)

    lats = [float(r['lat']) for r in pole_rows]
    lons = [float(r['lon']) for r in pole_rows]
    # Kulakov (master) DOI first, then Diehl, the age citation, and the
    # conglomerate field-test citation (ST-G)
    citations = ':'.join([MASTER_DOI, '10.1139/e94-034', AGE_CITATION,
                          CONGLOMERATE_CITATION])

    loc = {
        'location': LOCATION, 'location_type': 'Outcrop',
        'result_name': 'Lake Shore Traps ca. 1085 Ma pole',
        'result_type': 'a', 'result_quality': 'g',
        # ST-G = conglomerate test (positive); the passed test is conveyed by the
        # method code on the accepted (result_quality 'g') pole plus the note in
        # the description. See resources/field_test_codes.md and MagIC ST- codes.
        'method_codes': METHOD_CODES + ':DE-VGP:ST-G', 'citations': citations,
        'geologic_classes': 'Igneous', 'lithologies': 'Basalt',
        'lat_s': f'{min(lats):.3f}', 'lat_n': f'{max(lats):.3f}',
        'lon_w': f'{min(lons):.3f}', 'lon_e': f'{max(lons):.3f}',
        'age': AGE, 'age_sigma': AGE_SIGMA, 'age_unit': 'Ma',
        'dir_tilt_correction': '100',
        'pole_lat': f'{pole["inc"]:.1f}', 'pole_lon': f'{pole["dec"]:.1f}',
        'pole_alpha95': f'{pole["alpha95"]:.1f}', 'pole_k': f'{pole["k"]:.1f}',
        'pole_n_sites': str(pole['n']),
        'sites': ':'.join(r['site'] for r in pole_rows),
        'description': (
            'Lake Shore Traps mean pole built on the Kulakov et al. (2013) '
            'contribution with the Diehl & Haig (1994) sites merged in '
            '(duplicated flow LST28 counted once; N=50). Age constrained by a '
            'U-Pb (206Pb/238U, CA-ID-TIMS) weighted-mean date of 1085.57 +/- '
            '0.5 Ma on an andesitic flow within the Lake Shore Traps (Fairchild '
            'et al., 2017). Positive conglomerate test on Keweenawan '
            'conglomerates (Palmer, Halls & Pesonen, 1981; method code ST-G) '
            'supports a primary magnetization.'),
    }
    loc_path = d / 'locations.txt'
    with open(loc_path, 'w') as f:
        f.write('tab delimited\tlocations\n')
        f.write('\t'.join(LOC_COLS) + '\n')
        f.write('\t'.join(str(loc[c]) for c in LOC_COLS) + '\n')
    print(f'Wrote {loc_path.name}: 1 pole row')
    print(f'  {loc["location"]}: {loc["pole_lat"]}°N, {loc["pole_lon"]}°E, '
          f'A95={loc["pole_alpha95"]}°, N={loc["pole_n_sites"]}, k={loc["pole_k"]}')


def combine_and_validate(d):
    result = ipmag.upload_magic(dir_path=str(d), input_dir_path=str(d))
    print(f'Created {result[0]}' if result[0] else f'Validation issue: {result[1]}')


if __name__ == '__main__':
    main()
    combine_and_validate(Path(__file__).parent)
