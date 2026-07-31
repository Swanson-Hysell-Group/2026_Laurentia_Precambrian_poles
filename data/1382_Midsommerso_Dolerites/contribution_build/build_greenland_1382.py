"""Build the Midsommerso Dolerites + Zig-Zag Dal Basalt MagIC contributions (both 1382).

Source: a student MagIC contribution for Marcussen & Abrahamsen (1983),
"Palaeomagnetism of the Proterozoic Zig-Zag Dal Basalt and the Midsommerso
Dolerites, eastern North Greenland," Geophys. J. R. astr. Soc. 73, 367-387
(doi:10.1111/j.1365-246X.1983.tb03321.x), audited against the paper. One source
file carries both units; this script writes a separate contribution into each
unit's data directory.

Audited against Marcussen & Abrahamsen (1983):
  - directions are bedding-corrected (Table 2) -> dir_tilt_correction = 100;
  - demagnetization was AF only (thermal was used only for Curie-point
    determination), remanence carried by titanomagnetite / near-pure magnetite
    -> method_codes LP-DIR-AF:DE-FM (the source's invalid "LP-AF" is corrected;
    no LP-DIR-T);
  - lithology / class set from the text: Zig-Zag Dal = subaerial tholeiitic flood
    basalt (Extrusive / Lava Flow), Midsommerso = quartz-tholeiitic dolerite sills
    and dykes (Intrusive / Sill);
  - single polarity, no field test performed (no reversal / fold / baked-contact /
    conglomerate) -> no structured field-test CV column, pole_reversed_perc 0;
  - the pole is recomputed as the Fisher mean of the site VGPs and reported at the
    positive-latitude antipode (the recalculated-pole convention), reproducing the
    notebook / Nordic export (Midsommerso 10.0N/242.0E, Zig-Zag Dal 12.2N/242.8E);
    the source location rows carried the published southern pole (Marcussen's
    6.9S/62.0E, 12.2S/62.8E);
  - bounding-box fixed to the formation-mean sampling coordinates (Table 2
    footnote): Midsommerso 81.8N/327.8E (32.2 W), Zig-Zag Dal 81.2N/334.8E
    (25.2 W); the source had lat_s = -81.x and a mixed lon_w.

Age: 1382 +/- 2 Ma U-Pb baddeleyite on the associated intrusions (Upton et al.,
2005, Contrib. Mineral. Petrol. 149, 40-56; doi:10.1007/s00410-004-0634-7),
superseding the paper's Rb-Sr ca. 1250 Ma (old decay constant); the underlying
Independence Fjord Group gives a Rb-Sr clay-mineral age of ca. 1380 Ma.
"""
import os, io
from datetime import date
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
MID_OUT = os.path.dirname(HERE)
ZZ_OUT = os.path.abspath(os.path.join(HERE, '../../1382_ZigZag_Dal_Basalt'))
SRC = os.path.join(HERE, 'Marcussen1983_Greenland_magic_source.txt')

CIT = '10.1111/j.1365-246X.1983.tb03321.x'
CIT_LOC = '10.1111/j.1365-246X.1983.tb03321.x:10.1007/s00410-004-0634-7'
AGE, AGE_LOW, AGE_HIGH = 1382, 1380, 1384

# ---- read source (UTF-8; fixes the latin-1 mojibake) -------------------------
data = open(SRC, encoding='utf-8').read()
sites = None
for b in data.split('>>>>>>>>>>'):
    b = b.strip()
    if not b:
        continue
    if 'sites' in b.splitlines()[0]:
        sites = pd.read_csv(io.StringIO('\n'.join(b.splitlines()[1:])), sep='\t')

# recompute signed per-site VGPs from the bedding-corrected directions
sites = sites.drop(columns=[c for c in ['vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm'] if c in sites.columns])


def vgp_row(r):
    plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], r['dir_alpha95'], r['lat'], r['lon'])
    return pd.Series({'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
                      'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1)})


sites = pd.concat([sites, sites.apply(vgp_row, axis=1)], axis=1)

# canonical column fixes applied to every site row
sites['method_codes'] = 'LP-DIR-AF:DE-FM'
sites['citations'] = CIT
sites['dir_tilt_correction'] = 100
sites['age'] = AGE
sites['age_low'] = AGE_LOW
sites['age_high'] = AGE_HIGH
sites = sites.drop(columns=[c for c in ['age_sigma'] if c in sites.columns])


def north_pole(df):
    """Fisher mean of the site VGPs, reported at the positive-latitude antipode."""
    blk = pmag.flip(ipmag.make_di_block(df['vgp_lon'].tolist(), df['vgp_lat'].tolist()), combine=True)
    p = pmag.fisher_mean(blk)
    if p['inc'] < 0:                       # ensure northern hemisphere
        p['dec'] = (p['dec'] + 180) % 360
        p['inc'] = -p['inc']
    return p


def north_dir(df):
    """Fisher mean of the site directions in the normal-polarity (northern-pole) convention."""
    blk = pmag.flip(ipmag.make_di_block(df['dir_dec'].tolist(), df['dir_inc'].tolist()), combine=True)
    d = pmag.fisher_mean(blk)
    if d['inc'] < 0:
        d['dec'] = (d['dec'] + 180) % 360
        d['inc'] = -d['inc']
    return d


UNITS = {
    'Midsommerso Dolerites': dict(
        out=MID_OUT, gclass='Intrusive', gtype='Sill', lith='Diabase',
        lat=81.8, lon=327.8,
        result_name='Midsommerso Dolerites ca. 1382 Ma pole',
        description=(
            "Paleomagnetic pole for the Midsommerso Dolerites, quartz-tholeiitic "
            "dolerite sills, sheets and dykes intruding the Mesoproterozoic "
            "Independence Fjord Group sandstones of eastern North Greenland "
            "(Marcussen & Abrahamsen, 1983). The pole is the Fisher mean of the "
            "bedding-corrected site virtual geomagnetic poles, reported at the "
            "positive-latitude antipode. The single-polarity characteristic "
            "remanence is carried by titanomagnetite to near-pure magnetite "
            "(thermomagnetic Curie points) and was isolated by alternating-field "
            "demagnetization. The pole is used in the reconstructed Laurentia frame "
            "after rotating Greenland to North America. Age from U-Pb baddeleyite "
            "dating of the associated Zig-Zag Dal intrusions (1382 +/- 2 Ma; Upton "
            "et al., 2005), superseding the paper's Rb-Sr ca. 1250 Ma.")),
    'Zig-Zag Dal Basalt Formation': dict(
        out=ZZ_OUT, gclass='Extrusive', gtype='Lava Flow', lith='Basalt',
        lat=81.2, lon=334.8,
        result_name='Zig-Zag Dal Basalt Formation ca. 1382 Ma pole',
        description=(
            "Paleomagnetic pole for the Zig-Zag Dal Basalt Formation, a subaerial "
            "tholeiitic continental flood-basalt succession (up to ca. 1350 m of "
            "aphyric and porphyritic aa flows, with pillow lavas in the basal unit) "
            "conformably overlying the Independence Fjord Group sandstones of "
            "eastern North Greenland (Marcussen & Abrahamsen, 1983). The pole is the "
            "Fisher mean of the bedding-corrected site virtual geomagnetic poles, "
            "reported at the positive-latitude antipode. The single-polarity "
            "characteristic remanence is carried by titanomagnetite to near-pure "
            "magnetite (thermomagnetic Curie points) and was isolated by "
            "alternating-field demagnetization. The pole is used in the "
            "reconstructed Laurentia frame after rotating Greenland to North "
            "America, and is indistinguishable from the pole of the comagmatic "
            "Midsommerso Dolerites. Age from U-Pb baddeleyite dating of the "
            "associated intrusions (1382 +/- 2 Ma; Upton et al., 2005).")),
}


def write_magic(df, kind, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


def validate(path):
    import sys
    root = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
    sys.path.insert(0, os.path.join(root, 'scripts'))
    from validate_magic_contribution import validate_upload_file
    return validate_upload_file(path, tables=['locations', 'sites'])


stamp = date.today().strftime('%d.%b.%Y')
for name, u in UNITS.items():
    s = sites[sites['location'] == name].copy()
    s['geologic_classes'] = u['gclass']
    s['geologic_types'] = u['gtype']
    s['lithologies'] = u['lith']
    s['lat'] = u['lat']
    s['lon'] = u['lon']
    p, d = north_pole(s), north_dir(s)
    lo = pd.DataFrame([{
        'location': name, 'location_type': 'Region', 'result_name': u['result_name'],
        'result_type': 'a', 'result_quality': 'g',
        'sites': ':'.join(map(str, s['site'].tolist())),
        'method_codes': 'LP-DIR-AF:DE-FM:DE-VGP', 'citations': CIT_LOC,
        'geologic_classes': u['gclass'], 'lithologies': u['lith'],
        'lat_s': u['lat'], 'lat_n': u['lat'], 'lon_w': u['lon'], 'lon_e': u['lon'],
        'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
        'dir_tilt_correction': 100,
        'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
        'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1),
        'dir_n_sites': int(d['n']),
        'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
        'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
        'pole_n_sites': int(p['n']), 'pole_reversed_perc': 0,
        'continent_ocean': 'Greenland', 'country': 'Greenland',
        'description': u['description'],
    }])
    write_magic(s, 'sites', os.path.join(u['out'], 'sites.txt'))
    write_magic(lo, 'locations', os.path.join(u['out'], 'locations.txt'))
    base = 'Midsommerso_Dolerites' if 'Midsommerso' in name else 'ZigZag_Dal_Basalt'
    combined = os.path.join(u['out'], f'Marcussen1983_{base}_{stamp}.txt')
    with open(combined, 'w', encoding='utf-8') as f:
        f.write('tab delimited\tlocations\n'); lo.to_csv(f, sep='\t', index=False)
        f.write('>>>>>>>>>>\n')
        f.write('tab delimited\tsites\n'); s.to_csv(f, sep='\t', index=False)
    print(f'-I- {name}: sites {len(s)}, pole {p["inc"]:.1f}/{p["dec"]:.1f} '
          f'A95 {p["alpha95"]:.1f} K {p["k"]:.1f} N {int(p["n"])}; '
          f'dir {d["dec"]:.1f}/{d["inc"]:.1f} a95 {d["alpha95"]:.1f} k {d["k"]:.1f}')
    validate(combined)
