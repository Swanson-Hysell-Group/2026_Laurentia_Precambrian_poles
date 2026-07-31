"""Build the central North Greenland (Victoria Fjord) dolerite dyke pole (1382_Victoria).

Source: a student MagIC contribution for Abrahamsen & Van der Voo (1987),
"Palaeomagnetism of middle Proterozoic (c. 1.25 Ga) dykes from central North
Greenland," Geophys. J. R. astr. Soc. 91, 597-611
(doi:10.1111/j.1365-246X.1987.tb01660.x), audited against the paper. The build
reads a snapshot of the audited site table (Abrahamsen1987_Victoria_sites_source.txt);
the original student spreadsheet lived off-machine.

Audited against Abrahamsen & Van der Voo (1987), Table 1:
  - D6 (site 269) dir_dec set to 267 (the paper value; a prior build carried 275,
    which pushed the 9-dyke mean off the paper's stated D=265.2 to 266.1). With
    D6=267 the recalculated mean direction (265.2/25.2) reproduces the paper.
  - D10 (site 279) is the anomalous NE direction (compass-oriented under overcast);
    excluded from the pole (result_quality 'b').
  - method_codes: AF + thermal demagnetization, characteristic directions read from
    Zijderveld plots as stable end points (the paper predates routine PCA) ->
    LP-DIR-AF:LP-DIR-T:DE-FM (the prior build's DE-BFL best-fit-line code is
    dropped as no principal-component fitting was done).
  - lithology / class from the text: dolerite dykes intruding Archaean gneiss ->
    Intrusive / Volcanic Dike / Diabase. Near-vertical dykes -> in-situ
    (dir_tilt_correction 0).
  - a detailed, positive baked-contact test (baked gneiss and quartz-diorite carry
    the dyke direction; unbaked gneiss at 60 m retains its pre-dyke direction) ->
    ST-C, contact_test C+. Single polarity within the swarm -> pole_reversed_perc 0.

Age: adopted as 1382 +/- 2 Ma by correlation with the U-Pb baddeleyite-dated
Zig-Zag Dal / Midsommerso magmatism (Upton et al., 2005); the paper itself infers
ca. 1250 Ma from Rb-Sr on the correlative eastern North Greenland intrusives -- a
genuine age uncertainty.
"""
import os, io
from datetime import date
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Abrahamsen1987_Victoria_sites_source.txt')

CIT = '10.1111/j.1365-246X.1987.tb01660.x'
CIT_LOC = '10.1111/j.1365-246X.1987.tb01660.x:10.1007/s00410-004-0634-7'
AGE, AGE_LOW, AGE_HIGH = 1382, 1380, 1384
D6_DEC = 267.0        # paper value (Table 1, site 269); a prior build used 275

# ---- read snapshot source ----------------------------------------------------
raw = open(SRC, encoding='utf-8').read().splitlines()
sites = pd.read_csv(io.StringIO('\n'.join(raw[1:])), sep='\t')

# ---- audited fixes -----------------------------------------------------------
sites.loc[sites['site'] == 'D6', 'dir_dec'] = D6_DEC
sites['location'] = 'central North Greenland'
sites['geologic_classes'] = 'Intrusive'
sites['geologic_types'] = 'Volcanic Dike'
sites['lithologies'] = 'Diabase'
sites['method_codes'] = 'LP-DIR-AF:LP-DIR-T:DE-FM'
sites['citations'] = CIT
sites['result_type'] = 'i'
sites['dir_tilt_correction'] = 0
sites['age'] = AGE
sites['age_low'] = AGE_LOW
sites['age_high'] = AGE_HIGH
sites = sites.drop(columns=[c for c in ['age_sigma'] if c in sites.columns])

# recompute signed per-site VGPs at the corrected locality
sites = sites.drop(columns=[c for c in ['vgp_lat', 'vgp_lon', 'vgp_dp', 'vgp_dm'] if c in sites.columns])


def vgp_row(r):
    plon, plat, dp, dm = pmag.dia_vgp(r['dir_dec'], r['dir_inc'], r['dir_alpha95'], r['lat'], r['lon'])
    return pd.Series({'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
                      'vgp_dp': round(dp, 1), 'vgp_dm': round(dm, 1)})


sites = pd.concat([sites, sites.apply(vgp_row, axis=1)], axis=1)

good = sites[sites['result_quality'] != 'b']
blk = pmag.flip(ipmag.make_di_block(good['vgp_lon'].tolist(), good['vgp_lat'].tolist()), combine=True)
p = pmag.fisher_mean(blk)
if p['inc'] < 0:
    p['dec'] = (p['dec'] + 180) % 360
    p['inc'] = -p['inc']
dblk = pmag.flip(ipmag.make_di_block(good['dir_dec'].tolist(), good['dir_inc'].tolist()), combine=True)
d = pmag.fisher_mean(dblk)
if d['inc'] < 0:
    d['dec'] = (d['dec'] + 180) % 360
    d['inc'] = -d['inc']

locs = pd.DataFrame([{
    'location': 'central North Greenland', 'location_type': 'Region',
    'result_name': 'Central North Greenland (Victoria Fjord) dolerite dykes ca. 1382 Ma pole',
    'result_type': 'a', 'result_quality': 'g', 'sites': ':'.join(good['site'].tolist()),
    'method_codes': 'LP-DIR-AF:LP-DIR-T:DE-FM:DE-VGP:ST-C', 'citations': CIT_LOC,
    'geologic_classes': 'Intrusive', 'lithologies': 'Diabase',
    'lat_s': 81.5, 'lat_n': 81.5, 'lon_w': 315.3, 'lon_e': 315.3,
    'age': AGE, 'age_low': AGE_LOW, 'age_high': AGE_HIGH, 'age_unit': 'Ma',
    'dir_tilt_correction': 0,
    'dir_dec': round(d['dec'], 1), 'dir_inc': round(d['inc'], 1),
    'dir_alpha95': round(d['alpha95'], 1), 'dir_k': round(d['k'], 1), 'dir_n_sites': int(d['n']),
    'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
    'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1), 'pole_n_sites': int(p['n']),
    'pole_reversed_perc': 0, 'contact_test': 'C+',
    'continent_ocean': 'Greenland', 'country': 'Greenland',
    'description': (
        "Paleomagnetic pole for a swarm of ca. 1382 Ma dolerite dykes cutting "
        "Archaean crystalline basement on nunataks at the head of Victoria Fjord, "
        "central North Greenland (Abrahamsen & Van der Voo, 1987). The pole is the "
        "Fisher mean of the in-situ site virtual geomagnetic poles from the nine "
        "accepted dykes (a tenth dyke, oriented by magnetic compass under overcast, "
        "gave an anomalous direction and is excluded). The single-polarity "
        "characteristic remanence, carried by Ti-poor titanomagnetite, was isolated "
        "by alternating-field and thermal demagnetization. A detailed positive "
        "baked-contact test supports a primary origin: baked gneiss and "
        "quartz-diorite within roughly one dyke-width of the contact carry the dyke "
        "direction, whereas unbaked gneiss at 60 m retains a pre-dyke direction. The "
        "swarm is petrographically correlated with, and antiparallel in polarity to, "
        "the Midsommerso Dolerites and Zig-Zag Dal Basalts of eastern North "
        "Greenland, and the pole is used in the reconstructed Laurentia frame after "
        "rotating Greenland to North America. The age is adopted by correlation with "
        "the U-Pb baddeleyite-dated eastern North Greenland magmatism (Upton et al., "
        "2005); the paper infers ca. 1250 Ma from Rb-Sr on the correlative "
        "intrusives.")}])


def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a', encoding='utf-8')


def validate(path):
    import sys
    root = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
    sys.path.insert(0, os.path.join(root, 'scripts'))
    from validate_magic_contribution import validate_upload_file
    return validate_upload_file(path, tables=['locations', 'sites'])


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Abrahamsen1987_Victoria_Fjord_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    f.write('tab delimited\tlocations\n'); locs.to_csv(f, sep='\t', index=False)
    f.write('>>>>>>>>>>\n')
    f.write('tab delimited\tsites\n'); sites.to_csv(f, sep='\t', index=False)
print(f'-I- Victoria: sites {len(sites)} ({len(good)} accepted), pole '
      f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} K {p["k"]:.1f} N {int(p["n"])}; '
      f'dir {d["dec"]:.1f}/{d["inc"]:.1f} a95 {d["alpha95"]:.1f} k {d["k"]:.1f}')
validate(combined)
