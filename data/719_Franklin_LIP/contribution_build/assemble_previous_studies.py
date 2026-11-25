"""Assemble the previous-study Franklin site table (the inputs of pole D+09can of
Denyszyn et al., 2009) from transcriptions of each paper's site table read from page
images (source_transcriptions/paper_*.csv) and from the site list of the Denyszyn et al.
(2009) compilation (D. A. D. Evans, compilation workbook; source_transcriptions/evans_*.csv),
which sets which sites enter the grand mean and gives coordinates where the papers give
none. Prepared for Swanson-Hysell et al. (2026); writes Franklin_previous_studies_sites.csv,
read by build_franklin_grand_mean.py."""
import os, re
import numpy as np, pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(HERE, 'source_transcriptions')
OUT = os.path.join(HERE, 'Franklin_previous_studies_sites.csv')
P = lambda f: pd.read_csv(f'{S}/paper_{f}.csv', dtype={'site': str})
E = lambda f: pd.read_csv(f'{S}/evans_{f}.csv', dtype={'site': str})
rows = []


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def add(**kw):
    rows.append(kw)


def ncoll(note, pat):
    m = re.search(pat, str(note))
    return int(m.group(1)) if m else None


# ---------------- Fahrig et al. 1971 ---------------------------------------------
p, e = P('Fahrig1971'), E('Fahrig1971')
e['key'] = e['site'].str.split().str[0]
for r in p.itertuples():
    if str(r.polarity).startswith('anomalous'):
        continue                    # Group B: excluded by the authors
    w = e[e.key == r.site]
    status = 'accepted' if (len(w) and w.weight.iloc[0] == 1.0) else 'rejected'
    grp = r.group
    gtype = 'Sill:Volcanic Dike' if grp.startswith('Coronation') else 'Volcanic Dike'
    lat, coord_source = r.lat, 'Fahrig et al. (1971) Table 2'
    vgp_paper, note = (r.vgp_lat, r.vgp_lon), ''
    if r.site == '29':
        # Table 2 prints 54.5 N, 117.0 W. The site is the NE-trending dike near Hardisty
        # Lake (text, p. 457), and Table 1 ties it to K-Ar determination GSC 64-64, whose
        # sample FA-162-62(5) is located at 64 deg 27' N, 117 deg 00' W (Wanless et al.,
        # 1966, GSC Paper 65-17, p. 46). The tabulated VGP follows from 54.5 N, so the VGP
        # is calculated from the direction at the corrected latitude.
        lat, vgp_paper = 64.45, (np.nan, np.nan)
        coord_source = ('64 deg 27 min N, 117 deg 00 min W, the location of K-Ar sample '
                        'FA-162-62 (GSC 64-64; Wanless et al., 1966, Geological Survey of '
                        'Canada Paper 65-17) from this site')
        note = ('NE-trending dike near Hardisty Lake, south of Great Bear Lake. Table 2 of '
                'Fahrig et al. (1971) prints the latitude as 54.5 N and a VGP of 13 N, 153 E '
                'that follows from it; the VGP here is calculated from the direction at the '
                'latitude of the K-Ar sample.')
    add(study='Fahrig1971', citation='10.1139/e71-047', site_printed=r.site,
        site=f'F71-{r.site}', group=grp, lat=lat, lon=r.lon,
        coord_source=coord_source, geologic_types=gtype,
        lithologies='Diabase', n_samples=int(r.n_samples),
        n_total_samples=ncoll(r.note, r'NRM: N=(\d+)'), n_specimens=None,
        treatment=r.treatment, dec=r.dec, inc=r.inc, k=r.k, a95=r.a95,
        vgp_lat_paper=vgp_paper[0], vgp_lon_paper=vgp_paper[1], polarity=r.polarity,
        d09=status, method='AF blanket', note=note)

# ---------------- Fahrig & Schwarz 1973 ------------------------------------------
p, e = P('FahrigSchwarz1973'), E('FahrigSchwarz1973')
e['key'] = e['site'].str.split().str[0]
for r in p.itertuples():
    if r.site == '7':
        continue                    # rejected by the authors; no direction
    w = e[e.key == r.site].iloc[0]
    add(study='FahrigSchwarz1973', citation='10.1139/e73-057', site_printed=r.site,
        site=f'FS73-{r.site}', group='southern and central Baffin Island dikes',
        lat=w.lat, lon=w.lon,
        coord_source=('not given in Fahrig & Schwarz (1973); coordinates as used in the '
                      'compilation of Denyszyn et al. (2009)'),
        geologic_types='Volcanic Dike', lithologies='Diabase',
        n_samples=int(r.n_samples), n_total_samples=ncoll(r.note, r'N collected=(\d+)'),
        n_specimens=None, treatment=r.treatment, dec=r.dec, inc=r.inc, k=r.k, a95=np.nan,
        vgp_lat_paper=np.nan, vgp_lon_paper=np.nan, polarity=r.polarity,
        d09='accepted' if w.weight == 1.0 else 'rejected', method='AF blanket')

# ---------------- Park 1974 (Franklin sites; 'Mean' = dyke + contact rock) --------
p, e = P('Park1974'), E('Park1974')
# samples collected: dike + contact samples used, plus the discarded samples in brackets
COLLECTED = {'31': 15, '33': 15, '16': 8}
for nm, key, status in [('31 Mean', '31', 'accepted'), ('33 Mean', '33', 'accepted'),
                        ('16 Mean', '16', 'rejected')]:
    r = p[p.site == nm].iloc[0]
    reason = ('Not used in the grand mean of Denyszyn et al. (2009): this dike in southern '
              'Ontario lies outside the region interpreted by Denyszyn et al. (2009) to have '
              'had Franklin magmatism.') if key == '16' else ''
    add(study='Park1974', citation='10.1139/e74-098', site_printed=nm,
        site=f'P74-{key}', group='Franklin dikes (Park 1974)', lat=r.lat, lon=r.lon,
        coord_source='Park (1974) Table 1', geologic_types='Volcanic Dike:Baked Contact',
        lithologies='Diabase', n_samples=int(r.n_samples), n_total_samples=COLLECTED[key],
        n_specimens=None, treatment='AF blanket (dike and contact samples meaned)',
        dec=r.dec, inc=r.inc, k=r.k_calc, a95=r.a95_calc,
        vgp_lat_paper=r.vgp_lat, vgp_lon_paper=r.vgp_lon, polarity=r.polarity[:1].upper(),
        d09=status, method='AF blanket', reject_reason=reason)

# ---------------- Palmer & Hayatsu 1975 -------------------------------------------
p, e = P('PalmerHayatsu1975'), E('PalmerHayatsu1975')
for r in p.itertuples():
    if r.site in ('1L', '9L', 'Mean'):
        continue
    w = e[e.site.str.strip() == r.site]
    m = re.search(r'ca\. ([\d.]+)N, ([\d.]+)E', r.coord_note)
    lava = r.site.endswith('L')
    add(study='PalmerHayatsu1975', citation='10.1139/e75-130', site_printed=r.site,
        site=f'PH75-{r.site}', group=r.group, lat=float(m.group(1)), lon=float(m.group(2)),
        coord_source=('not given in Palmer & Hayatsu (1975) (map only); area coordinates '
                      'back-calculated from the tabulated VGPs'),
        geologic_types='Lava Flow' if lava else 'Sill',
        lithologies='Basalt' if lava else 'Diabase', n_samples=int(r.n_samples),
        n_total_samples=None, n_specimens=None, treatment=r.treatment, dec=r.dec,
        inc=r.inc, k=r.k, a95=r.a95, vgp_lat_paper=r.vgp_lat, vgp_lon_paper=r.vgp_lon,
        polarity='', d09='accepted' if (len(w) and w.weight.iloc[0] == 1.0) else 'rejected',
        method='AF blanket')

# ---------------- Park 1981 -------------------------------------------------------
p = P('Park1981')
for r in p.itertuples():
    if r.site not in ('1', '3', '6', '10', '14'):
        continue                    # site 5 (2 specimens) and summary rows
    # Table 1 alpha95 follows from the specimen count and R, except site 1 (printed 9;
    # N = 4, R = 3.91 give 16.1). The Denyszyn et al. (2009) compilation applies its filter
    # to the specimen count and the alpha95 recalculated from it.
    n_sp = r.n_specimens
    a95_spec = np.degrees(np.arccos(1 - (n_sp - r.R) / r.R * (20 ** (1 / (n_sp - 1)) - 1)))
    note = 'Directions are corrected for bedding (dips less than 6 deg).'
    if r.site == '1':
        note += (' The printed alpha95 (9 deg) differs from the 16.1 deg given by N = 4 '
                 'specimens and R = 3.91.')
    if r.site == '10':
        note += (' Enters the grand mean of Denyszyn et al. (2009) with its 4 specimens, '
                 'from 3 samples.')
    reason = ('Not used in the grand mean of Denyszyn et al. (2009), whose filter (N >= 4 and '
              'alpha95 <= 15 deg) is applied in their compilation to the specimen count and '
              'the alpha95 calculated from it'
              + ('.' if r.site == '1' else f' ({a95_spec:.1f} deg).')) if r.site in ('1', '3') else ''
    add(study='Park1981', citation='10.1139/e81-151', site_printed=r.site,
        site=f'P81-{r.site}', group='Brock Inlier sills', lat=68.37, lon=238.5,
        coord_source=('not given in Park (1981) (Fig. 1b only); common coordinate that '
                      'reproduces the tabulated VGPs'),
        geologic_types='Sill', lithologies='Gabbro:Diabase',
        n_samples=int(r.n_samples), n_total_samples=None, n_specimens=int(r.n_specimens),
        treatment='AF and thermal', dec=r.dec, inc=r.inc, k=r.k, a95=r.a95,
        vgp_lat_paper=r.vgp_lat, vgp_lon_paper=r.vgp_lon, polarity='R',
        d09='accepted' if r.site in ('6', '10', '14') else 'rejected', method='AF+TH PCA',
        tilt_correction=100, note=note, reject_reason=reason)

# ---------------- Christie & Fahrig 1983 (Franklin dykes) -------------------------
p, e = P('ChristieFahrig1983'), E('ChristieFahrig1983')
e['key'] = e['site'].str.split().str[0]
fr = p[p.group == 'Franklin dykes']
UNITS = {'03': '03/7922/7923/7924', '7922': '03/7922/7923/7924', '7923': '03/7922/7923/7924',
         '7924': '03/7922/7923/7924', '7917': '7917/7920', '7920': '7917/7920'}
for r in fr.itertuples():
    s = r.site
    s2 = s.lstrip('0') if s.startswith('0') else s
    unit = UNITS.get(s)
    w = e[e.key == (unit or s)].iloc[0]
    add(study='ChristieFahrig1983', citation='10.1139/e83-025', site_printed=s,
        site=f'CF83-{s}', group='Franklin dykes, northwest Baffin Island',
        lat=w.lat, lon=w.lon,
        coord_source=('not given in Christie & Fahrig (1983) (Fig. 1 map only); '
                      'coordinates as used in the compilation of Denyszyn et al. '
                      '(2009), consistent with Fig. 1 to ca. 0.1 deg'),
        geologic_types='Volcanic Dike', lithologies='Diabase',
        n_samples=int(r.n_samples), n_total_samples=ncoll(r.note, r'N cores collected=(\d+)'),
        n_specimens=None, treatment=r.treatment, dec=r.dec, inc=r.inc, k=r.k, a95=r.a95,
        vgp_lat_paper=np.nan, vgp_lon_paper=np.nan, polarity=r.polarity,
        d09='constituent' if unit else 'accepted', unit=unit or '', method='AF blanket')

# ---------------- Palmer et al. 1983 (upper-reverse + middle-normal lavas) --------
p, e = P('Palmer1983'), E('Palmer1983')
e['key'] = e['site'].str.split().str[0]
lav = p[p.group.astype(str).str.contains('upper reverse|middle normal')]
for r in lav.itertuples():
    s = r.site.split()[1]
    w = e[e.key == s].iloc[0]
    add(study='Palmer1983', citation='10.1139/e83-131', site_printed=r.site,
        site=f'P83-{s}', group=r.group.split(' (')[0] + ' zone, Natkusiak Formation',
        lat=w.lat, lon=w.lon,
        coord_source=('not given in Palmer et al. (1983) (Fig. 2 sections only); '
                      'coordinates as used in the compilation of Denyszyn et al. (2009)'),
        geologic_types='Lava Flow', lithologies='Basalt',
        n_samples=int(r.N_m2), n_total_samples=int(r.n_samples), n_specimens=None,
        treatment='AF and thermal; anchored PCA lines (method 2)',
        dec=r.dec_m2_anchored, inc=r.inc_m2_anchored, k=num(r.k_m2), a95=np.nan,
        vgp_lat_paper=np.nan, vgp_lon_paper=np.nan, polarity=r.polarity,
        d09='accepted' if w.weight == 1.0 else 'rejected',
        note=(f'Method (1) site mean {r.dec:.0f}/{r.inc:.0f} (N {r.n_samples}, k {r.k}); '
              f'mean of methods (1) and (2) {r.dec_mean12:.0f}/{r.inc_mean12:.0f}.'),
        method='anchored PCA')

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(df.groupby(['study', 'd09']).size())
print(len(df))
