"""Write Denyszyn2009_updated.csv: the unit VGP list of the Denyszyn et al. (2009) Franklin
grand means (D. A. D. Evans compilation, ../Denyszyn2009.csv) with updated values beside
each original column and a comment on what changed.

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026). The updated values come from the site table of the MagIC
contribution (../sites.txt, written by build_franklin_grand_mean.py), whose site values
were transcribed from the page images of each paper's site table.

Columns
- site ... common_inc: as in Denyszyn2009.csv.
- updated_vgp_lat_rotated, updated_vgp_lon_rotated: unit VGP in the polarity of the
  grand mean (Pacific pole), restored to North America for the Denyszyn et al. (2009)
  sites: Devon and Ellesmere by -20 deg about 72 N, 274 E (the Ellesmere microplate
  rotation of Denyszyn et al., 2009) and Greenland by -13.8 deg about 67.5 N, 241.5 E
  (Roest & Srivastava, 1989, as in the compilation workbook).
- updated_vgp_lat, updated_vgp_lon: the same VGP in present-day coordinates (differs from
  the rotated VGP only for the Denyszyn et al. sites).
- updated_polarity: polarity of the observed direction (westerly declination = n).
- magic_site: site name in ../sites.txt.
- change_deg: angle between the original and updated rotated VGPs.
- comment: what changed and why.
"""
import os, re
import numpy as np
import pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ev = pd.read_csv(os.path.join(OUT, 'Denyszyn2009.csv'), encoding='utf-8-sig')
sites = pd.read_csv(os.path.join(OUT, 'sites.txt'), sep='\t', skiprows=1).set_index('site')

PFX = {'Fahrig et al. 1971': 'F71-', 'Fahrig+Schwarz 1973': 'FS73-', 'Park 1974': 'P74-',
       'Palmer+Hayatsu 1975': 'PH75-', 'Park 1981': 'P81-', 'Christie+Fahrig 1983': 'CF83-',
       'Palmer et al. 1983': 'P83-', 'Denyszyn et al. 2009': ''}
EULER = {'Devon': [72.0, 274.0, -20.0], 'Greenland': [67.5, 241.5, -13.8]}
GRAND_MEAN = (162.0, 6.7)   # (lon, lat), for the polarity of the listed VGPs


def magic_site(r):
    prefix = [v for k, v in PFX.items() if r.reference.startswith(k)][0]
    site = str(r.site).strip()
    if prefix == 'CF83-' and '/' in site:
        return 'CF83-' + site.replace('/', '-')
    if prefix == '':
        return re.split(r'[ :*(]', site)[0]
    return prefix + site.split()[0]


def to_pacific(lat, lon):
    if pmag.angle([lon, lat], list(GRAND_MEAN))[0] > 90:
        return -lat, (lon + 180) % 360
    return lat, lon % 360


COMMENTS = {
    'GF': ("Denyszyn et al. (2009) Table 1 VGP (14.4 N, 177.0 E) used; the compilation value "
           "is calculated from the direction (14.3 N, 177.0 E)."),
    'F71-29': ("Relocated: Table 2 prints 54.5 N, 117.0 W (and the VGP 13 N, 153 E follows "
               "from it), but the site is the NE-trending dike near Hardisty Lake (text p. "
               "457), and Table 1 ties it to K-Ar determination GSC 64-64, whose sample "
               "FA-162-62 is at 64 deg 27 min N, 117 deg 00 min W (Wanless et al., 1966, GSC "
               "Paper 65-17, p. 46). VGP calculated from the direction (278/20) at 64.45 N."),
}
rows = []
for r in ev.itertuples():
    key = magic_site(r)
    s = sites.loc[key]
    lat, lon = to_pacific(s.vgp_lat, s.vgp_lon)
    rlat, rlon = lat, lon
    area = None
    if r.reference.startswith('Denyszyn'):
        area = 'Greenland' if 'Greenland' in r.reference else 'Devon'
        la, lo = pmag.pt_rot(EULER[area], [s.vgp_lat], [s.vgp_lon])
        rlat, rlon = to_pacific(la[0], lo[0])
    change = pmag.angle([r.vgp_lon_rotated, r.vgp_lat_rotated], [rlon, rlat])[0]
    pol = 'n' if np.cos(np.radians(s.dir_dec - 285)) > 0 else 'r'
    if key in COMMENTS:
        comment = COMMENTS[key]
    elif area == 'Devon':
        comment = ('Restored with the -20 deg microplate rotation of Denyszyn et al. (2009) '
                   'about 72 N, 274 E; the compilation used -35 deg. With -20 deg the 78 '
                   'units give the published D+09cg pole (8.4 N, 163.8 E, A95 2.8).')
    elif change < 0.1:
        comment = 'No change (agrees within rounding).'
    elif key.startswith('F71-'):
        comment = ('VGP calculated from the direction at the coordinates of Fahrig et al. '
                   '(1971) Table 2; the compilation uses their printed VGP, in whole degrees.')
    else:
        comment = 'Changed; see the MagIC site description.'
    if pol != r.polarity:
        comment += f' Polarity {r.polarity} -> {pol}.'
    rows.append({
        'site': r.site,
        'vgp_lat_rotated': r.vgp_lat_rotated, 'updated_vgp_lat_rotated': round(rlat, 2),
        'vgp_lon_rotated': r.vgp_lon_rotated, 'updated_vgp_lon_rotated': round(rlon, 2),
        'polarity': r.polarity, 'updated_polarity': pol,
        'vgp_lat': r.vgp_lat, 'updated_vgp_lat': round(lat, 2),
        'vgp_lon': r.vgp_lon, 'updated_vgp_lon': round(lon, 2),
        'reference': r.reference, 'note': r.note,
        'common_dec': r.common_dec, 'common_inc': r.common_inc,
        'magic_site': key, 'change_deg': round(change, 2), 'comment': comment})
out = pd.DataFrame(rows)
out.to_csv(os.path.join(OUT, 'Denyszyn2009_updated.csv'), index=False)


def mean(lat, lon):
    p = pmag.fisher_mean([[o, a, 1.] for a, o in zip(lat, lon)])
    return f"{p['inc']:.2f} N / {p['dec']:.2f} E, A95 {p['alpha95']:.2f}, K {p['k']:.1f}, N {p['n']}"


prev = ~out['reference'].str.startswith('Denyszyn')
for label, sel in [('D+09can (56 units)', prev), ('D+09cg (78 units)', slice(None))]:
    o, u = out[sel], out[sel]
    print(f'{label:20s} original {mean(o.vgp_lat_rotated, o.vgp_lon_rotated)}')
    print(f'{"":20s} updated  {mean(u.updated_vgp_lat_rotated, u.updated_vgp_lon_rotated)}')
print(out.loc[out.change_deg >= 0.1, ['site', 'change_deg', 'comment']].to_string(max_colwidth=70))
