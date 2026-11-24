"""Map the site coordinates of the Long Range dykes contribution (Murthy et al., 1992).

Murthy et al. (1992) give no site coordinates, so the coordinates in ``../sites.txt`` were
read from the maps of the paper (see the docstring of ``build_long_range.py``):

- PR84-191, PR84-199 and PR84-206 are at the U-Pb sample locality of Kamo & Gower (1994;
  53 deg 46.2' N, 56 deg 37.1' W), which they state corresponds to these sites.
- The other Sandwich Bay sites are read from the Fig. 2 sketch map, whose frame is labelled
  54 deg 00' N, 57 deg 35' W (upper left) and 53 deg 22' N, 56 deg 14' W (lower right). The
  sketch map is close to a plate carree projection of that frame: its 10 km scale bar
  agrees with the frame dimensions to ca. 2%.
- Dyke 6 (PR86-008, PR86-015) lies outside Fig. 2 and is read from Fig. 1, using the outline
  of the Fig. 2 frame drawn on Fig. 1 as the reference (ca. 5 km).

The figure shows the sites by dyke, the Fig. 2 frame, the Kamo & Gower (1994) locality and
the single reference locality (53.7 N, 303.3 E) used for the pole in earlier compilations.

Usage: python plot_site_coordinates.py
Writes Long_Range_site_coordinates.png next to this script.

Optional: python plot_site_coordinates.py --paper <path to the Murthy et al. (1992) PDF>
also writes Long_Range_site_coordinates_on_Fig2.png, the coordinates plotted on page 1225
of the paper (Fig. 2, georeferenced from its labelled frame corners) as a check of the map
reading. That image reproduces the published figure, so it is not part of the repository.
"""
import argparse
import os
import subprocess
import tempfile

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SITES = os.path.join(HERE, '..', 'sites.txt')

# Fig. 2 frame (deg N, deg E)
FIG2_N, FIG2_S = 54.0, 53 + 22 / 60
FIG2_W, FIG2_E = 360 - (57 + 35 / 60), 360 - (56 + 14 / 60)
KG94 = (53 + 46.2 / 60, 360 - (56 + 37.1 / 60))
REF = (53.7, 303.3)
PRIMARY = {'2', '4', '6'}
COLORS = {'1': '#9e9e9e', '2': '#1f77b4', '3': '#bdbdbd', '4': '#d62728', '5': '#757575',
          '6': '#2ca02c'}


def load_sites():
    s = pd.read_csv(SITES, sep='\t', skiprows=1)
    s = s[s['result_type'] == 'i'].copy()
    s['dyke'] = s['description'].str.extract(r'Sampling site in dyke (\d)')[0]
    s['source'] = 'Fig. 2'
    s.loc[s['description'].str.contains('Fig. 1 map'), 'source'] = 'Fig. 1'
    s.loc[s['description'].str.contains('Kamo & Gower'), 'source'] = 'Kamo & Gower (1994)'
    return s


# label offsets in points (x, y) for sites close together in the zoom panel
OFFSETS = {'PR84-171': (-8, 4), 'PR84-154': (-8, -4), 'PR84-191': (8, 10),
           'PR84-178': (8, -8), 'PR84-162': (8, 0), 'PR85-565': (8, -6)}


def label_groups(ax, s, extent, offsets=None, default=(6, 4)):
    """Label each distinct site position (names of co-located sites stacked)."""
    w, e, so, n = extent
    for (lat, lon), sub in s.groupby([s['lat'].round(3), s['lon'].round(3)], sort=False):
        if not (w < lon < e and so < lat < n):
            continue
        names = list(sub['site'])
        dx, dy = (offsets or {}).get(names[0], default)
        ax.annotate('\n'.join(names), xy=(lon, lat), xycoords=ccrs.PlateCarree(),
                    xytext=(dx, dy), textcoords='offset points', fontsize=6.5,
                    ha='right' if dx < 0 else 'left', va='center', zorder=6)


def draw(ax, s, extent, frame_label=True):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m'),
                   facecolor='#f2efe9', edgecolor='none')
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '10m'),
                   facecolor='none', edgecolor='0.4', linewidth=0.5)
    ax.plot([FIG2_W, FIG2_E, FIG2_E, FIG2_W, FIG2_W], [FIG2_N, FIG2_N, FIG2_S, FIG2_S, FIG2_N],
            'k--', linewidth=0.8, transform=ccrs.PlateCarree())
    if frame_label:
        ax.text(FIG2_W + 0.02, FIG2_N - 0.03, 'Fig. 2 frame', fontsize=7, va='top',
                transform=ccrs.PlateCarree())
    for (dyke, sub) in s.groupby('dyke'):
        prim = dyke in PRIMARY
        ax.scatter(sub['lon'], sub['lat'], s=42, marker='o' if prim else 's',
                   facecolor=COLORS[dyke], edgecolor='k', linewidth=0.5, zorder=4,
                   transform=ccrs.PlateCarree(),
                   label=f"Dyke {dyke} ({'primary' if prim else 'anomalous'})")
    ax.scatter([KG94[1]], [KG94[0]], marker='*', s=170, facecolor='none', edgecolor='k',
               linewidth=0.8, zorder=5, transform=ccrs.PlateCarree(),
               label='U-Pb locality, Kamo & Gower (1994)')
    ax.scatter([REF[1]], [REF[0]], marker='x', s=50, color='k', zorder=5,
               transform=ccrs.PlateCarree(), label='Single reference locality 53.7 N, 303.3 E')
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='0.6')
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'size': 7}


def main_figure(s, out):
    fig = plt.figure(figsize=(13, 4.4))
    a = fig.add_subplot(1, 3, 1, projection=ccrs.PlateCarree())
    draw(a, s, [301.9, 304.7, 52.9, 54.25])
    a.text(304.16, 53.30, 'Dyke 6\n(from Fig. 1)', fontsize=7, va='center',
           transform=ccrs.PlateCarree())
    a.set_title('(a) Sandwich Bay area, southeast Labrador', fontsize=9)
    b = fig.add_subplot(1, 3, 2, projection=ccrs.PlateCarree())
    zoom = [303.24, 303.50, 53.64, 53.82]
    draw(b, s, [FIG2_W, FIG2_E, FIG2_S, FIG2_N])
    outside_zoom = s[~((s['lon'] > zoom[0]) & (s['lon'] < zoom[1])
                       & (s['lat'] > zoom[2]) & (s['lat'] < zoom[3]))]
    label_groups(b, outside_zoom, [FIG2_W, FIG2_E, FIG2_S, FIG2_N])
    b.plot([zoom[0], zoom[1], zoom[1], zoom[0], zoom[0]],
           [zoom[2], zoom[2], zoom[3], zoom[3], zoom[2]], color='0.3', linewidth=0.8,
           transform=ccrs.PlateCarree())
    b.set_title('(b) Fig. 2 frame of Murthy et al. (1992)', fontsize=9)
    c = fig.add_subplot(1, 3, 3, projection=ccrs.PlateCarree())
    draw(c, s, zoom, frame_label=False)
    label_groups(c, s, zoom, OFFSETS)
    c.set_title('(c) Dykes 3, 4 and 5 (box in b)', fontsize=9)
    h, l = a.get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=4, fontsize=7.5, frameon=False)
    fig.suptitle('Long Range dykes (Murthy et al., 1992): site coordinates read from the '
                 'maps of the paper', fontsize=10)
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.16, top=0.88, wspace=0.22)
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)


def paper_overlay(s, pdf, out):
    """Plot the site coordinates on Fig. 2 of the paper (page 2 of the PDF, 300 dpi)."""
    from PIL import Image
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(['pdftoppm', '-r', '300', '-f', '2', '-l', '2', '-png', pdf,
                        os.path.join(tmp, 'p')], check=True)
        page = Image.open(os.path.join(tmp, sorted(os.listdir(tmp))[0])).convert('L')
        page.load()
    # pixel positions of the Fig. 2 frame lines on the 300 dpi page image
    x0, x1, y0, y1 = 1216, 2254, 1658, 2492
    crop = (1150, 1600, 2350, 2600)
    im = page.crop(crop)

    def px(lat, lon):
        return (x0 + (lon - FIG2_W) / (FIG2_E - FIG2_W) * (x1 - x0) - crop[0],
                y0 + (FIG2_N - lat) / (FIG2_N - FIG2_S) * (y1 - y0) - crop[1])

    fig, ax = plt.subplots(figsize=(10, 8.4))
    ax.imshow(im, cmap='gray')
    inside = s[(s['lat'] >= FIG2_S) & (s['lat'] <= FIG2_N)]
    for key, sub in inside.groupby([inside['lat'].round(3), inside['lon'].round(3)]):
        x, y = px(*key)
        ax.plot(x, y, 'o', mfc='none', mec='red', ms=11, mew=1.8)
        ax.text(x + 10, y - 10, '/'.join(sub['site'].str[-3:]), color='red', fontsize=7,
                weight='bold')
    x, y = px(*REF)
    ax.plot(x, y, 'x', color='blue', ms=9, mew=1.8)
    ax.axis('off')
    ax.set_title('Site coordinates of the contribution (red) on Fig. 2 of Murthy et al. '
                 '(1992); blue x: 53.7 N, 303.3 E', fontsize=9)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--paper', help='path to the Murthy et al. (1992) PDF (optional check)')
    args = p.parse_args()
    sites = load_sites()
    main_figure(sites, os.path.join(HERE, 'Long_Range_site_coordinates.png'))
    if args.paper:
        paper_overlay(sites, args.paper,
                      os.path.join(HERE, 'Long_Range_site_coordinates_on_Fig2.png'))
