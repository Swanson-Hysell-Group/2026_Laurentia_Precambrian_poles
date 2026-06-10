"""Convert the Whitmeyer & Karlstrom (2007) basement-province shapefile to a
committed GeoJSON used by the interactive pole map.

Source: ``Whitmeyer2007.shp`` from the sibling repository
``Laurentia_Paleogeography`` (Whitmeyer, S.J. & Karlstrom, K.E., 2007, Tectonic
model for the Proterozoic growth of North America, Geosphere 3, 220-259,
https://doi.org/10.1130/GES00055.1). The shapefile is NAD83 geographic
(lon/lat), so coordinates are written straight to GeoJSON (Leaflet/WGS84-ready).

Polygons are grouped into the six age categories used on the static
compilation map, dissolved per category (unary_union) and simplified so the
result is small enough to embed in the map. The output
``data/geologic_provinces/Whitmeyer2007_provinces.geojson`` is committed so the
book build is self-contained (it does not need the sibling repo or geopandas).

Run once (or when the source changes), in an env with pyshp + shapely::

    python scripts/build_province_geojson.py
"""

import json
import os
from collections import defaultdict

import shapefile  # pyshp
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SHP = ("/Users/penokean/0000_GitHub/Laurentia_Paleogeography/Data/"
       "Whitmeyer2007_GIS/Whitmeyer2007.shp")
OUT = os.path.join(ROOT, "data", "geologic_provinces",
                   "Whitmeyer2007_provinces.geojson")

SIMPLIFY_TOL = 0.03  # degrees (~3 km); coarse boundaries are fine at this scale

# Age category and fill color for each Whitmeyer 'Unit' value, following the
# coloring of the static compilation map.
CATEGORY = {
    "> 2.5 Ga Archean crust": "Archean provinces",
    "1.9 - 1.8 Ga reworked Archean crust": "Reworked Archean",
    "2.0 - 1.8 Ga juvenile orogens": "Paleoproterozoic",
    "2.0 - 1.8 Ga juvenile arcs": "Paleoproterozoic",
    "1.80 - 1.76 Ga juvenile arcs": "Paleoproterozoic",
    "1.76 - 1.72 Ga juvenile crust": "Paleoproterozoic",
    "1.72 - 1.68 Ga granitoids": "Paleoproterozoic",
    "1.72 - 1.68 Ga juvenile arcs": "Paleoproterozoic",
    "1.69 - 1.65 Ga juvenile crust": "Paleoproterozoic",
    "1.65 - 1.60 Ga granitoids": "Paleoproterozoic",
    "approx 1.70 Ga quartzite deposits": "Paleoproterozoic",
    "approx. 1.65 Ga quartzite deposits": "Paleoproterozoic",
    "1.55 - 1.35 Ga juvenile crust": "Mesoproterozoic",
    "1.45 - 1.35 Ga granitoids": "Mesoproterozoic",
    "1.3 - 1.0 Ga collisional orogens": "Mesoproterozoic",
    "1.3 - 0.95 Ga granitoids": "Mesoproterozoic",
    "1.2 - 1.1 Ga Midcontinent rift system": "Mesoproterozoic",
    "Basin": "Belt-Purcell Supergroup",
    "2.5 - 2.0 Ga miogeoclinal sediments": "Belt-Purcell Supergroup",
    "Eastern rift basins": "Neoprot–Cambrian rift",
    "< 0.78 Ga Windermere Supergroup": "Neoprot–Cambrian rift",
}

# Order controls draw order (Archean drawn under, younger over) and legend order.
CATEGORY_COLOR = {
    "Archean provinces": "#d9d9d9",
    "Reworked Archean": "#969696",
    "Paleoproterozoic": "#f4a6c0",
    "Mesoproterozoic": "#fda07a",
    "Belt-Purcell Supergroup": "#fbf3b0",
    "Neoprot–Cambrian rift": "#fed8b1",
}


def main():
    reader = shapefile.Reader(SHP)
    geoms = defaultdict(list)
    skipped = set()
    for rec in reader.iterShapeRecords():
        cat = CATEGORY.get(rec.record["Unit"])
        if cat is None:
            skipped.add(rec.record["Unit"])
            continue
        geom = shape(rec.shape.__geo_interface__)
        if not geom.is_valid:
            geom = geom.buffer(0)
        geoms[cat].append(geom)

    features = []
    for cat in CATEGORY_COLOR:
        if cat not in geoms:
            continue
        merged = unary_union(geoms[cat]).simplify(
            SIMPLIFY_TOL, preserve_topology=True)
        features.append({
            "type": "Feature",
            "properties": {"category": cat, "color": CATEGORY_COLOR[cat]},
            "geometry": mapping(merged),
        })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump({"type": "FeatureCollection", "features": features}, fh)

    size_kb = os.path.getsize(OUT) / 1024
    print(f"Wrote {os.path.relpath(OUT, ROOT)} "
          f"({len(features)} categories, {size_kb:.0f} KB)")
    if skipped:
        print(f"Unmapped Unit values skipped: {sorted(skipped)}")


if __name__ == "__main__":
    main()
