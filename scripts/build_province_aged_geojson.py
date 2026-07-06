"""Build an age-resolved Whitmeyer province GeoJSON for the reconstruction figure.

The committed ``Whitmeyer2007_provinces.geojson`` (built by
``build_province_geojson.py``) dissolves the 21 Whitmeyer & Karlstrom (2007)
``Unit`` values into six coarse color categories -- fine for the static map, but
it discards the per-unit ages, so every Mesoproterozoic (and Paleoproterozoic)
domain would appear at once in the through-time reconstruction figure.

This script instead emits one feature per unit carrying an ``appearance_age`` --
the reconstruction age (Ma) at or before which that domain is drawn -- while
keeping the six-category ``category``/``color`` so the figure looks the same.
Two ideas set the appearance ages:

- *Accretionary* units (new crust added to Laurentia) appear at the older
  (first-formation) bound of their age range, so the continent grows by
  accretion through the sequence.
- *Reworking / intraplate* units (the 1.2-1.1 Ga Midcontinent Rift and the
  1.3-0.95 Ga granitoids) were emplaced into lithosphere that had already
  accreted. Making them appear at their magmatic age would punch holes into the
  older continent and imply fake growth, so instead each such polygon inherits
  the appearance age of the nearest accretionary unit (its host crust). The
  Midcontinent Rift accordingly appears piecewise with its Archean,
  Paleoproterozoic, and early-Mesoproterozoic hosts.

Reworking units are kept as individual polygons (not dissolved) so the host age
can vary along their length; accretionary units are dissolved per unit.

Output: ``data/geologic_provinces/Whitmeyer2007_provinces_aged.geojson``
(committed so the book build stays self-contained). Run in an env with
pyshp + shapely::

    python scripts/build_province_aged_geojson.py
"""

import json
import os
from collections import defaultdict

import shapefile  # pyshp
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

from build_province_geojson import CATEGORY, CATEGORY_COLOR, SHP, SIMPLIFY_TOL

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "geologic_provinces",
                   "Whitmeyer2007_provinces_aged.geojson")

# Appearance age (Ma) for each accretionary unit: the older (first-formation)
# bound of its age range, i.e. the age by which that crust had accreted to
# Laurentia. Reworking units (below) are not listed here -- they inherit a host
# age. Grenville collisional crust uses its older 1.3 Ga bound (a user choice);
# raise it toward ~1350 to make it enter one panel earlier.
APPEARANCE_AGE = {
    "> 2.5 Ga Archean crust": 2500,
    "1.9 - 1.8 Ga reworked Archean crust": 1900,
    "2.0 - 1.8 Ga juvenile orogens": 2000,
    "2.0 - 1.8 Ga juvenile arcs": 2000,
    "1.80 - 1.76 Ga juvenile arcs": 1800,
    "1.76 - 1.72 Ga juvenile crust": 1760,
    "1.72 - 1.68 Ga granitoids": 1720,
    "1.72 - 1.68 Ga juvenile arcs": 1720,
    "1.69 - 1.65 Ga juvenile crust": 1690,
    "1.65 - 1.60 Ga granitoids": 1650,
    "approx 1.70 Ga quartzite deposits": 1700,
    "approx. 1.65 Ga quartzite deposits": 1650,
    "1.55 - 1.35 Ga juvenile crust": 1550,
    "1.45 - 1.35 Ga granitoids": 1450,
    "1.3 - 1.0 Ga collisional orogens": 1300,
    "2.5 - 2.0 Ga miogeoclinal sediments": 1470,  # Belt-Purcell basin fill
    "Basin": 1470,
    "Eastern rift basins": 780,
    "< 0.78 Ga Windermere Supergroup": 780,
}

# Units emplaced into already-accreted lithosphere: appearance age is inherited
# from the nearest accretionary unit (host crust), and until their own magmatic
# age (older bound of their range) they are rendered in the host's color -- so the
# Midcontinent Rift, say, reads as the Archean/Paleoproterozoic crust it cuts
# across until ca. 1.2 Ga, then as its own rift color.
REWORKING_MAGMATIC_AGE = {
    "1.2 - 1.1 Ga Midcontinent rift system": 1200,
    "1.3 - 0.95 Ga granitoids": 1300,
}


def feature(geom, unit, appearance_age, extra=None):
    cat = CATEGORY[unit]
    props = {
        "unit": unit,
        "category": cat,
        "color": CATEGORY_COLOR[cat],
        "appearance_age": appearance_age,
    }
    if extra:
        props.update(extra)
    return {
        "type": "Feature",
        "properties": props,
        "geometry": mapping(geom.simplify(SIMPLIFY_TOL, preserve_topology=True)),
    }


def main():
    reader = shapefile.Reader(SHP)
    by_unit = defaultdict(list)
    skipped = set()
    for rec in reader.iterShapeRecords():
        unit = rec.record["Unit"]
        if unit not in CATEGORY:
            skipped.add(unit)
            continue
        geom = shape(rec.shape.__geo_interface__)
        if not geom.is_valid:
            geom = geom.buffer(0)
        by_unit[unit].append(geom)

    # Dissolved accretionary geometries, used both as output features and as the
    # host candidates whose appearance age the reworking polygons inherit.
    accretion = {u: unary_union(g) for u, g in by_unit.items()
                 if u not in REWORKING_MAGMATIC_AGE}

    features = []
    for unit, geom in accretion.items():
        features.append(feature(geom, unit, APPEARANCE_AGE[unit]))

    for unit, magmatic_age in REWORKING_MAGMATIC_AGE.items():
        for poly in by_unit.get(unit, []):
            probe = poly.representative_point()
            host_unit = min(accretion,
                            key=lambda u: accretion[u].distance(probe))
            features.append(feature(poly, unit, APPEARANCE_AGE[host_unit], {
                "reworking": True,
                "magmatic_age": magmatic_age,
                "host_category": CATEGORY[host_unit],
            }))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump({"type": "FeatureCollection", "features": features}, fh)

    size_kb = os.path.getsize(OUT) / 1024
    print(f"Wrote {os.path.relpath(OUT, ROOT)} "
          f"({len(features)} features, {size_kb:.0f} KB)")
    ages = sorted({f["properties"]["appearance_age"] for f in features},
                  reverse=True)
    print(f"appearance ages present: {ages}")
    if skipped:
        print(f"Unmapped Unit values skipped: {sorted(skipped)}")


if __name__ == "__main__":
    main()
