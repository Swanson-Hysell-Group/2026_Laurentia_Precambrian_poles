# Baie des Moutons syenite — contribution build

How the site-level MagIC data for the two ca. 583 Ma Baie des Moutons (Mutton
Bay) syenite poles are assembled. The notebook
`pole_notebooks/583_Baie_des_Moutons.ipynb` loads the finished `../sites.txt` /
`../locations.txt`; everything here is provenance.

## Source

McCausland, P. J. A., Hankard, F., Van der Voo, R., & Hall, C. M. (2011).
Ediacaran paleogeography of Laurentia: Paleomagnetism and ⁴⁰Ar–³⁹Ar geochronology
of the 583 Ma Baie des Moutons syenite, Quebec. *Precambrian Research*, 187,
58–78. doi:[10.1016/j.precamres.2011.02.004](https://doi.org/10.1016/j.precamres.2011.02.004).

There is **no measurement-/specimen-level MagIC contribution** for this study, so
the contribution is built from the published **Table 1** site means (transcribed
into `BaieDesMoutons_build_contribution.py`).

## `BaieDesMoutons_build_contribution.py` → `../sites.txt`, `../locations.txt`

The Baie des Moutons complex carries two characteristic remanence components
that are nowhere observed to coexist and give two paleopoles:

- **ChRM A** — 8 steep, easterly syenite sites (PSD magnetite). Site means →
  Fisher pole **42.9°N, 331.9°E, A95 11.5°, N 8** (published 42.6°N/332.7°E,
  dp 11.7/dm 12.4). GPMDB 9364. High-paleolatitude interpretation.
- **ChRM B** — 6 shallow sites (feldspar-porphyry, aplite, carbonatite dykes +
  the late red-syenite site MB20, which is reversed and brought to common
  polarity). Fisher pole **−33.9°N, 321.9°E, A95 15.7°, N 6** (published
  −34.2°N/321.5°E, dp 10.9/dm 21.8). GPMDB 9365. Low-paleolatitude
  interpretation; resembles the younger ca. 565 Ma Sept-Îles pole.

Build details:

- **Coordinates** are geographic (`dir_tilt_correction = 0`): the nested-cone
  complex is interpreted as untilted since emplacement (varied dyke
  orientations), so no tilt correction is applied.
- **Site lat/lon** from Table 1, given as decimal minutes added to 50°N and
  58°W (NAD-27): `lat = 50 + min/60`, `lon (0–360 E) = 360 − (58 + min/60)`.
- `dir_n_samples` is the number of ChRM-bearing (endpoint + great-circle)
  specimens used in each site mean (Table 1 `e+g`). Site MB15 is great-circle
  only; its α95 is estimated as `140/√(k·n)`.
- VGP and dp/dm are computed from each site mean with `pmag.dia_vgp`; the pole
  is the polarity-unified Fisher mean of the site VGPs (reported in the
  hemisphere McCausland et al. chose for each component).
- Two **locations** are written, one per component/pole (the A and B site sets
  are disjoint, and each location carries its own `pole_*` result).
- **Age**: 583.4 ± 2.0 Ma, weighted-mean hornblende ⁴⁰Ar–³⁹Ar plateau age (sites
  MB44, MB20). Baked-contact tests at six sites were inconclusive.

`upload_magic` validates the result. This replaces the earlier stub notebooks
`583_Baie_A.ipynb` / `583_Baie_B.ipynb`, which carried only the compilation-level
values, with a single site-level notebook covering both poles.
