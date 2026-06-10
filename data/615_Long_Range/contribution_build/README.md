# Long Range dykes — contribution build

How the site-level MagIC data for the ca. 615 Ma Long Range dyke pole are
assembled. The notebook `pole_notebooks/615_Long_Range.ipynb` loads the finished
`../sites.txt` / `../locations.txt`; everything here is provenance.

## Source

Murthy, G. S., Gower, C., Tubrett, M., & Patzold, R. (1992). Paleomagnetism of
Eocambrian Long Range dykes and Double Mer Formation from Labrador, Canada.
*Canadian Journal of Earth Sciences*, 29(6), 1224–1234.
doi:[10.1139/e92-098](https://doi.org/10.1139/e92-098).

There is **no measurement-level MagIC contribution** for this study, so the
contribution is built from the published **Table 1** dyke means, digitized into
`Murthy1992_dyke_means.csv`.

## `LongRange_build_contribution.py` → `../sites.txt`, `../locations.txt`

- All six dykes are written. The three **primary** dykes (2, 4, 6) are
  `result_quality='g'`; the three **anomalous** dykes (1, 3, 5) are
  `result_quality='b'` (dropped by `pt.load_magic_sites` unless `drop_bad=False`).
  Murthy et al. interpret dykes 2, 4, 6 as carrying the primary ca. 615 Ma
  remanence (combined mean D = 124.8°, I = 55.5°, k = 48, α95 = 18.0°; paleopole
  10.8°S, 164.3°E).
- Geographic coordinates (`dir_tilt_correction = 0`); locality 53.7°N, 56.7°W
  (303.3°E); per-dyke VGPs + dp/dm from `pmag.dia_vgp`. `dir_n_samples` is the
  total samples per dyke (Table 1 `N` summed over its sampling sites).
- `locations.txt` carries the **primary-dyke (2, 4, 6)** pole
  (12.4°N, 345.4°E, A95 20.8°, N 3 — the northern antipode of the published
  10.8°S, 164.3°E), which reproduces the paper's primary-remanence direction
  exactly. The notebook computes its **pole of record** from the data directly
  (see below). GPMDB 6934–6936.

## Pole of record (in the notebook and the compilation)

This version of the compilation adopts **Murthy et al.'s three primary dykes
(2, 4, 6)** as the Long Range pole: **12.4°S, 165.4°E (A95 20.8°, N = 3)** —
reported in Murthy et al.'s published (southern) polarity, close to their
paleopole at 10.8°S, 164.3°E — reproducing the published primary-remanence
direction (124.8°/55.5°, k 48, α95 18.0°). `data/Laurentia_poles.csv` has been
**updated** to this value.

It replaces the prior Nordic Paleomagnetic Workshop value ("Dykes 1, 2, 3, 4 and
6", 19°N/355.3°E, A95 17.4°), which does not trace to Murthy's dyke means: it
follows the **Hodych et al. (2004)** recalculation — an eight-site mean
(D = 110.3°, I = 56.8°; "omitting only dyke #5 and combining sites from the same
dyke segment") converted to a pole at the reference site, with A95 taken as the
geometric mean of dp and dm. The NPW databases carried that mean forward
(consistent since the 2010 Luleå meeting) without redoing the VGP averaging, so
the listed value is not reproducible from the published site directions —
recomputing the five-dyke mean instead gives ~29°N/354°E with A95 ~31° (shown in
the notebook for comparison). The choice of the three primary dykes is reinforced
by **Kamo & Gower (1994)**, whose ca. 615 Ma U-Pb baddeleyite age on **dyke 4 (a
primary dyke)** directly dates the characteristic primary remanence.

## Correction applied

The earlier extraction listed **dyke 6 as 124.8°/55.5°** — that is actually the
paper's *combined* 3-dyke mean, not dyke 6. Dyke 6's own mean (its two sites
PR86-008 = 117.0°/66.0° and PR86-015 = 132.6°/68.6°) is **124.6°/67.5°**. With
the correct value, the Fisher mean of dykes 2, 4, 6 reproduces the published
124.8°/55.5° (k 48, α95 18.0°). The previous extraction also omitted the
anomalous **dyke 5** (279.1°/28.1°), now included.

## Change from the prior compilation

The pre-2022 Nordic compilation entry ("Dykes 1, 2, 3, 4 and 6", 19°N/355.3°E,
A95 17.4° — the Hodych et al. 2004 eight-site recalculation) has been **replaced**
in `data/Laurentia_poles.csv` by Murthy et al.'s three-primary-dyke pole
(12.4°N/345.4°E, B = 3, N = 69). The R-score is unchanged (5/7, Grade B).
