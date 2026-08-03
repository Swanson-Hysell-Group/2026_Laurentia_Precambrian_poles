# Callander Complex — contribution build

How the site-level MagIC data for the ca. 575 Ma Callander alkaline complex pole
are assembled. The notebook `pole_notebooks/575_Callander.ipynb` loads the
finished `../sites.txt` / `../locations.txt`; everything here is provenance.

## Source

Symons, D. T. A., & Chiasson, A. D. (1991). Paleomagnetism of the Callander
Complex and the Cambrian apparent polar wander path for North America.
*Canadian Journal of Earth Sciences*, 28(3), 355–363.
doi:[10.1139/e91-032](https://doi.org/10.1139/e91-032).

There is **no measurement-level MagIC contribution** for this study, so the
contribution is built from the published **Table 1** site means, digitized into
`Symons1991_site_means.csv` (columns: `site, rock_type, dec, inc, a_95, k, n`,
where `n` = endpoint + great-circle specimens used in each site mean).

## `Callander_build_contribution.py` → `../sites.txt`, `../locations.txt`

- **Component:** the characteristic **A** magnetization (magnetite + hematite),
  shown to be primary by a positive baked-contact test and a (weak) reversals
  test.
- **Coordinates:** geographic (`dir_tilt_correction = 0`) — the 3 km circular
  complex "has not been subsequently metamorphosed or tilted."
- **Locality:** 46.2°N, 79.4°W (280.6°E); per-site VGPs + dp/dm from
  `pmag.dia_vgp`.
- **Site selection:** sites with α95 ≤ 15° are accepted (`result_quality='g'`);
  site 3 (carbonatite dike, α95 = 25°) and site 29 (the contact-test host
  gneiss, α95 = 31°) are flagged `result_quality='b'` and dropped by
  `pt.load_magic_sites`. Site 22 (< 4 reliable specimens) was not digitized.
  The **26** accepted sites = 16 lamprophyre dikes + 5 syenitic-core sites
  (mesocratic/leucocratic nepheline syenite + potassium trachyte) + 5
  fenite-aureole sites; their `n` values sum to **205** specimens.
- **Result:** polarity-unified Fisher mean of the 26 site directions reproduces
  the published mean (D = 82.2°, I = 82.7°, k = 83, α95 = 3.1°, N = 26); the
  site-VGP Fisher mean gives the pole **46.5°N, 300.8°E (A95 5.8°)** — reported
  in the northern hemisphere as the antipode of the published 46.3°S, 121.4°E
  (dp 5.9, dm 6.1). GPMDB 6458.

`upload_magic` validates the result. This replaces the earlier extraction
(`Symons1991.csv` + the generic `Symons1991_csv_to_magic.py`), which lacked
locality coordinates, VGPs, the location-level pole, and the site-exclusion
flags.
