# Mamainse Point Volcanics — contribution build (provenance)

Reproducible recipe for the four Mamainse Point paleomagnetic poles.

## Source (copied here for self-containment)

- `SH2014a_sites_source.txt` — the Mamainse Point site table of Swanson-Hysell et
  al. (2014a), *Geology*, doi 10.1130/G35271.1 (MagIC contribution **16333**),
  copied from `APWP_StratModels/data/pmag_compiled/Swanson-Hysell2014a/`. 99
  tilt-corrected site means with stratigraphic `height`, `vgp_lon`/`vgp_lat`.

## Four poles (reproduce `Keweenawan_pole_means.csv` and the prior Nordic compilation)

Poles are calculated from stratigraphic subsets of the VGPs (Swanson-Hysell et
al., 2009, 2014a; `APWP_StratModels/code/01_VGP_compilation.ipynb`):

| Pole | Stratigraphic height | Polarity zone | N | PLat/PLon | GPMDB |
|---|---|---|---|---|---|
| Mamainse lower R1 | < 600 m | Alona Bay reversed (older) | 24 | 49.5 / 227.0 | 9510 |
| Mamainse lower R2 | 1070–1350 m | Alona Bay reversed (younger) | 14 | 37.5 / 205.2 | 9511 |
| Mamainse Flour Bay | 1350–1810 m + 1860–2100 m | Flour Bay normal + reversed | 24 | 36.1 / 189.7 | 9512 |
| Mamainse upper N | > 2400 m | Portage Lake normal (above the Great Conglomerate) | 34 | 31.2 / 183.2 | 9513 |

Notes:
- All four poles are existing Nordic-compilation poles (GPMDB 9510–9513) and are
  exported to the Nordic summaries.
- The succession records **three geomagnetic reversals**; the four poles trace the
  progressive equatorward decrease in inclination interpreted as rapid plate
  motion (Swanson-Hysell et al., 2014a).
- The **Flour Bay** pole combines a normal and a reversed zone (antipodal VGPs);
  polarities are unified before averaging (`pmag.flip(combine=True)` /
  `compute_mean_pole(unify_polarity=True)`).
- The three transitional flows at ~1810–1860 m (heights 1820, 1833, 1858 — the
  Flour Bay normal→reversed transition) are excluded from the stable-polarity
  poles.

## Regenerate

```bash
python build_mamainse_contribution.py
```

writes `../sites.txt` (96 sites) and `../locations.txt` (4 poles) and validates
each pole against its target.
