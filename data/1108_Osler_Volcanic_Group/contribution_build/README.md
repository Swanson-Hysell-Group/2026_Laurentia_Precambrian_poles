# Osler Volcanic Group — contribution build (provenance)

Reproducible recipe for the three Osler Volcanic Group paleomagnetic poles. This
follows the user's suggested approach of enhancing the Swanson-Hysell et al.
(2014b) Osler contribution with the additional Swanson-Hysell et al. (2019) sites
and the Halls (1974) sites, with the three poles in the locations table.

## Sources (copied here for self-containment)

Copied from `APWP_StratModels/data/pmag_compiled/`, which mirrors the published
MagIC contributions:

- `Halls1974_sites_source.txt` — Halls (1974), CJES, doi 10.1139/e74-113;
  Nipigon Strait reversed (25 sites) and normal (5 sites) Osler flows.
- `SH2014b_sites_source.txt` — Swanson-Hysell et al. (2014b), G-cubed,
  doi 10.1002/2013gc005180; Simpson Island stratigraphic section reversed flows.
- `SH2019a_AgatePoint_sites_source.txt` (+ `..._samples_source.txt`) —
  Swanson-Hysell et al. (2019), GSA Bulletin, doi 10.1130/b31944.1; Agate Point
  reversed flows and Puff Island normal flows.

## Three poles (reproduce APWP_StratModels `Keweenawan_pole_means.csv`)

Following `APWP_StratModels/code/01_VGP_compilation.ipynb`:

| Pole | Construction | N | PLat/PLon | GPMDB |
|---|---|---|---|---|
| Osler reverse lower (R1) | SH2014b flows height < 1041 m | 30 | 40.9 / 218.6 | 9515 |
| Osler reverse middle | SH2014b flows 1041–2082 m | 20 | 42.7 / 211.3 | 10017 |
| Osler reverse upper (R2) | SH2014b flows height > 2082 m + Halls (1974) reversed + Agate Point (2019) | 64 | 42.3 / 203.4 | 9514 (updated) |
| Osler normal | Halls (1974) 5 normal sites combined into 2 flows (1/2/5; 3/4) + 2 Puff Island flows (2019) | 4 | 32.0 / 171.9 | new |

Notes:
- All **four** poles are written to `locations.txt`. The **three
  reversed-polarity poles** (lower 9515, middle 10017, upper 9514) are existing
  Nordic-compilation poles and are exported to the Nordic summaries
  (`1108_/1107_/1105_Osler_reverse_*.csv`).
- The **Osler normal** pole is computed and plotted in the notebook and kept in
  `locations.txt`, but is **not** exported to the Nordic compilation because it
  has too few sites (only four cooling units; N = 4).
- The Halls (1974) normal pole combines its 5 "sites" into the **2 cooling-unit
  flows** they actually represent (field mapping of Swanson-Hysell & Fairchild,
  2014), as `result_type='a'` with the combination documented per site.
- The Agate Point reversed VGPs are antipodal to the rest; polarities are unified
  before averaging (`pmag.flip(combine=True)` / `compute_mean_pole(unify_polarity=True)`).
- `location` is set to the pole grouping so the notebook can split sites by pole.

## Regenerate

```bash
python build_osler_contribution.py
```

writes `../sites.txt` (98 sites) and `../locations.txt` (3 poles) and validates
each pole against the target.
