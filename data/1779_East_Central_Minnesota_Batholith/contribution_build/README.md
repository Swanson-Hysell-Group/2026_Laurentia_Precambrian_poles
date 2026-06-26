# ECMB pole — contribution build (provenance)

Reproducible recipe for the East-Central Minnesota Batholith (ECMB) ca. 1779 Ma
paleomagnetic pole and the enhancement of MagIC contribution **17072**
(Swanson-Hysell et al., 2021, *Tectonics*, DOI 10.1029/2021TC006751).

## Files

- `sites_17072_source.txt`, `locations_17072_source.txt`,
  `samples_17072_source.txt`, `specimens_17072_source.txt` — the published
  contribution 17072 tables, copied from
  `2021_ECMB/data/rockmag_pmag_MagIC/` for self-containment. The large
  `measurements.txt` (~4.7 MB) is **not** copied; pull it from the published
  contribution (or the source repo) if a full re-upload with measurements is
  needed.
- `ECMB_build_contribution.py` — selects the pole sites and writes the runtime
  data one level up: `../sites.txt` (full multi-component site table, unchanged)
  and `../locations.txt` (source location row with the `pole_*` result added).
- `build_updated_17072.py` — downloads the **full** published contribution 17072
  from MagIC (incl. the measurement table) and writes a single combined upload
  file `../Swanson-Hysell2021_updated_ECMB_<date>.txt` that is **identical to the
  published 17072 except the pole is added to the `locations` table**. This is the
  artifact to upload to MagIC as an update of 17072 (analogous to the
  `Fairchild2017_updated_…txt` in `data/1084_Michipicoten_Island_Formation/`).
  Needs network for the download; falls back to the published tables in the
  source repo if MagIC is unreachable. The online validation step is optional —
  if it can't reach MagIC the combined file is still written (validate it at
  upload time). The `citations` on the pole location row is set to the paper DOI
  (10.1029/2021TC006751); everything else, including all measurements, is
  unchanged from the published contribution.

## Pole definition (reproduces the paper)

The site table carries several demagnetization components per site
(lc = low-coercivity overprint, mc = medium-coercivity ChRM, hc = high-coercivity
Midcontinent Rift overprint, with thermal lt/mt/ht fits). The pole is the Fisher
mean of the **medium-coercivity (mc) VGPs** of the **NE-trending diabase dikes**,
keeping only site means with mc `dir_alpha95 < 8°` (the paper's cut). This drops
NED3, NED17, NED19, NED21; the ca. 1096 Ma NW-trending Midcontinent Rift dike
NWD1 and the granite sites are excluded. Result: **N = 23 sites, 148 samples,
pole 265.8°E / 20.5°N, A95 = 4.5°, K = 45.6** — matching Swanson-Hysell et al.
(2021) and the prior Nordic compilation (GPMDB 9970).

The dikes are near-vertical intrusions reported in geographic coordinates (no
tilt correction; `dir_tilt_correction = 0`). A positive **inverse baked-contact
test** (NWD1 bakes NE-trending dike NED17) records that the dike magnetization
predates ca. 1096 Ma — encoded on the location result as method code `ST-C`.

## To regenerate

```bash
python ECMB_build_contribution.py
```
