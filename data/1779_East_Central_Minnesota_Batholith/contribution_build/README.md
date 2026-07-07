# ECMB pole — contribution build (provenance)

Reproducible recipe for the East-Central Minnesota Batholith (ECMB) ca. 1779 Ma
paleomagnetic pole and the update of MagIC contribution **20213**
(Swanson-Hysell et al., 2021, *Tectonics*, DOI 10.1029/2021TC006751).

The published contribution 20213 archives the full site / sample / specimen /
measurement data of the study, but its `locations` table carries **no
location-level pole result**. This build produces an updated contribution that is
**identical to 20213 except that the paleomagnetic pole is added to the
`locations` table**.

## Authoritative source

`../previous_MagIC/magic_contribution_20213.txt` — the full published
contribution 20213, downloaded from MagIC. Everything is read from this single
file, so the sites, samples, specimens, and measurements tables in the update are
guaranteed to be exactly those of the prior contribution (their raw text blocks
are copied through byte-for-byte).

## Build

```bash
python build_ECMB_contribution.py
```

`build_ECMB_contribution.py` reads contribution 20213 and writes:

- `../locations.txt` — the source location row with the ECMB `pole_*` result
  added (notebook runtime data).
- `../Swanson-Hysell2021_updated_ECMB_<DD.Mon.YYYY>.txt` — the combined MagIC
  upload file: the pole-bearing `locations` table followed by the `sites`,
  `samples`, `specimens`, and `measurements` tables copied verbatim from 20213.
  The `contribution` table is omitted; MagIC assigns the id, version, and
  contributor on upload. Upload this as an update of the 20213 contribution
  (DOI 10.1029/2021TC006751).

The script reproduces the pole from the 20213 site table as a self-check
(N = 23, 20.5°N / 265.8°E, A95 = 4.5, K = 45.6) before writing.

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
test** (the ca. 1096 Ma NW-trending dike NWD1 bakes NE-trending dike NED17)
records that the dike magnetization predates ca. 1096 Ma — encoded on the
location result as MagIC method code `ST-C-I` (Inverse contact test).

The pole age is 1779.1 ± 2.3 Ma (95% CI), from U-Pb dates bracketing the dikes
between the St. Cloud Granite (1781.44 ± 0.51 Ma) they intrude and the Richmond
Granite (1776.76 ± 0.49 Ma) they do not. This value is a Monte Carlo 95% CI over a
*uniform* emplacement window between the two granites (Swanson-Hysell et al., 2021),
not a Gaussian estimate — so it is encoded as `age` = 1779.1 with
`age_low`/`age_high` = 1776.8/1781.4 (the published 95% CI bounds) and **no
`age_sigma`**, since the MagIC data model defines `age_sigma` as a 1σ uncertainty
that does not apply to this distribution.

## Rich pole metadata on the location row

Beyond the `pole_*` columns, the location result carries the pole-relevant
metadata in structured, data-model-compliant columns (following the richness of
legacy location-only contributions such as Cleaver Dykes, MagIC 13653):

- **Location mean direction** (geographic; underlies the pole): `dir_dec` 179.9,
  `dir_inc` 76.7, `dir_alpha95` 2.5, `dir_k` 150.3, `dir_n_sites` 23,
  `dir_n_samples` 148, at `dir_tilt_correction` 0. (`dir_n_samples` is the
  data-model-valid home for the sample count; there is no `pole_n_samples`
  column.)
- **Field-test result** as a controlled-vocabulary column, complementing the
  `method_codes`. Only the test actually performed is recorded; columns for
  tests that were not done are left out rather than filled with `ND`:
  - `contact_test` = `IC+` — positive **inverse** contact test (the ca. 1096 Ma
    NW-trending dike NWD1 bakes NE-trending dike NED17), the structured
    counterpart of method code `ST-C-I`.

## History

An earlier build assembled the update from a different, less complete MagIC
contribution (**17072**) and was missing 5,018 measurement rows and 5 specimens
relative to 20213; it also used the non-canonical `tab \t` table marker, the
`ST-C` (baked contact) rather than `ST-C-I` (inverse contact) method code, and a
`pole_n_samples` column that is not in the MagIC data model. This build supersedes
it by sourcing everything from 20213 and correcting those issues.
