# Lake Shore Traps — contribution build

How the site-level MagIC data for the Lake Shore Traps (LST) ca. 1085 Ma pole
were assembled and published. The notebook
`pole_notebooks/1086_Lake_Shore_Traps.ipynb` loads the finished data directly
from the published contribution [MagIC 20696](https://earthref.org/MagIC/20696)
(`pt.fetch_magic_contribution('20696', ...)`); everything here is provenance.

## Source contributions (this folder)

| File | Source study | DOI |
|---|---|---|
| `magic_contribution_16334_Diehl_Haig_1994.txt` | Diehl & Haig (1994) | `10.1139/e94-034` |
| `magic_contribution_16335_Kulakov_2013.txt` | Kulakov et al. (2013) | `10.1139/cjes-2013-0003` |

Both are legacy site + location-level contributions (no measurement-level data),
retrieved from MagIC by DOI.

## `LST_build_contribution.py` → MagIC 20696

A one-time build that merged the Diehl & Haig (1994) sites **into** the Kulakov
et al. (2013) contribution (the follow-up "master" study), producing the single
enhanced contribution now published as 20696, which carries the location-level
LST mean pole. The notebook pulls solely from 20696. Enhancements applied to the
legacy tables (see the script docstring for the full list):

- `citations` set to each study's DOI (the legacy tables said "This study"), so
  per-site provenance is preserved;
- `method_codes` corrected to `LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM` (AF + thermal,
  PCA, Fisher means), confirmed from both papers;
- site-level `geologic_classes`/`geologic_types`/`lithologies`/`age` added;
- `vgp_dp`/`vgp_dm` computed; longitudes converted to 0–360°E;
- a location-level pole result written to `locations.txt`.

**Duplicate flow:** site `LST28` was measured by both studies. The Diehl & Haig
(1994) measurement is used in the pole (chronological priority); the Kulakov
et al. (2013) re-measurement is retained as a separate site `LST28_K2013` flagged
`result_quality='b'` and excluded from the pole. The pole therefore uses **50
unique cooling-unit site means** (51 site rows: 50 good + 1 excluded duplicate).

Result: mean pole **23.1°N, 185.7°E, A95 3.3, k 38.8, N 50**, reproducing the
previously compiled Lake Shore Traps pole. Age adopted from the U-Pb
(206Pb/238U, CA-ID-TIMS) andesite weighted-mean date of 1085.57 ± 0.5 Ma
(Fairchild et al., 2017). The pole carries a positive Keweenawan conglomerate
test (`ST-G`; Palmer, Halls & Pesonen, 1981).
