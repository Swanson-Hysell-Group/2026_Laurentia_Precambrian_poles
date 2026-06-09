# Schroeder-Lutsen basalts — contribution build

How the site-level MagIC data for the Schroeder-Lutsen basalts (SLB) ca. 1090 Ma
pole are assembled, and how the pole is published into MagIC. The notebook
`pole_notebooks/1090_Schroeder_Lutsen_Basalts.ipynb` loads only the finished
`../sites.txt` / `../locations.txt`; everything here is provenance.

## Source tables (this folder)

| File | Source | Used for |
|---|---|---|
| `Fairchild2017_19680_sites.txt` | MagIC contribution [19680](https://earthref.org/MagIC/19680) (Fairchild et al., 2017), full sites table | the 40 Two Island River flow sites (`SLB01`–`SLB40`) |
| `TauxeKodama2009_sites.txt` | compiled sites of Tauxe & Kodama (2009) | the 10 Schroeder-Lutsen flow sites (`ns006`–`ns015`) |

## `SLB_build_contribution.py` → `../sites.txt`, `../locations.txt`

Selects the SLB sites and writes the clean per-pole contribution:

- **40 Fairchild et al. (2017)** sites: `location == 'Two Island River'`,
  `dir_comp_name == 'mag'` (magnetite, low-T component), `dir_tilt_correction == 100`.
- **10 Tauxe & Kodama (2009)** sites: `ns006`–`ns015` (the `nsl` / above-NSVG
  sequence), tilt-corrected; VGPs recomputed from the site mean directions.
- Each site is an individual lava flow → `result_type = 'i'`, one VGP per cooling
  unit. Per-site `citations` carry the two source DOIs (`10.1130/L580.1`,
  `10.1016/j.pepi.2009.07.006`).
- The 15 Books (1972) Schroeder sites used in the looser APWP `pole_means.csv`
  (N = 65 → 28.3/187.6) are **excluded**, following the published selection
  (Michipicoten precedent: reproduce the paper, not the broader compilation).
- The Fisher mean of the 50 VGPs reproduces the published pole: **26.9°N,
  188.0°E, A95 3.0, k 45.4, N 50** (Fairchild et al., 2017: 27.1°N / 187.8°E /
  A95 3.0). `upload_magic` validates the result.

The site selection and the non-Fisherian (two-cluster) VGP distribution are
documented in the `../locations.txt` pole `description`.

## Publishing the pole into MagIC 19680

Contribution 19680 archives Fairchild et al. (2017) sites/samples/specimens/
measurements but has **no location-level pole** for either of its two poles
(Michipicoten Island Formation, on the `Michipicoten Island` location; and the
Schroeder-Lutsen basalts, on the `Two Island River` location). The single updated
19680 carrying **both** poles is built by

```
data/1084_Michipicoten_Island_Formation/contribution_build/build_updated_19680.py
```

That script downloads 19680, appends the literature site means each pole needs
that are not already in 19680 (here: the 10 Tauxe & Kodama flow sites; the 40
Fairchild Two Island River sites are already present), sets the pole on each
location, passes all other tables through unchanged, and writes one validated
merged file into the `1084` folder. It reads the SLB pole and sites from this
pole's `../sites.txt` / `../locations.txt`. See that script's docstring for
details. The merged measurements file is large (~2.5 MB) and is not committed.
