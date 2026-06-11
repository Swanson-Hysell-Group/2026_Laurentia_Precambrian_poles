# Schroeder-Lutsen basalts — contribution build

How the site-level MagIC data for the Schroeder-Lutsen basalts (SLB) ca. 1090 Ma
pole are assembled, and how the pole is published into MagIC. The notebook
`pole_notebooks/1090_Schroeder_Lutsen_Basalts.ipynb` loads only the finished
`../sites.txt` / `../locations.txt`; everything here is provenance.

## Source tables (this folder)

| File | Source | Used for |
|---|---|---|
| `Fairchild2017_published_site_means.csv` | Fairchild et al. (2017) study analysis file (`pmag_results.csv` in `Swanson-Hysell-Group/2017_Late_Rift`) | the **published site mean directions** of the 40 Two Island River flows (`SLB01`–`SLB40`; mag/hem × tilt-0/100) |
| `Fairchild2017_19680_sites.txt` | MagIC contribution [19680](https://earthref.org/MagIC/19680) (Fairchild et al., 2017), full sites table | the site **coordinates and method codes** for the 40 flow sites |
| `TauxeKodama2009_sites.txt` | compiled sites of Tauxe & Kodama (2009) | the 10 Schroeder-Lutsen flow sites (`ns006`–`ns015`) |

## `SLB_build_contribution.py` → `../sites.txt`, `../locations.txt`

Selects the SLB sites and writes the clean per-pole contribution:

- **40 Fairchild et al. (2017)** sites: magnetite (`mag`), tilt-corrected. The
  site **mean directions** are the **published study values** (from
  `Fairchild2017_published_site_means.csv`); coordinates and method codes are
  carried from the 19680 sites table.
- **10 Tauxe & Kodama (2009)** sites: `ns006`–`ns015` (the `nsl` / above-NSVG
  sequence), tilt-corrected; VGPs recomputed from the site mean directions.
- Each site is an individual lava flow → `result_type = 'i'`, one VGP per cooling
  unit. Per-site `citations` carry the two source DOIs (`10.1130/L580.1`,
  `10.1016/j.pepi.2009.07.006`).
- The 15 Books (1972) Schroeder sites used in the looser APWP `pole_means.csv`
  (N = 65 → 28.3/187.6) are **excluded**, following the published selection
  (Michipicoten precedent: reproduce the paper, not the broader compilation).
- The Fisher mean of the 50 VGPs reproduces the published pole exactly: **27.1°N,
  187.8°E, A95 3.0, k 46.4, N 50** (Fairchild et al., 2017: 27.1°N / 187.8°E /
  A95 3.0). `upload_magic` validates the result.

### Why the published study means (not the 19680-archived means)

The site means archived in MagIC contribution 19680 differ from the published
study values for **8 of the 40 Two Island River sites** (`SLB08`, `SLB10`,
`SLB15`, `SLB23`, `SLB27`, `SLB28`, `SLB31`, `SLB32`). At each of those sites the
19680 mean **retains one additional specimen** (n + 1) that the study rejected as
a directional outlier — a near-antipodal / excursional sample direction (e.g.
`SLB08.5a` = 153°/−13° versus the ~282°/55° flow mean). Including it collapses the
site precision (e.g. `SLB08` *k* 617 → 6, `SLB32` *k* 271 → 8) and shifts the site
mean by ~0.6–8°. Built from the 19680-archived means the pole comes out at
**26.9°N / 188.0°E** (a ~0.2° offset); built from the published study means it
reproduces **27.1°N / 187.8°E**. The excluded sample at each site is identified
reproducibly as the one whose removal restores the published site mean (direction
+ *k* + *n*) exactly (off by < 0.05°).

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
location, and writes one validated merged file into the `1084` folder. It reads
the SLB pole and sites from this pole's `../sites.txt` / `../locations.txt`. See
that script's docstring for details. The merged measurements file is large
(~2.5 MB) and is not committed.

As part of building the update it also **corrects the 8 Two Island River site
means** described above to the published study values (all four mag/hem ×
tilt-0/100 rows per site, from `Fairchild2017_published_site_means.csv`) and
**flags the 8 rejected outlier specimens** `result_quality='b'`, so the archived
site means are consistent with the specimen data and reproduce the published pole.
This makes the updated contribution the corrected version of 19680 to upload.
