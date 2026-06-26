# Lower Freda Formation — data provenance

The lower Freda pole notebook (`pole_notebooks/1075_Lower_Freda.ipynb`) loads its
specimen-level data **directly from the published MagIC contribution**
[earthref.org/MagIC/16450](https://earthref.org/MagIC/16450) via
`pt.fetch_magic_contribution('16450', ...)`, which caches the contribution as
`../magic_contribution_16450.txt` and unpacks `../sites.txt`, `../samples.txt`,
and `../specimens.txt` so the notebook runs offline.

MagIC 16450 archives the paleomagnetic directions of **Henry, Mauk & Van der Voo
(1977)** — 58 samples / 127 specimens of Nonesuch Shale and Freda Sandstone,
treated by thermal, alternating-field, and chemical demagnetization. The
specimen-level directions were originally reported in the **Henry (1976)** MSc
thesis (University of Michigan). This notebook uses only the **Freda Sandstone**
high-temperature specularite DRM specimens (`description == 'Freda sandstone'`,
`dir_comp == 'PRIM'`, tilt-corrected, `dec > 240°`; n = 85).

There is no local build script: the contribution was assembled on MagIC, and the
lower Freda pole was analyzed by Fuentes et al. (2025) in the source repository
**`github.com/Swanson-Hysell-Group/Upper_Freda_Pole`** (`code/Freda_pmag.ipynb`).
This notebook reproduces that workflow from the published data.

## Sources

- Henry, S. G., Mauk, F. J., & Van der Voo, R. (1977). Paleomagnetism of the
  upper Keweenawan sediments: Nonesuch Shale and Freda Sandstone. *Canadian
  Journal of Earth Sciences*, 14, 1128–1138.
  doi:[10.1139/e77-103](https://doi.org/10.1139/e77-103). MagIC 16450.
- Henry, S. G. (1976). *Paleomagnetism of the upper Keweenawan sediments: the
  Nonesuch Shale and Freda Sandstone.* MSc thesis, University of Michigan
  (specimen-level data source).
- Age model and SVEI re-analysis: Fuentes, A. J., Fairchild, L. M., Hodgin,
  E. B., Alemu, T., & Swanson-Hysell, N. (2025). Termination of Laurentia's rapid
  plate motion at the start of the Grenvillian Orogeny. *JGR Solid Earth*, 130,
  e2025JB031794. doi:[10.1029/2025JB031794](https://doi.org/10.1029/2025JB031794).

## Pole

The reported pole is the **high-temperature specularite detrital remanent
magnetization (DRM)** from the Presque Isle River exposures (Henry's HR1/HR2
localities): the Fisher mean of the 85 tilt-corrected Freda DRM specimen
directions (D = 271.3°, I = +0.7°, k = 31, α95 = 2.8°) gives a pole at
**1.2°N, 179.7°E (A95 2.4°)** — reproducing Henry et al.'s (1977) published
179.5°E/1.2°N and the legacy Nordic-compilation "Freda Sandstone" pole
(2.2°N/179.0°E, GPMDB 2051). The DRM passes Fuentes et al.'s (2025) regional
bootstrap fold test on the combined Freda directions (95–102% unfolding),
demonstrating pre-tilting acquisition; because Grenvillian inversion tilted the
Freda at ca. 1050–1000 Ma (Douglas Fault calcite 1,049 ± 13 Ma; Hodgin et al.,
2024), this tightly brackets the specularite DRM as the primary ca. 1075 Ma
depositional magnetization.

The lower Freda lies at the **paleo-equator**, so inclination-shallowing
correction is negligible: Fuentes et al. (2025) applied the same SVEI flattening
uncertainty propagation as for the upper Freda and report the corrected pole at
**1.3°N, 179.8°E (A95 2.6°)** — reproduced here deterministically by unsquishing
the DRM inclinations (the operation behind Fuentes' `Henry_ht_pole`).

This is the *lower Freda* member of the Oronto Group APWP sequence (Nonesuch →
lower Freda → upper Freda → Jacobsville) along the Keweenawan Track.
(`data/Laurentia_poles.csv` is not modified by this notebook.)
