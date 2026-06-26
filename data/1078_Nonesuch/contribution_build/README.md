# Nonesuch Formation — data provenance

The Nonesuch DRM pole notebook (`pole_notebooks/1078_Nonesuch.ipynb`) loads its
site-level data **directly from the published MagIC contribution**
[earthref.org/MagIC/20614](https://earthref.org/MagIC/20614) (Slotznick et al.,
2024) via `pt.fetch_magic_contribution('20614', ...)`, which caches the
contribution as `../magic_contribution_20614.txt` so the notebook runs offline.

There is no local build script in this repository: contribution 20614 was
assembled and the pole computed in the source repository
**`github.com/Swanson-Hysell-Group/Nonesuch_Formation`** (`Code/Nonesuch_pmag.ipynb`).
This notebook reproduces that workflow from the published data.

## Source

- Slotznick, S. P., Swanson-Hysell, N. L., Zhang, Y., Clayton, K. E., Wellman,
  C. H., Tosca, N. J., & Strother, P. K. (2024). Reconstructing the
  paleoenvironment of an oxygenated Mesoproterozoic shoreline and its record of
  life. *GSA Bulletin*, 136(7–8), 2842–2864.
  doi:[10.1130/B37100.1](https://doi.org/10.1130/B37100.1). MagIC 20614.
- Age model: Fuentes, A. J., et al. (2025). Termination of Laurentia's rapid
  plate motion at the start of the Grenvillian Orogeny. *JGR Solid Earth*, 130,
  e2025JB031794. doi:10.1029/2025JB031794.

## Pole

The reported pole is the **detrital remanent magnetization (DRM)** from the
Potato River Falls section: combined detrital hematite (`hdt`, n = 105) +
detrital magnetite (`mt`, n = 77) specimen-site directions (tilt-corrected),
sharing a common mean, corrected for inclination shallowing by the E/I method
(Tauxe & Kent, 2004) with Kent-ellipse uncertainty propagated through the Pierce
et al. (2022) Monte Carlo → **6.6°N, 182.9°E (A95 ≈ 2.4°)**. The pigmentary
hematite (`hct`) carries a distinct CRM that defines a separate, nearby pole and
is excluded from the DRM pole. The notebook reproduces both from the MagIC
`sites`/`locations` tables.

This pole supersedes the legacy Henry et al. (1977) Nonesuch pole carried in the
prior compilation. (`data/Laurentia_poles.csv` is not modified by this notebook.)
