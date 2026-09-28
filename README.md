# Precambrian Paleomagnetic Poles for Laurentia

A curated, openly documented compilation of Precambrian paleomagnetic poles for Laurentia (the Precambrian craton at the core of North America and Greenland) from its ca. 1800 Ma amalgamation through the Ediacaran. Wherever possible, each pole is rebuilt from published site-level data, recomputed as a Fisher mean of site virtual geomagnetic poles (VGPs), and assessed against Nordic Paleomagnetism Workshop grades (A/B), the R-criteria of Meert et al. (2020), and the key-pole criteria of Buchan (2013).

**Rendered site:** <https://swanson-hysell-group.github.io/2026_Laurentia_Precambrian_poles/>

## Associated manuscript

This repository accompanies the manuscript:

> Swanson-Hysell, N.L., Zhang, Y., Slotznick, S.P., Zielinski, L.A., Eyster, A., and Piispa, E.J. *From sites to supercontinents: reconstructing Proterozoic Laurentia with a curated paleomagnetic pole compilation.* In review, *Canadian Journal of Earth Sciences*.

The pole table, apparent polar wander path (APWP) figures, paleolatitude figures, and paleomagnetic reconstructions in the manuscript are generated from the data and scripts in this repository.

## Context

The compilation builds on the synthesis of Laurentia's Precambrian paleogeography by [Swanson-Hysell (2021)](https://doi.org/10.1016/B978-0-12-818533-9.00009-6) and the quality-assessed global pole list of [Evans et al. (2021)](https://doi.org/10.1016/B978-0-12-818533-9.00007-2), both of which came out of the Nordic Paleomagnetism Workshops leading up to 2017. It incorporates the assessments made at subsequent workshops in Kringerdalen, Norway (2022) and Iloranta, Finland (2026).

Compiling site-level data in [MagIC](https://earthref.org/MagIC) format lets sites from multiple studies be combined, statistical field tests and secular-variation analyses be applied, and spatial and temporal uncertainties be propagated into APWP analyses. The most consequential updates relative to the 2017 compilation are in the late Mesoproterozoic to mid-Neoproterozoic: new Midcontinent Rift sedimentary poles and chronostratigraphy at the close of the Keweenawan Track, an updated Jacobsville Formation pole, a thermochronologic recalibration of the Grenville Loop apex to ca. 887 Ma, and revisions to the ca. 780 Ma Gunbarrel LIP pole.

The compilation is intended to be a living document that is updated as new paleomagnetic and geochronologic data become available.

## Repository contents

| Path | Contents |
| --- | --- |
| [pole_notebooks/](pole_notebooks/) | One Jupyter notebook per pole (named `<nominal age>_<unit>.ipynb`) that assembles the site-level data, recomputes the pole, applies field and statistical tests, and documents the geochronology, grade, and R-score. [pole_tools.py](pole_notebooks/pole_tools.py) holds the shared analysis and plotting functions. |
| [data/](data/) | Per-pole directories (e.g. `data/1078_Nonesuch/`) with MagIC-format tables (`locations.txt`, `sites.txt`, and, where available, `samples.txt`, `specimens.txt`, `measurements.txt`) along with any source files and conversion scripts. |
| [data/nordic_summaries/](data/nordic_summaries/) | Per-pole summary rows in the Nordic-workshop 71-column format written by each notebook; the combined compilation (`nordic_summaries_combined.csv`); Kent-distribution statistics for sedimentary poles; APWP spline inputs and outputs; and the LaTeX pole tables used in the manuscript. |
| `data/*.csv`, `data/*.xlsx` | Reference compilations: Evans et al. (2021), older (pre-1780 Ma) Laurentia poles, and the Torsvik et al. (2012) Phanerozoic compilation. |
| [data/geologic_provinces/](data/geologic_provinces/) | Laurentia basement provinces of Whitmeyer and Karlstrom (2007) as GeoJSON, used in the pole map and reconstructions. |
| [scripts/](scripts/) | Scripts that build the combined tables, figures, interactive pole map, Nordic-format Excel workbook, and API documentation, and that validate MagIC contributions before upload. |
| [scripts/sphereude/](scripts/sphereude/) | A self-contained Julia project that fits a spherical-spline APWP to the poles using [SphereUDE.jl](https://github.com/ODINN-SciML/SphereUDE.jl) (Sapienza et al., 2025). See its [README](scripts/sphereude/README.md). |
| [_static/](_static/) | Generated figures (APWP, paleolatitude, reconstruction ladder) and site assets. |
| [resources/](resources/) | The R-criteria scoring framework, a MagIC data-entry guide, the `pole_tools` API reference, and MagIC data-model and controlled-vocabulary files. |
| `index.md`, `compilation.md`, `paleolatitude.md`, `changes.md`, `pole_map.ipynb` | Top-level pages of the JupyterBook site. |
| `pole_checklist.md`, `geochronology_checklist.md` | Working trackers of per-pole progress (geochronology verification, notebook, MagIC contribution). |

## Reproducing the compilation and figures

Create the environment (Python 3.11 with PmagPy, cartopy, and JupyterBook 2 / MyST):

```bash
make install
mamba activate laurentia-poles
```

Then:

```bash
make tables    # rebuild the combined pole CSVs and manuscript LaTeX tables
make figures   # regenerate the pole map, compilation table, and paleolatitude figures
make build     # build the static HTML book into _build/html
make start     # serve the book locally with live reload
```

The APWP, paleolatitude-compilation, and reconstruction figures are produced by `scripts/build_apwp_figure.py`, `scripts/build_paleolatitude_compilation_figure.py`, and `scripts/build_reconstruction_figure.py`. The spline fits these rely on are committed in `data/nordic_summaries/`, so the figures can be regenerated without Julia; refitting the spline requires the Julia project in `scripts/sphereude/`.

Pushes to `main` rebuild and deploy the site to GitHub Pages via [.github/workflows/deploy.yml](.github/workflows/deploy.yml).

## License

Content is released under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/). Please cite the associated manuscript and the original data sources credited in each pole notebook.
