# Coldwell Complex — contribution build

How the site-level MagIC data for the ca. 1107 Ma Coldwell Complex pole are
assembled. The notebook `pole_notebooks/1107_Coldwell.ipynb` loads the finished
`../sites.txt`; `../locations.txt` carries the pole. Everything here is provenance.

## Source

Kulakov, E. V., Smirnov, A. V., & Diehl, J. F. (2014). Paleomagnetism of the
~1.1 Ga Coldwell Complex (Ontario, Canada): Implications for Proterozoic
geomagnetic field morphology and plate velocities. *Journal of Geophysical
Research: Solid Earth*, 119(12), 8633–8654.
doi:[10.1002/2014JB011463](https://doi.org/10.1002/2014JB011463). (`Kulakov2014a`.)

There is **no measurement-level MagIC contribution** for this study, so the
contribution is built from the published **Table 2** site mean directions,
digitized in `Coldwell_build_contribution.py`.

Age: U-Pb zircon/baddeleyite on the eastern gabbro, 1108 ± 1 Ma and 1107 +5/−1 Ma
(Heaman & Machado, 1992). GPMDB 9838.

## `Coldwell_build_contribution.py` → `../sites.txt`, `../locations.txt`

- Writes the **40 accepted sites (238 samples)**, grouped (after Lewchuk & Symons,
  1990) into three intrusive centers via the `description` field and `dir_polarity`:
  - **Center A** (eastern gabbro / ferroaugite syenite): 14 sites, **reversed** (`dir_polarity='r'`).
  - **Center B** (western gabbro and syenite): 10 sites, **normal** (`dir_polarity='n'`).
  - **Center C** (central biotite gabbro / nepheline syenite / syenite, Geordie Lakes area): 16 sites, **reversed** (`'r'`).
- Polarity convention: Table 2 lists inclinations as positive magnitudes. The
  reversed Center A/C ChRM is stored with **negative inclination** (the prior
  compilation lists I = −63.7); Center B normal with positive inclination. Per-site
  VGPs are computed with `pmag.dia_vgp` from the signed directions (reversed →
  southern VGP), and the pole is the polarity-unified Fisher mean.
- Geographic coordinates (`dir_tilt_correction = 0`; intrusive, no bedding
  correction). Longitudes converted from the published °W to 0–360 °E.

## Pole of record

The preferred pole (**Pole CCr** of Kulakov et al., 2014) is the Fisher mean of
the **30 reversed-polarity site VGPs of Centers A and C** (whose group mean
directions are statistically indistinguishable: combined D = 114.8°, I = −63.7°,
α95 = 3.6°, k = 54): **47.2°N, 206.5°E (A95 = 4.8°, K = 31, N = 30)**, reported in
the northern (polarity-unified) hemisphere — reproducing the published value
exactly. `locations.txt` carries this pole. The 10 normal-polarity Center B sites
(D = 298.0°, I = 56.9°) are retained in the contribution and used in the
notebook's reversal test. `data/Laurentia_poles.csv` (which already lists
47.2/206.7) is not modified by this notebook.
