# Revisions and additions

## Additions

The following poles were added at the 2026 Iloranta Workshop:
- [ca. 780 Ma Gunbarrel LIP pole](pole_notebooks/780_Gunbarrel.ipynb) *This pole replaces previous Gunbarrel LIP results as it comes from localities that are not rotated relative to Laurentia.*
- [ca. 887 Ma Adirondack metamorphic anorthosite](pole_notebooks/887_Adirondack.ipynb) *This pole replaces the previously included Haliburton pole as the pole that is representative of the Grenville Loop. Given that this poles*
- [ca. 990 Ma Jacobsville Formation pole](pole_notebooks/990_Gunbarrel.ipynb)
- [ca. 1045 Ma upper Freda Formation pole](pole_notebooks/1045_Upper_Freda.ipynb)
- [ca. 1078 Ma Nonesuch Formation pole](pole_notebooks/pole_notebooks/1078_Nonesuch.ipynb)
- [ca. 1082 Cardenas Lava pole](pole_notebooks/1082_Cardenas.ipynb) *supersedes previous Cardenas pole with new data and many more sites from the lavas rather than also including undated intrusions that could vary in age*

The following poles were added at the 2022 Kringlerdalen Workshop:
- [ca. 757 Ma Chuar Group pole](pole_notebooks/755_Chuar_Group.ipynb)
- [ca. 1144 NW Ontario lamprophyre dikes pole](pole_notebooks/1144_Lamprophyre.ipynb)
- [ca. 1779 Ma East Central Minnesota Batholith pole](pole_notebooks/1779_East_Central_Minnesota_Batholith.ipynb)

## Changes

For the [ca. 1141 Ma Abitibi Dykes pole](pole_notebooks/1141_Abitibi.ipynb), the additional Abitibi site of Halls et al. (2015) was added to the mean calculation that brings the pole from N=7 to N=8. That was the same study that showed that the originally included A1 dike of Ernst and Buchan (1993) is actually part of the ca. 2167 Ma Biscotasing dike swarm.

In Elston et al. (2002), there was a decision for the Purcell lavas to calculate the pole only for lavas where the remanence was dominated by hematite. In doing so, dual-polarity remanences held by magnetite are discarded. It seems preferable to include all of the sites rather than only including the hematite ones that, while potentially early, are demonstrably secondary. We modify the Purcell Lavas pole to be calculated from all sites as documented in [ca. 1427 Purcell Lavas pole](pole_notebooks/1427_Purcell_Lava.ipynb). This changes the position and increases the uncertainty due to increased dispersion among the VGPs. It was also determined that the pole should be downgraded from an B to an A pole. More work should be done on the Purcell lavas to improve this pole and ideally isolate a magnetite-held thermal remanent magnetization.

Given the evidence presented in Ding et al. (2026), that the Tobacco Root Mountains are rotated relative to stable Laurentia, the previous Tobacco Root Gunbarrel dikes pole is no longer included in the compilation. The previous 1448±49 Ma Tobacco Root Dykes pole is no longer either as it.

Previously the Haliburton Highlands pole of the Grenville Province was included as a representative pole for the Grenville orogen with the thought that its cooling history was well-constrained. However, evidence presented in Zhang et al. (2026) reveal that there was both protracted high-grade metamorphism in portions of the orogen interior and that magnetite remanence in the Marcy Massif was acquired much later ca. 887 ± 23 Ma. Given that it is unlikely that Haliburton cooled below magnetic blocking temperatures by ca. 1015 Ma this pole is no longer included within the compilation.

The Euler rotation of Svalbard to Laurentia is uncertain to do multiple tectonic episodes including Caledonian translation and North Atlantic opening. As a result, at the 2026 Iloranta Workshop, these were moved to a separate off-craton block category being compiled by Prof. D.A.D. Evans.

For the [ca. 759 Ma Uinta Mountain Group pole](pole_notebooks/759_Uinta.ipynb), the preferred pole is updated to the site mean rather than the mean of sampling-locality means. The site mean is the Fisher mean of the 74 primary characteristic-remanence site virtual geomagnetic poles (1.9°N/160.6°E, A95 2.1°, N=74), recomputed at the site level from each site's tilt-corrected direction and coordinates. The previously reported pole was the mean of the 9 locality means of Weil et al. (2006) (0.8°N/161.3°E, A95 4.6°, N=9). The two agree within ~1.3°, but the site mean weights all sites equally and is better resolved. The depositional age was also refined following Dehler et al. (2023): the older age bound is the CA-ID-TIMS maximum depositional age of 766.3 ± 0.5 Ma, and the younger bound of ~730 Ma comes from correlation of the upper Uinta Mountain Group (Red Pine Formation) with the upper Chuar Group via shared vase-shaped microfossils.

For the [ca. 1092 Ma Portage Lake Volcanics pole](pole_notebooks/1092_Portage_Lake_Volcanics.ipynb), a typo in the published Swanson-Hysell et al. (2019) code inadvertently dropped two Books (1972) flows (PL157 and PL154) from the site set: adjacent string literals `'PL157' 'PL154'` were concatenated by Python into the single non-existent site id `'PL157PL154'`, so neither flow was matched. The published pole was therefore computed from N=78 cooling units (27.5°N/182.5°E, A95 2.3°). Restoring the two flows, and adding the Hnat et al. (2006) Greenstone-top flow H_PL05 that was missing from the earlier site table, yields the corrected pole of 27.1°N/183.0°E (A95 2.4°, N=80). The position shifts by only ~0.6°, but this notebook is now the version of record for the Portage Lake Volcanics pole.

## Thoughts

The discrepancy between the upper Belt and the ca. 1382 ± 2 Ma Zig-Zag Dal igneous province persists. The geochronology is sparse for the province with the only date being two multi-grain baddeleyite fractions. It would be worthwhile to develop more geochronology to see if this age is too young or too old.