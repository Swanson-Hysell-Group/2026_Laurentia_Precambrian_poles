"""Build the MagIC contribution for the Coldwell Complex (Kulakov et al., 2014).

Prepared for the compilation of Laurentia's Precambrian paleomagnetic poles of
Swanson-Hysell et al. (2026), which recalculates each pole from site-level data.

Kulakov, E. V., Smirnov, A. V., & Diehl, J. F. (2014). Paleomagnetism of the ~1.1 Ga
Coldwell Complex (Ontario, Canada): Implications for Proterozoic geomagnetic field
morphology and plate velocities. J. Geophys. Res. Solid Earth 119, 8633-8654.
doi:10.1002/2014JB011463

Geochronology: Heaman, L. M. & Machado, N. (1992). Timing and origin of midcontinent
rift alkaline magmatism, North America: evidence from the Coldwell Complex. Contrib.
Mineral. Petrol. 110, 289-303. doi:10.1007/BF00310744. Uncertainties are 2 sigma
(p. 295 and the Fig. 3 caption); age_sigma is 1 sigma.

Only site means are published ('available by request from the corresponding author';
no supporting information), so this is a site-level contribution. There is no earlier
MagIC contribution for this paper; 12757/17396 are Lewchuk & Symons (1990), a different
study.

Sources and checks
------------------
- Table 2 (site means) and Table 1 (lithology and rock magnetism) were transcribed
  twice, independently, from the rendered pages; the two transcriptions agree on every
  coordinate, count, direction, k and alpha95 of all 58 sites.
- Site VGPs are calculated from each site's direction and coordinates and stored to
  0.01 deg, in the polarity of the site direction; the VGP listed in Table 2 (to 0.1
  deg, in normal polarity) is given in each site description. Averaging the calculated
  VGPs reproduces every pole of the paper exactly at its printed precision, whereas
  averaging the listed values, rounded to 0.1 deg, shifts the Center A and Center B
  pole latitudes by 0.01 deg across a rounding boundary. The listed and calculated VGPs
  agree within 0.3 deg except at GLG and CLD8 (the calculated VGPs reproduce the
  authors' Center C and A+C poles) and CLD11 (0.4 deg); the check is asserted below.
  CLD1's VGP is calculated at the coordinates the listed VGP follows from (those of
  the first author's dissertation), as the Center A pole requires.
  Pole rows report poles in the northern hemisphere, as in the paper.
- CLD3: the site longitude is read from Fig. 2 (-86.33, geographic_precision 0.01),
  which places it with the other eastern Center A sites and reproduces the listed VGP
  within 0.03 deg (Table 2 lists -86.673).
- Poles and mean directions on the location rows are recalculated from the site rows;
  they reproduce the published values within 0.1 deg (asserted below).

Structure and interpretation
----------------------------
- One location, 'Coldwell Complex': an umbrella row, the five poles of Table 4 (CCA,
  CCn, CCC, CCr = A+C, and A+B+C), and the two Center B variants the paper reports.
  Pole CCr is the pole used in the Swanson-Hysell et al. (2026) compilation.
- Reversal test: recalculated from the Table 2 site means with simulated critical
  angles, Center C vs Center B passes (classification B), while the angles between the
  Center A and the combined A+C means and the Center B mean (9.3 and 7.0 deg, as in
  section 3.3) exceed their critical angles. The positive C vs B test shows a symmetrically reversing
  field while the complex was emplaced, so reversal_test = 'Rb' is recorded on all
  five Table 4 poles, with the recalculated numbers in the descriptions.
- Site and center ages are ages of the rocks (U-Pb, Heaman & Machado 1992), not
  interpreted magnetization ages: Center A 1108 +/- 1 Ma (also CLD9, whose direction is
  a present-field overprint), Center C 1109 +8/-4 Ma, Center B 1107 +9/-5 Ma (the
  granite of the western part). Pole CCr and the A+B+C pole: 1107 Ma (1105-1109 Ma).
  CCr averages reversed-polarity bodies that record paleolatitudinal change within the
  Alona Bay reversed-polarity chron, so it is given the midpoint of an interval spanning
  their dates rather than the 1108 +/- 1 Ma date of a single body (which Table 4
  assigns to it).
"""
import contextlib
import io
import os
import sys
from datetime import date

import matplotlib
matplotlib.use('Agg')  # ipmag.common_mean_watson calls plt.show(), which blocks on a GUI backend
import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(OUT))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

LOCATION = 'Coldwell Complex'
HM1992 = '10.1007/BF00310744'
GPMDB_RESULT = 9838

# Table 2: site, center, lat (N), lon (printed negative = W), N (samples used),
# n (samples measured), method, D, I, a95, k, VGP lat, VGP lon (printed in normal
# polarity). '-' = not given.
TABLE_2 = """
CLD1    A 48.713 -86.322  3  7 AF+T 116.0 -61.5  7.4  187   46.2 200.6
CLD3    A 48.719 -86.673  8  9 AF+T 125.7 -73.8  8.6   37   57.5 224.3
CLD4    A 48.727 -86.328 14 17 AF+T 103.7 -66.4  4.4   76   42.0 214.4
CLD5    A 48.729 -86.329  7 11 AF+T 116.9 -73.3  4.4   76   52.9 224.1
CLD6    A 48.742 -86.334  9 10 AF+T 106.3 -67.2  8.5   34   44.0 214.5
NMG1    A 48.847 -86.484  9 10 AF+T  99.1 -51.9 12.3   32   29.5 200.4
SL10    A 48.768 -86.374  5  6 AF+T 118.5 -70.7  6.4  114   52.8 217.2
SL14    A 48.780 -86.420  6 11 AF+T 125.5 -49.2 11.0   32   45.1 179.7
SL15    A 48.773 -86.393  5 12 T    118.4 -68.5 11.0   43   51.7 212.2
SL33    A 48.785 -86.428  8  9 AF+T 118.6 -66.7  3.9  174   50.9 208.3
SL34    A 48.726 -86.327  4  7 AF+T 138.5 -75.6  9.1   73   63.4 231.1
SL5-6   A 48.732 -86.330  8  8 AF+T  98.0 -63.5  6.7   61   36.7 213.0
SL7     A 48.756 -86.310  5  6 AF+T 113.7 -68.8  6.6  110   49.2 214.6
SL7-8   A 48.755 -86.311  4  6 AF   127.4 -62.4 11.0   52   54.0 196.0
CLD2    A 48.708 -86.303  0 13 AF+T     -     -    -     -     -     -
CLD9    A 48.768 -86.373  3 12 AF+T 351.6  63.4 13.6   55     -     -
SL8     A 48.754 -86.310  0  7 AF+T     -     -    -     -     -     -
NMG2    A 48.800 -86.314  8  8 AF+T 129.4 -48.1 35.2  2.7     -     -
NMG3    A 48.787 -86.316  7  8 AF+T 171.6 -44.0 58.4  1.7     -     -
CCW1    B 48.814 -86.673  5 11 T    267.0  53.2 14.3   23   22.9 209.0
CCW2    B 48.817 -86.684  8 11 AF+T 293.3  51.0  5.9   77   38.1 190.1
CLD7    B 48.770 -86.549  4 12 AF   318.0  57.3  9.3   73   58.1 179.6
SL13E   B 48.796 -86.453  8  8 T    303.5  62.6  4.8  120   51.6 198.3
SL18    B 48.799 -86.651  6  8 AF+T 299.2  49.8  6.0  107   41.3 184.8
SL19    B 48.796 -86.646  5  7 AF+T 295.5  58.5 10.0   46   44.1 196.7
SL20    B 48.800 -86.636  6  8 AF+T 308.2  59.3 14.0   21   52.8 189.9
SL37    B 48.802 -86.626  6  8 T    321.3  54.4  4.9  154   58.5 171.9
SL47    B 48.817 -86.682  7  8 AF   287.7  58.7 10.8   28   39.2 201.5
SL49    B 48.797 -86.651  6  8 AF   289.6  55.1  5.2  142   38.2 196.5
MD1     B 48.798 -86.654  8  8 AF+T   3.5   0.2 62.8  1.5     -     -
MD2     B 48.797 -86.656 12 15 AF+T 313.3   8.4 43.3  1.6     -     -
MD3     B 48.795 -86.658  9 12 AF+T  54.2 -43.5 53.2  1.7     -     -
MD4     B 48.792 -86.661  3 12 AF+T 349.0  12.2 84.5  1.7     -     -
SL26    B 48.807 -86.662  0 12 AF+T     -     -    -     -     -     -
SL48    B 48.798 -86.640  4  8 AF+T 316.3  27.9 27.7    9     -     -
CLD11   C 48.784 -86.596  4  7 AF   129.4 -49.2  7.7  107   47.9 176.7
GLS     C 48.824 -86.485  5  5 AF+T  92.1 -57.2 14.9   23   28.7 209.3
GLG     C 48.821 -86.485  4  6 AF+T 131.5 -43.8  8.4   90   46.0 229.8
GLS2    C 48.824 -86.485  3  4 AF   110.2 -49.6  7.9  162   35.3 191.2
CLD8    C 48.765 -86.525  3  6 AF   131.2 -64.5  6.3  258.7 44.0 214.5
SL11    C 48.792 -86.487  5 10 AF+T  77.1 -68.9  4.8  205   30.4 229.9
SL13    C 48.793 -86.463  6  6 AF+T 120.4 -60.5  4.4  191   48.4 196.7
SL14-13 C 48.779 -86.420  8  9 AF+T  89.8 -69.1  3.1  284   36.6 224.5
SL21    C 48.789 -86.611  6  8 AF+T 120.7 -57.5  7.0   77   46.8 192.2
SL25    C 48.793 -86.501  4  6 AF   110.6 -57.0  1.7 2315   40.0 198.0
SL28    C 48.802 -86.626  6 12 AF+T 118.6 -45.4  7.7   64   38.4 181.7
SL32    C 48.795 -86.495  8  8 AF+T  99.4 -75.2  6.4   67   45.7 232.2
SL41    C 48.797 -86.445  6  8 AF   132.2 -67.9  8.3   56   59.7 206.0
SL44    C 48.764 -86.526  3  9 AF+T 114.0 -59.7  7.5   56   43.8 199.3
SL45    C 48.772 -86.513  8 10 AF   108.5 -69.7  9.1   38   46.8 218.0
SL46    C 48.795 -86.497  3  9 AF   114.1 -76.3  5.9  292   52.6 232.3
GLG2    C 48.824 -86.485  0  6 AF+T     -     -    -     -     -     -
CLD10   C 48.780 -86.581  6  8 AF+T 262.8 -66.7 58.8  1.8     -     -
DM1     C 48.789 -86.480 13 13 AF+T 208.3 -28.1 51.3  1.5     -     -
SL23    C 48.772 -86.556  2  8 AF+T     -     -    -     -     -     -
SL30    C 48.794 -86.489  0  7 AF+T     -     -    -     -     -     -
SL31    C 48.792 -86.491  2  8 AF+T     -     -    -     -     -     -
SL40    C 48.792 -86.468  0  8 AF+T     -     -    -     -     -     -
"""

# Table 1: site | lithology as printed | NRM range (A/m) | Hc (mT) | Hcr (mT) | Mrs/Ms.
# Table 1 (and section 3.1, Fig. 2) name GLS1 the site that Table 2 names GLS2.
TABLE_1 = """
CLD1    | Rhyolite       | 0.018-0.355   |  4.77 | 14.22 | 0.058
CLD3    | E. gabbro      | 0.069-3.211   | 17.43 | 33.4  | 0.210
CLD4    | E. gabbro      | 0.074-28.80   |  7.99 | 18.44 | 0.071
CLD5    | E. gabbro      | 0.122-1.13    |  9.43 | 22.34 | 0.086
CLD6    | Aug. syenite   | 0.104-0.271   |  8.66 | 24.65 | 0.072
NMG1    | E. gabbro      | 0.305-45.88   |  6.38 | 18.2  | 0.057
SL10    | Aug. syenite   | 0.044-0.184   |  4.56 | 11.7  | 0.059
SL14    | Aug. syenite   | 0.0034-0.396  |  7.17 | 14.86 | 0.087
SL15    | Aug. syenite   | 0.146-0.643   |  6.78 | 16.25 | 0.072
SL33    | Basalt         | 2.305-21.58   | 14.31 | 22.59 | 0.183
SL34    | Aug. syenite   | 0.591-36.43   |  7.32 | 16.94 | 0.043
SL5-6   | E. gabbro      | 6.12-29.90    | 13.56 | 24.4  | 0.168
SL7     | E. gabbro      | 1.621-2.418   |  7.51 | 13.07 | 0.097
SL7-8   | E. gabbro      | 0.607-3.06    |  7.88 | 29.06 | 0.091
CLD2    | E. gabbro      | 0.0009-0.0019 | 16.6  | 23.18 | 0.252
CLD9    | Aug. syenite   | 0.405-0.757   | 11.43 | 19.32 | 0.237
SL8     | E. gabbro      | 5.95-10.11    |  1.39 |  4.24 | 0.013
NMG2    | E. gabbro      | 0.031-2.586   | 11.56 | 27.63 | 0.110
NMG3    | E. gabbro      | 0.36-1.98     |  3.38 | 13.25 | 0.036
CCW1    | W. gabbro      | 0.004-0.098   | 10.92 | 13.44 | 0.131
CCW2    | Pillowed lava  | 0.002-0.024   | 10.23 | 17.26 | 0.154
CLD7    | Mafic dyke     | 0.005-0.229   |  4.87 | 15.89 | 0.207
SL13E   | Mafic dyke     | 0.977-13.25   |  9.92 | 19.21 | 0.112
SL18    | W. gabbro      | 4.748-30.06   |  4.93 |  9.38 | 0.101
SL19    | W. gabbro      | 0.314-2.00    |  5.05 | 10.46 | 0.051
SL20    | W. gabbro      | 0.604-59.47   |  7.33 | 12.43 | 0.124
SL37    | W. gabbro      | 1.326-33.72   |  3.44 |  6.93 | 0.047
SL47    | Pillowed lava  | 0.834-9.591   |  9.54 | 23.10 | 0.154
SL49    | W. gabbro      | 6.864-88.68   |  4.75 | 28.56 | 0.025
MD1     | Gabbro         | 0.564-7.24    |  4.51 |  8.94 | 0.058
MD2     | Gabbro         | 0.03-30.12    |  3.00 |  7.84 | 0.034
MD3     | Gabbro         | 6.32-70.01    |  5.83 | 13.52 | 0.019
MD4     | Gabbro         | 0.777-2.826   |  7.45 | 14.86 | 0.095
SL26    | Aug. syenite   | 0.0001-0.343  | 10.55 | 20.83 | 0.114
SL48    | Neph. syenite  | 0.061-1.993   |  9.43 | 23.14 | 0.097
CLD11   | Q-tz syenite   | 0.017-0.04    |  1.88 |  5.09 | 0.016
GLG     | Gabbro         | 4.787-34.15   |  6.57 | 16.29 | 0.074
GLS     | Q-tz syenite   | 0.097-0.438   |  4.36 | 10.12 | 0.051
GLS2    | Q-tz syenite   | 0.01-0.067    |  5.18 | 11.49 | 0.042
CLD8    | Neph.syenite   | 0.102-0.994   |  2.42 |  7.45 | 0.044
SL11    | Amph. syenite  | 0.026-1.11    |  3.94 |  8.01 | 0.051
SL13    | Amph. syenite  | 0.043-0.391   |  4.53 | 13.61 | 0.099
SL14-13 | Augite syenite | 0.665-5.797   |  6.99 | 13.77 | 0.082
SL21    | Q-tz syenite   | 0.059-1.511   |  8.75 | 33.2  | 0.062
SL25    | Neph. syenite  | 0.066-0.196   |  2.18 |  5.72 | 0.023
SL28    | Amph. syenite  | 0.021-2.724   |  2.9  | 13.43 | 0.026
SL32    | Diabase lava   | 0.012-9.104   |  8.27 | 17.9  | 0.145
SL41    | Amph. syenite  | 0.025-1.362   |  2.47 |  7.13 | 0.023
SL44    | Neph. syenite  | 0.148-3.857   |  3.57 | 10.8  | 0.034
SL45    | Neph. syenite  | 0.184-1.59    |  4.89 | 16.11 | 0.053
SL46    | Neph. syenite  | 0.081-9.537   |  1.38 |  4.91 | 0.013
CLD10   | Q-tz syenite   | 0.012-0.412   | 10.43 | 36.4  | 0.096
GLG2    | Gabbro         | 4.291-18.15   |  7.32 | 17.26 | 0.061
DM1     | Amph. syenite  | 2.98-96.63    |  7.55 | 18.43 | 0.1045
SL23    | Q-tz syenite   | 0.038-0.261   |  5.58 | 26.44 | 0.059
SL30    | Amph. syenite  | 0.118-1.012   |  7.57 | 16.94 | 0.108
SL31    | Diabase Lava   | 0.011-0.632   | 10.86 | 18.6  | 0.156
SL40    | Q-tz syenite   | 0.031-0.781   |  3.72 |  9.51 | 0.039
"""

# lithology as printed -> geologic_classes, geologic_types, lithologies, phrase.
# The paper does not say how the rhyolite (CLD1) and basalt (SL33) occur.
ROCK = {
    'E. gabbro': ('Intrusive', 'Pluton', 'Gabbro', 'eastern gabbro'),
    'W. gabbro': ('Intrusive', 'Pluton', 'Gabbro', 'western gabbro'),
    'Gabbro': ('Intrusive', 'Pluton', 'Gabbro', 'gabbro'),
    'Aug. syenite': ('Intrusive', 'Pluton', 'Syenite', '(ferro)augite syenite'),
    'Augite syenite': ('Intrusive', 'Pluton', 'Syenite', 'augite syenite'),
    'Amph. syenite': ('Intrusive', 'Pluton', 'Syenite', 'amphibole syenite'),
    'Q-tz syenite': ('Intrusive', 'Pluton', 'Quartz Syenite', 'quartz syenite'),
    'Neph. syenite': ('Intrusive', 'Pluton', 'Nepheline syenite', 'nepheline syenite'),
    'Neph.syenite': ('Intrusive', 'Pluton', 'Nepheline syenite', 'nepheline syenite'),
    'Mafic dyke': ('Intrusive', 'Volcanic Dike', 'Mafic Dike', 'mafic dyke'),
    'Pillowed lava': ('Extrusive', 'Volcanic Pillow', 'Basalt',
                      'pillowed lava (pillow basalt, section 2)'),
    'Diabase lava': ('Extrusive', 'Lava Flow', 'Diabase', 'diabase lava'),
    'Diabase Lava': ('Extrusive', 'Lava Flow', 'Diabase', 'diabase lava'),
    'Rhyolite': ('Igneous', 'Not Specified', 'Rhyolite', 'rhyolite'),
    'Basalt': ('Igneous', 'Not Specified', 'Basalt', 'basalt'),
}
METHOD = {'AF': 'LP-DIR-AF', 'T': 'LP-DIR-T', 'AF+T': 'LP-DIR-AF:LP-DIR-T'}
METHOD_TEXT = {
    'AF': 'alternating-field demagnetization only',
    'T': 'thermal demagnetization only',
    'AF+T': ('alternating-field and thermal demagnetization of sister specimens of each '
             'core, the two ChRM directions averaged per sample'),
}
CENTER_TEXT = {
    'A': ('Center A (eastern part of the complex: eastern gabbro and ferroaugite syenite; '
          'intrusive-center naming of Lewchuk and Symons, 1990)'),
    'B': ('Center B (western part of the complex: western gabbros, syenites and minor '
          'pillowed lavas; intrusive-center naming of Lewchuk and Symons, 1990)'),
    'C': ('Center C (central part of the complex: biotite gabbros, nepheline and amphibole '
          'syenites, quartz syenites and lavas; intrusive-center naming of Lewchuk and '
          'Symons, 1990)'),
}
UNUSED = {  # Table 2 footnote b and section 3.3
    'CLD2': 'all measured specimens gave erratic demagnetization paths (0 of 13 samples)',
    'SL8': 'all measured specimens gave erratic demagnetization paths (0 of 7 samples)',
    'SL26': 'all measured specimens gave erratic demagnetization paths (0 of 12 samples)',
    'GLG2': 'all measured specimens gave erratic demagnetization paths (0 of 6 samples)',
    'SL30': 'all measured specimens gave erratic demagnetization paths (0 of 7 samples)',
    'SL40': 'all measured specimens gave erratic demagnetization paths (0 of 8 samples)',
    'SL23': 'too few samples with a ChRM (N = 2 of 8 < 3); no direction is tabulated',
    'SL31': 'too few samples with a ChRM (N = 2 of 8 < 3); no direction is tabulated',
    'CLD9': ('the site mean (D 351.6, I 63.4) is interpreted by the authors as a '
             'remagnetization by the present-day geomagnetic field'),
}
for _s in ('NMG2', 'NMG3', 'MD1', 'MD2', 'MD3', 'MD4', 'SL48', 'CLD10', 'DM1'):
    UNUSED[_s] = 'scattered sample ChRM directions (alpha95 > 15)'

# ---- ages: ages of the rocks (Heaman and Machado, 1992; 2 sigma) ----------------------
CENTER_AGE = {
    'A': dict(age=1108, age_sigma=0.5),
    'B': dict(age=1107, age_low=1102, age_high=1116),
    'C': dict(age=1109, age_low=1105, age_high=1117),
}
CENTER_AGE_TEXT = {
    'A': ('Age of the rock: U-Pb 1108 +/- 1 Ma (2 sigma; age_sigma is 1 sigma), the '
          'upper-intercept age of the Two Duck Lake olivine gabbro of the eastern part of the '
          'complex (Heaman and Machado, 1992), which the paper assigns to Center A (Table 4).'),
    'B': ('Age of the rock: U-Pb zircon upper-intercept age 1107 +9/-5 Ma (2 sigma) of the '
          'granite of the western part of the complex (Heaman and Machado, 1992), given as '
          'age_low/age_high 1102-1116 Ma. The paper infers from the normal polarity of Center '
          'B that it was emplaced a few million years after Centers A and C, at about '
          '1100-1105 Ma (section 4.3).'),
    'C': ('Age of the rock: U-Pb upper-intercept age 1109 +8/-4 Ma (2 sigma) of the nepheline '
          'syenite of the central part of the complex (Heaman and Machado, 1992), given as '
          'age_low/age_high 1105-1117 Ma.'),
}


def read_table2(text):
    cols = ['site', 'center', 'lat', 'lon', 'n_used', 'n_measured', 'method', 'dec', 'inc',
            'a95', 'k', 'vgp_lat', 'vgp_lon']
    df = pd.DataFrame([line.split() for line in text.strip().splitlines()], columns=cols)
    for c in cols[2:]:
        if c != 'method':
            df[c] = pd.to_numeric(df[c].replace('-', None))
    return df


def angle(lat1, lon1, lat2, lon2):
    return pmag.angle([lon1, lat1], [lon2, lat2])[0]


table2 = read_table2(TABLE_2)
table1 = pd.DataFrame([[x.strip() for x in line.split('|')]
                       for line in TABLE_1.strip().splitlines()],
                      columns=['site', 'lith', 'nrm', 'hc', 'hcr', 'mrs_ms']).set_index('site')
assert set(table1.index) == set(table2['site'])
assert len(table2) == 58  # '504 core samples from 58 sites' (section 3.1)

CLD3_PRINTED_LON = -86.673
table2.loc[table2['site'] == 'CLD3', 'lon'] = -86.33

# transcription check: each tabulated VGP against the VGP of its direction (reversed
# sites against the antipode of the printed, normal-polarity VGP)
table2['reversed'] = table2['inc'] < 0
table2['vgp_lat_site'] = table2['vgp_lat'].where(~table2['reversed'], -table2['vgp_lat'])
table2['vgp_lon_site'] = table2['vgp_lon'].where(~table2['reversed'],
                                                 (table2['vgp_lon'] + 180) % 360)
# CLD1: the VGP listed in Table 2 follows from the coordinates given for the site in the
# first author's 2014 dissertation; its VGP is calculated there so that the Center A,
# A+C and A+B+C poles are reproduced at their printed precision.
VGP_COORDS = {'CLD1': (48.769, -86.547)}
VGP_CALC = {}
for r in table2[table2['vgp_lat'].notna()].itertuples():
    vlat, vlon = VGP_COORDS.get(r.site, (r.lat, r.lon))
    plon, plat, dp, dm = pmag.dia_vgp(r.dec, r.inc, r.a95, vlat, vlon % 360)
    VGP_CALC[r.site] = (round(plat, 2), round(plon, 2), round(dp, 1), round(dm, 1),
                        angle(plat, plon, r.vgp_lat_site, r.vgp_lon_site))
mismatch = {s for s, v in VGP_CALC.items() if v[4] > 0.3}
assert mismatch == {'CLD11', 'GLG', 'CLD8'}, mismatch

# ---- sites ------------------------------------------------------------------------------
sites = []
for r in table2.itertuples():
    rock = table1.loc[r.site]
    gclass, gtype, lith, rock_text = ROCK[rock['lith']]
    used = r.site not in UNUSED
    text = f"{CENTER_TEXT[r.center]}; {rock_text} (Table 1 lithology '{rock['lith']}'). "
    if r.site == 'SL13E':
        text += ('Thin mafic dyke cutting syenite in the central part of the complex (Table 2 '
                 'footnote c); the authors include its normal-polarity direction in the '
                 'Center B mean (section 3.3) and suggest it may be site 13 of Lewchuk and '
                 'Symons (1990) (section 4.2). ')
    if r.site == 'SL13':
        text += ('Same location as site 13 of Lewchuk and Symons (1990), which gave a '
                 'normal-polarity direction (D 300, I 57); SL13 gives an almost antipodal, '
                 'reversed direction (section 4.2). ')
    if r.site == 'CLD3':
        text += (f'Site longitude read from the paper\'s map (Fig. 2), which plots CLD3 on '
                 f'Highway 17 between SL34/CLD4 and CLD1 at the eastern margin of the complex: '
                 f'86.331W from the map graticule, ~86.325W after allowing for the 0.005-0.007 '
                 f'deg offset of the neighbouring sites CLD1, CLD2, CLD5-6 and CLD6 from their '
                 f'listed positions; entered as -86.33 (geographic_precision 0.01 deg). Table 2 '
                 f'lists {CLD3_PRINTED_LON}, which would place the site ~25 km west of the other '
                 f'Center A sites; the map position reproduces the listed VGP (57.5N, 224.3E) '
                 f'within 0.03 deg. Latitude as listed in Table 2. ')
    if r.site == 'NMG1':
        text += ('Coordinates as listed in Table 2 (48.847N, -86.484), which reproduce the '
                 'listed VGP; Fig. 2 plots NMG1 at the eastern margin of the complex next to '
                 'NMG2 and NMG3, at about 48.80N, 86.32W. ')
    if r.site == 'GLG':
        text += ('Coordinates as listed in Table 2 (48.821N, -86.485); Fig. 2 plots GLG about '
                 '2.5 km north of GLS, at about 48.844N, 86.486W. ')
    if r.site == 'CLD1':
        text += ("The VGP listed in Table 2 (46.2N, 200.6E) follows from the site coordinates "
                 "given for CLD1 in the first author's 2014 Michigan Tech dissertation (48.769N, "
                 "-86.547), and vgp_lat/vgp_lon are calculated at those coordinates, with which "
                 "the authors' Center A pole is reproduced; the coordinates of this row are those "
                 "of the paper (with which Fig. 2 agrees), which give a VGP 0.2 deg away. ")
    if r.site == 'GLS2':
        text += ('Table 1 and section 3.1 name this quartz syenite site north of Neys '
                 'Provincial Park GLS1; Table 2 names it GLS2 (there is no GLS1 in Table 2 and '
                 'no GLS2 in Table 1). The rock magnetic values below are those of GLS1 in '
                 'Table 1. ')
    if r.site == 'CLD11':
        text += ('The VGP listed in Table 2 is 0.4 deg from the VGP calculated from the '
                 'tabulated direction and coordinates. ')
    if pd.notna(r.dec):
        text += (f'Site mean of the characteristic remanence (ChRM: principal component lines '
                 f'through the origin, >= 5 steps, MAD < 10, unblocked by 575-590 C or 60-70 '
                 f'mT) from {METHOD_TEXT[r.method]}; {int(r.n_used)} of the '
                 f'{int(r.n_measured)} measured samples used (Table 2). ')
    else:
        text += f'No site mean: {UNUSED[r.site]}; {int(r.n_measured)} samples measured (Table 2). '
    if not used and pd.notna(r.dec):
        text += f'Not used by the authors (Table 2 footnote b): {UNUSED[r.site]}. '
    text += (f"Rock magnetism (Table 1): NRM {rock['nrm']} A/m, Hc {rock['hc']} mT, Hcr "
             f"{rock['hcr']} mT, Mrs/Ms {rock['mrs_ms']}. ")
    vgp = {}
    if used and pd.notna(r.vgp_lat):
        plat, plon, dp, dm, _ = VGP_CALC[r.site]
        if r.site != 'CLD1':
            text += (f'vgp_lat/vgp_lon are calculated from the site direction and coordinates; '
                     f'Table 2 lists the VGP as {r.vgp_lat}N, {r.vgp_lon}E (in normal polarity). ')
        if r.site in ('GLG', 'CLD8'):
            text += ('The authors\' Center C and A+C poles are reproduced with the calculated '
                     'VGP. ')
        vgp = dict(vgp_lat=plat, vgp_lon=plon, vgp_dp=dp, vgp_dm=dm)
    if r.site == 'CLD9':
        text += ('The direction is a present-day-field overprint, not a magnetization of '
                 'Coldwell age; the age given is that of the rock. ')
    text += CENTER_AGE_TEXT[r.center]
    row = dict(site=r.site, location=LOCATION, result_type='i',
               result_quality='g' if used else 'b', citations='This study',
               geologic_classes=gclass, geologic_types=gtype, lithologies=lith,
               method_codes=METHOD[r.method] + (':DE-BFL:DE-FM' if pd.notna(r.dec) else ''),
               lat=r.lat, lon=round(r.lon % 360, 3), age_unit='Ma', **CENTER_AGE[r.center],
               dir_n_total_samples=int(r.n_measured), description=text.strip())
    if r.site == 'GLS2':
        row['site_alternatives'] = 'GLS1'
    if pd.notna(r.dec):
        row.update(dir_tilt_correction=0, dir_dec=r.dec, dir_inc=r.inc, dir_alpha95=r.a95,
                   dir_k=r.k, dir_n_samples=int(r.n_used),
                   dir_comp_name='present-day field overprint' if r.site == 'CLD9'
                   else 'characteristic')
        if used:
            row.update(dir_polarity='r' if r.reversed else 'n', **vgp)
    elif r.n_used > 0:
        row.update(dir_n_samples=int(r.n_used))
    if r.site == 'CLD3':
        row['geographic_precision'] = 0.01
    sites.append(row)
sites = pd.DataFrame(sites)
for c in ('dir_n_samples', 'dir_n_total_samples', 'dir_tilt_correction'):
    sites[c] = sites[c].astype('Int64')
good = sites[sites['result_quality'] == 'g']
assert len(good) == 40 and good['dir_n_samples'].sum() == 238  # abstract, section 3.3

# ---- statistics ---------------------------------------------------------------------------
center = dict(zip(table2['site'], table2['center']))


def pole_of(sub):
    """Fisher mean of polarity-unified site VGPs, reported in the northern hemisphere."""
    blk = [[lon, lat] for lat, lon in zip(sub['vgp_lat'], sub['vgp_lon'])]
    blk = [[(lon + 180) % 360, -lat] if lat < 0 else [lon, lat] for lon, lat in blk]
    p = pmag.fisher_mean(blk)
    if p['inc'] < 0:
        p['dec'], p['inc'] = (p['dec'] + 180) % 360, -p['inc']
    return p


def direction_of(sub):
    """Fisher mean of site directions in the observed majority polarity."""
    up = (sub['dir_inc'] < 0).sum() >= len(sub) / 2
    blk = [[(d + 180) % 360, -i] if (i < 0) != up else [d, i]
           for d, i in zip(sub['dir_dec'], sub['dir_inc'])]
    return pmag.fisher_mean(blk)


def reversal(rev, nor):
    """Angle between means and simulated critical angle (McFadden & McElhinny, 1990)."""
    a = [[d, i] for d, i in zip(rev['dir_dec'], rev['dir_inc'])]
    b = [[(d + 180) % 360, -i] for d, i in zip(nor['dir_dec'], nor['dir_inc'])]
    with contextlib.redirect_stdout(io.StringIO()):
        out = ipmag.common_mean_watson(a, b, NumSims=10000, plot='no', random_seed=23)
    return out[1], out[2]


grp = {c: good[good['site'].map(center) == c] for c in 'ABC'}
TESTS = {k: reversal(rev, grp['B']) for k, rev in
         [('A', grp['A']), ('C', grp['C']), ('AC', pd.concat([grp['A'], grp['C']]))]}
assert TESTS['C'][0] < TESTS['C'][1] and TESTS['A'][0] > TESTS['A'][1] \
    and TESTS['AC'][0] > TESTS['AC'][1], TESTS
REVERSAL_TEXT = (
    'Reversal test (McFadden and McElhinny, 1990), recalculated from the Table 2 site means: '
    'the reversed Center C and the normal Center B site means share a common mean (angle '
    f'between the means {TESTS["C"][0]:.1f} deg, critical angle {TESTS["C"][1]:.1f} deg; '
    'classification B), consistent with a symmetrically reversing field during emplacement '
    'of the complex; this test is the basis of the reversal_test entry. For the reversed '
    'Center A means and the combined Center A and C means, the angles between their mean '
    f'and the Center B mean ({TESTS["A"][0]:.1f} and {TESTS["AC"][0]:.1f} deg, as reported '
    'in section 3.3) exceed the recalculated critical angles '
    f'({TESTS["A"][1]:.1f} and {TESTS["AC"][1]:.1f} deg). Critical angles from a Watson V '
    'simulation (10000 draws).')
METHODS_TEXT = (
    'ChRM from alternating-field and/or thermal demagnetization (12-15 steps; 2G 760-R '
    'magnetometer), principal component lines (>= 5 steps, MAD < 10); site means accepted '
    'with N >= 3 samples and alpha95 <= 15. Geographic coordinates: the paper rules out '
    'tilting of the complex (flat-lying Coubran-Geordie lakes lavas; section 4.1).')

GEO = dict(location=LOCATION, location_type='Region', continent_ocean='North America',
           country='Canada', state_province='Ontario',
           region='North shore of Lake Superior near Marathon, Midcontinent Rift',
           terranes='Laurentia:Superior Province', tectonic_settings='Continental Rift',
           age_unit='Ma')
POLES = [
    dict(key='A', name='Pole CCA (Center A)', members='A', test=True,
         published=(48.9, 209.8, 6.0, 45, 115.1, -66.1, 4.4, 82), age=CENTER_AGE['A'],
         sb='13.2 +/- 5.2',
         text=('Mean VGP of the 14 accepted sites of Center A (eastern gabbro and ferroaugite '
               'syenite; 19 sampled), pole CCA (section 3.3, Tables 2 and 4). Table 3 lists '
               'k 89 for Center A and Table 2 lists k 82, which the site means give. ' + CENTER_AGE_TEXT['A'])),
    dict(key='B', name='Pole CCn (Center B)', members='B', test=True,
         published=(44.9, 193.2, 8.0, 37, 298.0, 56.9, 5.8, 70), age=CENTER_AGE['B'],
         sb='12.1 +/- 5.5',
         text=('Mean VGP of the 10 accepted normal-polarity sites of Center B (western part: 9 '
               'of 15 sampled sites, plus the mafic dyke SL13E of the central part), pole CCn '
               '(section 3.3, Tables 2 and 4). The authors note that 10 sites may not fully '
               'average secular variation (sections 4.1, 4.3). ' + CENTER_AGE_TEXT['B'])),
    dict(key='C', name='Pole CCC (Center C)', members='C', test=True,
         published=(45.6, 203.7, 7.6, 24, 114.6, -61.6, 5.8, 42), age=CENTER_AGE['C'],
         sb='15.9 +/- 3.6',
         text=('Mean VGP of the 16 accepted sites of Center C (central part; 23 sampled), pole '
               'CCC (section 3.3, Tables 2 and 4). Table 3 lists the Center C mean as D 114.5, '
               'k 48 and Table 2 as D 114.6, k 42, which the site means give. Table 4 gives no age for this '
               'pole. ' + CENTER_AGE_TEXT['C'])),
    dict(key='AC', name='Pole CCr (Centers A and C)', members='AC', test=True,
         published=(47.2, 206.5, 4.8, 31, 114.8, -63.7, 3.6, 54),
         age=dict(age=1107, age_low=1105, age_high=1109),
         gpmdb=GPMDB_RESULT, pole_vv_q=6,
         text=('The paper\'s preferred pole for the Coldwell Complex (pole CCr; abstract, '
               'section 3.3, Table 2): Fisher mean of the VGPs of the 30 reversed-polarity '
               'sites of Centers A (eastern part) and C (central part), whose group mean '
               'directions the authors find statistically indistinguishable. The abstract, the '
               'text and Table 2 give the pole longitude as 206.5E, which the site directions '
               'reproduce (Table 4 lists 206.7E). Age: 1107 Ma, the midpoint of 1105-1109 '
               'Ma (age_low/age_high). The pole averages reversed-polarity bodies of Centers A '
               'and C, dated by U-Pb at 1108 +/- 1 Ma (Center A gabbro) and 1109 +8/-4 Ma '
               '(Center C nepheline syenite; 2 sigma; Heaman and Machado, 1992), that record '
               'paleolatitudinal change within the Alona Bay reversed-polarity chron, so they '
               'need not share one age; Table 4 of the paper assigns the 1108 +/- 1 Ma date '
               'of the Center A gabbro to this pole. Reliability: six of seven points of Van der Voo (1990) '
               '(section 4.3). GPMDB gives B = 40 and N = 238, the numbers of sites and samples '
               'of Centers A, B and C combined (this pole uses 30 sites and 177 samples).')),
    dict(key='ABC', name='Centers A, B and C combined pole', members='ABC', test=True,
         published=(46.8, 203.1, 4.2, None, None, None, None, None),
         age=dict(age=1107, age_low=1105, age_high=1109),
         text=('Fisher mean of the polarity-unified VGPs of all 40 accepted sites of Centers A, '
               'B and C, normal and reversed (Table 4; no mean direction or K is printed). The '
               'authors prefer pole CCr because the 10 Center B sites may not fully average '
               'secular variation (section 4.3). Age: 1107 Ma (1105-1109 Ma) for the pole that '
               'averages the reversed Centers A and C, dated at 1108 +/- 1 Ma and 1109 +8/-4 Ma '
               '(2 sigma), with the normal Center B, which the paper infers to be a few million '
               'years younger; the western granite gives 1107 +9/-5 Ma (Heaman and Machado, '
               '1992).')),
    dict(key='B-CCW1', name='Center B without site CCW1', members='B', drop='CCW1',
         published=(47.2, 190.9, 6.7, 60, 301.7, 56.8, 5.0, 108), age=CENTER_AGE['B'],
         sb='9.1 +/- 4.8',
         text=('Center B mean without site CCW1, whose more westerly declination plots away from '
               'the other sites (section 3.3). Not the authors\' preferred result: they attribute '
               'the CCW1 direction to secular variation and use all 10 sites (pole CCn). '
               + CENTER_AGE_TEXT['B'])),
    dict(key='B-SL13E', name='Center B without site SL13E', members='B', drop='SL13E',
         direction_only=True, published=(None, None, None, None, 297.5, 56.3, 6.4, 66),
         age=CENTER_AGE['B'],
         text=('Center B mean of the 9 sites of the western part, without the mafic dyke SL13E '
               'of the central part (Table 2, mean without SL13E; section 3.3). The paper gives '
               'no pole for this mean. ' + CENTER_AGE_TEXT['B'])),
]


def r1(x):
    return None if x is None else round(float(x), 1)


locations = [dict(
    **GEO, sites=':'.join(sites['site']), lat_s=sites['lat'].min(), lat_n=sites['lat'].max(),
    lon_w=sites['lon'].min(), lon_e=sites['lon'].max(), age_low=1102, age_high=1117,
    geologic_classes=':'.join(sorted(sites['geologic_classes'].unique())),
    lithologies=':'.join(sorted(sites['lithologies'].unique())),
    method_codes='LP-DIR-AF:LP-DIR-T', citations='This study:' + HM1992,
    description=(
        'Coldwell Complex, a ca. 1108 Ma multiphase alkaline intrusion ~25 km across on the '
        'north shore of Lake Superior (Midcontinent Rift), intruding Archean rocks of the '
        'Superior Province. 58 sites (Table 2), grouped by the authors into the intrusive '
        'centers of Lewchuk and Symons (1990): Center A (eastern part; 19 sites, 14 accepted), '
        'Center B (western part; 16 sites incl. the mafic dyke SL13E of the central part, 10 '
        'accepted), Center C (central part; 23 sites, 16 accepted). 18 sites were not used '
        '(result_quality b): 6 with erratic demagnetization, 11 with alpha95 > 15 or N < 3, '
        'and CLD9 (present-field remagnetization). Mostly road-cut sites along Highway 17, plus '
        'sites north of Marathon, near Middleton, north of Neys Provincial Park and north of '
        'the Red Sucker Point Provincial Nature Reserve (section 3.1); site coordinates as '
        'listed in Table 2, except CLD3, whose longitude is read from Fig. 2. Age: age_low/age_high span the U-Pb ages of the three centers (Heaman '
        'and Machado, 1992; 2 sigma): 1108 +/- 1 Ma (Center A), 1109 +8/-4 Ma (Center C) and '
        '1107 +9/-5 Ma (Center B). ' + METHODS_TEXT + ' VGP angular dispersion corrected for '
        'within-site dispersion (section 4.1, jackknife 95% limits): Center A Sb = 13.2 +/- '
        '5.2 (N 14), Center B 12.1 +/- 5.5 (N 10; 9.1 +/- 4.8 without CCW1), Center C 15.9 '
        '+/- 3.6 (N 16).'))]

CHECK = []
for p in POLES:
    sub = good[good['site'].map(center).isin(list(p['members']))]
    if p.get('drop'):
        sub = sub[sub['site'] != p['drop']]
    d = direction_of(sub)
    plat, plon, a95, kk, pdec, pinc, pa95, pk = p['published']
    if pdec is not None:
        assert angle(d['inc'], d['dec'], pinc, pdec) < 0.2, (p['name'], d)
    row = dict(**GEO, result_name=p['name'], result_type='a', result_quality='g',
               sites=':'.join(sub['site']), lat_s=sub['lat'].min(), lat_n=sub['lat'].max(),
               lon_w=sub['lon'].min(), lon_e=sub['lon'].max(), **p['age'],
               geologic_classes=':'.join(sorted(sub['geologic_classes'].unique())),
               lithologies=':'.join(sorted(sub['lithologies'].unique())),
               method_codes='LP-DIR-AF:LP-DIR-T:DE-BFL:DE-FM'
               + ('' if p.get('direction_only') else ':DE-VGP'),
               citations='This study:' + HM1992,
               dir_tilt_correction=0, dir_dec=r1(d['dec']), dir_inc=r1(d['inc']), dir_alpha95=r1(d['alpha95']),
               dir_k=r1(d['k']), dir_n_sites=int(d['n']))
    text = p['text']
    if not p.get('direction_only'):
        q = pole_of(sub)
        sep = angle(q['inc'], q['dec'], plat, plon)
        assert (r1(q['inc']), r1(q['dec']), r1(q['alpha95'])) == (plat, plon, a95), \
            (p['name'], q)
        n_rev = int((sub['dir_inc'] < 0).sum())
        row.update(pole_comp_name='characteristic', pole_lat=r1(q['inc']),
                   pole_lon=r1(q['dec']), pole_alpha95=r1(q['alpha95']), pole_k=r1(q['k']),
                   pole_n_sites=int(q['n']),
                   pole_reversed_perc=round(100 * n_rev / len(sub), 1))
        CHECK.append((p['name'], q['inc'], q['dec'], q['alpha95'], q['k'], int(q['n']),
                      plat, plon, a95, sep))
        text += (' Pole: Fisher mean of the polarity-unified site VGPs, reported in the '
                 'northern hemisphere as in the paper.')
    if p.get('test'):
        row['reversal_test'] = 'Rb'
        text += ' ' + REVERSAL_TEXT
    if p.get('sb'):
        text += (' VGP angular dispersion corrected for within-site dispersion (section 4.1): '
                 f'Sb = {p["sb"]} (jackknife 95% limits).')
    if p.get('gpmdb'):
        row['external_database_ids'] = f'GPMDB[{p["gpmdb"]}]'
    if p.get('pole_vv_q'):
        row['pole_vv_q'] = p['pole_vv_q']
    row['description'] = text
    locations.append(row)
locations = pd.DataFrame(locations)
for c in ('dir_n_sites', 'pole_n_sites', 'pole_vv_q', 'dir_tilt_correction'):
    locations[c] = locations[c].astype('Int64')

# ---- ages table: the U-Pb dates (Heaman and Machado, 1992; 2 sigma) ------------------
AGE_ROWS = [
    (dict(age=1108, age_sigma=0.5),
     'U-Pb baddeleyite and zircon upper-intercept age, Two Duck Lake olivine gabbro (sample '
     '86CL1), eastern part of the complex (Center 1 of Heaman and Machado; Center A of the '
     'paper), 1108 +/- 1 Ma (2 sigma; age_sigma is 1 sigma); taken by Heaman and Machado as '
     'the best estimate of the emplacement age of the complex.'),
    (dict(age=1107, age_low=1106, age_high=1112),
     'U-Pb baddeleyite and zircon upper-intercept age, layered gabbro (sample CL7), eastern '
     'part of the complex (Center 1 of Heaman and Machado; Center A of the paper), 1107 +5/-1 '
     'Ma (2 sigma) in their text; their concordia diagram (Fig. 3A) is labelled 1107 +2/-1 Ma.'),
    (dict(age=1109, age_low=1105, age_high=1117),
     'U-Pb zircon and sphene upper-intercept age, nepheline syenite (sample CL22), central '
     'part of the complex (Center 2 of Heaman and Machado; Center C of the paper), 1109 +8/-4 '
     'Ma (2 sigma).'),
    (dict(age=1107, age_low=1102, age_high=1116),
     'U-Pb zircon upper-intercept age, granite (sample CL2), western margin of the complex '
     '(Center 3 of Heaman and Machado; the western part is Center B of the paper), 1107 +9/-5 '
     'Ma (2 sigma; one of four fractions omitted from the regression). Granites were not '
     'sampled for paleomagnetism.'),
]
ages = pd.DataFrame([dict(location=LOCATION, **a, age_unit='Ma',
                          method_codes='GM-UPB:GM-UPB-CC-T0', citations=HM1992,
                          timescale_eon='Proterozoic', timescale_era='Mesoproterozoic',
                          timescale_period='Stenian', description=d) for a, d in AGE_ROWS])


# ---- write ----------------------------------------------------------------------------
def write_magic(df, kind, path):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False)


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locations, 'locations', os.path.join(OUT, 'locations.txt'))
write_magic(ages, 'ages', os.path.join(OUT, 'ages.txt'))
stamp = date.today().strftime('%d.%b.%Y')
combined = os.path.join(OUT, f'Kulakov2014_Coldwell_Complex_{stamp}.txt')
with open(combined, 'w', encoding='utf-8') as f:
    for i, (kind, df) in enumerate([('locations', locations), ('sites', sites),
                                    ('ages', ages)]):
        if i:
            f.write('>>>>>>>>>>\n')
        f.write(f'tab delimited\t{kind}\n')
        df.to_csv(f, sep='\t', index=False)

print('reversal tests (angle, critical):', {k: (round(a, 2), round(c, 2))
                                            for k, (a, c) in TESTS.items()})
for name, lat, lon, a95, kk, n, plat, plon, pa95, sep in CHECK:
    print(f'  {name:36s} {lat:5.2f}N {lon:6.2f}E A95 {a95:4.2f} K {kk:5.1f} N {n:2d}'
          f'  | published {plat}/{plon} A95 {pa95}  ({sep:.2f} deg)')
print(f'wrote {combined}')

from validate_magic_contribution import validate_upload_file  # noqa: E402
if not validate_upload_file(combined, tables=['locations', 'sites', 'ages']):
    raise SystemExit('validation failed')
