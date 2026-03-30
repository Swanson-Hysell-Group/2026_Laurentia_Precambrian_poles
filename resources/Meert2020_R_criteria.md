# Meert et al. (2020) R-Criteria for Paleomagnetic Pole Reliability

Meert, J.G., Pivarunas, A.F., Evans, D.A.D., Pisarevsky, S.A., Pesonen, L.J., Li, Z.-X., Elming, S.-A., Miller, S.R., Zhang, S., and Salminen, J.M. (2020). The magnificent seven: A proposal for modest revision of the Van der Voo (1990) quality index. *Tectonophysics*, 790, 228549. https://doi.org/10.1016/j.tecto.2020.228549

## Overview

The "R" (Reliability) criteria are a revision of the Van der Voo (1990) Q-criteria. Each pole is scored 0 or 1 on seven criteria (R1-R7), yielding a total R-score out of 7. The R-score is a checklist, not a pass/fail gate; decisions on how to apply it are left to the individual researcher.

## The Seven R-Criteria

### R1: Age of the rock constrained to within +/- 15 Ma; magnetization presumed to be the same age

- Radiometric age constrained to within +/- 15 Ma.
- This is stricter than the original Q1 (which allowed +/- 4% or +/- 40 Ma for the Precambrian).
- For well-defined remagnetizations, the age constraint applies to the age of remagnetization, not the rock.
- A demonstrably synfolding magnetization qualifies if there are independent age constraints on the deformation.

### R2: Techniques and statistical analysis

Three sub-criteria (a-c), all of which should ideally be met:

**(a) Demagnetization:** At least two methods of stepwise demagnetization (e.g., AF and thermal) should be attempted on a pilot suite to demonstrate that individual vector components are being separated effectively.

**(b) Component analysis:** Directional data should be analyzed using PCA (Kirschvink, 1980) or great circle intersections (Halls, 1976, 1978; McFadden and McElhinny, 1988) to separate overlapping unblocking temperature/coercivity components.

**(c) Adequate PSV sampling:** A VGP scatter that adequately averages paleosecular variation assessed by:
- **Deenen et al. (2011) test:** A95 of the *VGP distribution* should fall within N-dependent bounds: `12 * N^(-0.40) <= A95 <= 82 * N^(-0.63)`
- **Statistical thresholds:** N >= 25 (samples), 10 <= K <= 70, B >= 8 sites (minimum 3 samples per site).
- K > 70 warrants suspicion of inadequate PSV averaging.
- K < 10 suggests data are too dispersed.
- A sample is an independently oriented core or block (may consist of one or more specimens).

### R3: Characterization of magnetic mineralogy / rock magnetism

- A reasonable attempt to identify and comment on the magnetic carriers in the study.
- Methods include: rock magnetic tests (IRM acquisition, hysteresis, Day plots, FORC diagrams, 3-axis IRM, low-temperature treatment), magnetic fabric studies (AMS, AIR, AAR), and/or petrographic/microscopic examination (reflected light, SEM, TEM).
- Identification of magnetic carriers aids in determining primary vs. secondary nature of remanence.
- Particularly important for sedimentary redbeds (DRM vs. CRM issues).

### R4: Field tests that constrain the age of magnetization

Any of the following statistically robust field tests:

**(a) Baked contact / inverse baked contact test:**
- Positive (C+): intrusion and baked host have same direction; unbaked host has different direction; no hybrid zone.
- Also positive (C+): baked host matches intrusion, unbaked host is unstable.
- Inconclusive (Co): baked host matches intrusion direction but unbaked host also shows same direction or unstable behavior.
- Negative (C-): all directions are similar, suggesting remagnetization.
- Inverse baked contact tests also qualify.

**(b) Fold/tilt/slump test:**
- Should pass rigorous statistical analyses (McFadden, 1990; McFadden and Jones, 1981; Watson and Enkin, 1983; Tauxe and Watson, 1994; Enkin, 2003).
- Age of folding should be close to the age of the rocks.
- Should be applied stepwise with optimal grouping within 90-110% unfolding.
- Syn-folding magnetizations do not meet R4 unless demonstrably syn-folding remagnetizations.

**(c) Conglomerate test:**
- Positive if clasts show statistically random directions (Watson, 1956; Shipunov et al., 1998; Heslop and Roberts, 2018a).
- N sufficiently large to test randomness (n >= 10 recommended for field practicality).
- Ideal: intraformational conglomerate with clasts from underlying unit.
- Age of conglomerate and relation to bounding units is critical.

**(d) Unconformity test (Kirschvink, 1978):**
- Positive: polarity sequence below unconformity is truncated (discontinuous across unconformable surface).
- Negative: polarity zonation is continuous across unconformity.

### R5: Structural control and tectonic coherence

- Poles from allochthonous or parautochthonous terranes, non-stratified rocks, and regions that have undergone internal vertical axis rotations will not meet R5.
- Results from intrusive rocks younger than the last deformational event may meet this criterion.
- The region must have been a rigid part of the craton since the time the magnetization was acquired.
- For Precambrian: authors should specify their definition of "craton" along with present and paleogeographic bounds.
- Poles based on flattening corrections will not meet R5 unless corroborated by intercalated volcanic rocks or other sedimentary rocks with similar R-value that do not require flattening corrections.

### R6: Presence of magnetic reversals

- Statistically significant antipodal normal and reverse directions.
- Graded using McFadden and McElhinny (1990) reversal test: R_A (gamma_o < 5 deg), R_B (gamma_c < 10 deg), R_C (gamma_c < 20 deg and "Indeterminate"), or R_0 (gamma_c >= 20 deg).
- Also acceptable: support for a common mean via Heslop and Roberts (2018b) test.
- A negative reversal test (gamma_o > gamma_c) does not qualify.
- The test is based on isolated observations from one polarity grouping; R_AI, R_BI, R_CI designations are used.
- A positive reversal test is supportive but not conclusive of primary magnetization (dual-polarity remagnetization is possible).

### R7: No resemblance to younger poles by more than a period

- Paleomagnetic poles with overlapping A95 envelopes with younger poles (of R >= 3) will not meet this criterion.
- Comparison should be made only with poles from stable regions within the connected craton.
- Poles from orogenic belts should not be used for comparison.
- Field tests that constrain the magnetization to be older than the younger pole(s) it resembles will satisfy R7.
- This criterion is contentious for the Precambrian (APWP self-intersection is statistically expected); some authors favor abolishing it.

## Summary Table (Table 3 from paper)

| R | Brief Description | Limits |
|---|---|---|
| 1 | Age constrained; magnetization presumed same age | Radiometric age within +/- 15 Ma |
| 2 | Techniques and statistical analysis | Stepwise demag confirmed by multiple methods; PSV test: N >= 25, 10 <= K <= 70, B >= 8 sites (min 3 samples/site) |
| 3 | Evaluation of remanence carriers | Rock magnetic and/or microscopic examination and identification of magnetic carriers |
| 4 | Field tests that constrain age of magnetization | Fold/tilt test; baked contact test; conglomerate test or other field tests |
| 5 | Structural control and tectonic coherence | Data from thrust sheets or intrusives must be younger than last tectonic deformation; detrital sedimentary rocks that do not require inclination corrections will meet this |
| 6 | Presence of magnetic reversals | Statistically significant antipodal directions: R_A, R_B, or R_C rated (M&M 1990) or support for common mean (H&R 2018b) |
| 7 | No resemblance to younger poles (> period) based on overlapping A95 | Field tests that constrain magnetization to be older than resembled pole(s) |

## Key Equations

**Deenen et al. (2011) A95 envelope for PSV:**
```
A95_min = 12 * N^(-0.40)
A95_max = 82 * N^(-0.63)
```

**S_T (total angular dispersion of VGPs; Cox, 1970):**
```
S_T = sqrt( (1/(n-1)) * sum(delta_i^2) )    for i = 1...n
```
where delta_i = angle between i-th VGP and the mean VGP, n = number of sites (> 5).

**95% deviation angle (approximation for alarm band):**
```
psi_95 = 140 / sqrt(k)
```

## Notes for Implementation

- R2c (Deenen test) is already implemented in `pole_assessment.py` as `Deenen_test()`.
- R2c statistical thresholds (N, K, B) are implemented as `R2_test()`.
- R4 and R6 require case-by-case evaluation from the source publication.
- R5 requires knowledge of the tectonic setting relative to the craton.
- R7 requires comparison against the Laurentia APWP compilation.
