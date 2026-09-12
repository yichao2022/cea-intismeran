# Canonical prose values (do not hand-type; regenerate with results_pipeline.py)

- Combo: 937,504 cost / 12.24 QALY / 16.30 discounted LY
- Pembro: 740,150 cost / 7.59 QALY / 10.66 discounted LY
- Incremental: 197,354 cost / 4.65 QALY / 5.64 discounted LY / 9.68 undiscounted LY
- ICER: 42,412 per QALY
- NMB @100K: 267,974   @150K: 500,638
- RF/LR/DM QALY (combo): 9.02 / 1.67 / 1.56
- Incremental RF/LR/DM QALY: 4.44 / 0.36 / -0.15
- Administration cost per arm: 2,241

## Analytic threshold prices
- 100,000/QALY -> 467,973.84 per course
- 150,000/QALY -> 700,637.88 per course

## Scenario ICERs
- Base case: 42,412 (dCost 197,354, dQALY 4.65)
- Weibull OS: 56,911 (dCost 368,253, dQALY 6.47)
- Rank #1 joint model: 61,128 (dCost 393,705, dQALY 6.44)
- Unconstrained (no GP): 51,740 (dCost 376,788, dQALY 7.28)
- GP floor 48 mo: 42,395 (dCost 196,847, dQALY 4.64)
- GP floor 60 mo (base): 42,412 (dCost 197,354, dQALY 4.65)
- GP floor 72 mo: 42,464 (dCost 198,738, dQALY 4.68)
- GP floor 96 mo: 42,681 (dCost 203,620, dQALY 4.77)
- 10-year horizon: 73,791 (dCost 96,674, dQALY 1.31)
- 20-year horizon: 40,517 (dCost 129,268, dQALY 3.19)
- 30-year horizon: 41,159 (dCost 175,562, dQALY 4.27)
- 40-year horizon (base): 42,412 (dCost 197,354, dQALY 4.65)
- Waning V1 (OS-survival convergence): 18,391 (dCost 36,514, dQALY 1.99)
- Waning V2 (hazard convergence): 104,590 (dCost 240,389, dQALY 2.30)
- No additional OS hazard benefit beyond 5 years: 20,824 (dCost 38,751, dQALY 1.86)
- 0% discount: 37,614 (dCost 292,233, dQALY 7.77)
- 5% discount: 48,726 (dCost 168,485, dQALY 3.46)
