# Canonical prose values (do not hand-type; regenerate with results_pipeline.py)

- Combo: 937,567 cost / 12.24 QALY / 16.30 discounted LY
- Pembro: 740,209 cost / 7.59 QALY / 10.66 discounted LY
- Incremental: 197,358 cost / 4.65 QALY / 5.64 discounted LY / 9.68 undiscounted LY
- ICER: 42,413 per QALY
- NMB @100K: 267,970   @150K: 500,634
- RF/LR/DM QALY (combo): 9.02 / 1.67 / 1.56
- Incremental RF/LR/DM QALY: 4.44 / 0.36 / -0.15
- Administration cost per arm: 2,303

## Analytic threshold prices
- 100,000/QALY -> 467,970.35 per course
- 150,000/QALY -> 700,634.39 per course

## Scenario ICERs
- Base case: 42,413 (dCost 197,358, dQALY 4.65)
- Weibull OS: 56,911 (dCost 368,256, dQALY 6.47)
- Rank #1 joint model: 61,129 (dCost 393,709, dQALY 6.44)
- Unconstrained (no GP): 51,741 (dCost 376,791, dQALY 7.28)
- GP floor 48 mo: 42,396 (dCost 196,851, dQALY 4.64)
- GP floor 60 mo (base): 42,413 (dCost 197,358, dQALY 4.65)
- GP floor 72 mo: 42,465 (dCost 198,742, dQALY 4.68)
- GP floor 96 mo: 42,682 (dCost 203,624, dQALY 4.77)
- 10-year horizon: 73,794 (dCost 96,677, dQALY 1.31)
- 20-year horizon: 40,518 (dCost 129,271, dQALY 3.19)
- 30-year horizon: 41,159 (dCost 175,566, dQALY 4.27)
- 40-year horizon (base): 42,413 (dCost 197,358, dQALY 4.65)
- Waning V1 (OS-survival convergence): 18,393 (dCost 36,517, dQALY 1.99)
- Waning V2 (hazard convergence): 104,591 (dCost 240,392, dQALY 2.30)
- No additional OS hazard benefit beyond 5 years: 20,826 (dCost 38,754, dQALY 1.86)
- 0% discount: 37,614 (dCost 292,236, dQALY 7.77)
- 5% discount: 48,727 (dCost 168,489, dQALY 3.46)
