# Canonical prose values (do not hand-type; regenerate with results_pipeline.py)

- Combo: 948,604 cost / 12.24 QALY / 16.30 discounted LY
- Pembro: 746,904 cost / 7.59 QALY / 10.66 discounted LY
- Incremental: 201,699 cost / 4.65 QALY / 5.64 discounted LY / 9.68 undiscounted LY
- ICER: 43,346 per QALY
- NMB @100K: 263,629   @150K: 496,293
- RF/LR/DM QALY (combo): 9.02 / 1.67 / 1.56
- Incremental RF/LR/DM QALY: 4.44 / 0.36 / -0.15
- Administration cost per arm: 13,340

## Analytic threshold prices
- 100,000/QALY -> 463,628.73 per course
- 150,000/QALY -> 696,292.77 per course

## Scenario ICERs
- Base case: 43,346 (dCost 201,699, dQALY 4.65)
- Weibull OS: 57,582 (dCost 372,598, dQALY 6.47)
- Rank #1 joint model: 61,803 (dCost 398,050, dQALY 6.44)
- Unconstrained (no GP): 52,337 (dCost 381,133, dQALY 7.28)
- GP floor 48 mo: 43,331 (dCost 201,192, dQALY 4.64)
- GP floor 60 mo (base): 43,346 (dCost 201,699, dQALY 4.65)
- GP floor 72 mo: 43,393 (dCost 203,083, dQALY 4.68)
- GP floor 96 mo: 43,592 (dCost 207,965, dQALY 4.77)
- 10-year horizon: 77,108 (dCost 101,019, dQALY 1.31)
- 20-year horizon: 41,879 (dCost 133,613, dQALY 3.19)
- 30-year horizon: 42,177 (dCost 179,907, dQALY 4.27)
- 40-year horizon (base): 43,346 (dCost 201,699, dQALY 4.65)
- Waning V1 (OS-survival convergence): 20,580 (dCost 40,859, dQALY 1.99)
- Waning V2 (hazard convergence): 106,480 (dCost 244,734, dQALY 2.30)
- No direct OS benefit: -15,725 (dCost -11,458, dQALY 0.73)
- 0% discount: 38,177 (dCost 296,608, dQALY 7.77)
- 5% discount: 49,977 (dCost 172,811, dQALY 3.46)
