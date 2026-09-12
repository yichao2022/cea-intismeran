# Canonical prose values (do not hand-type; regenerate with results_pipeline.py)

- Combo: 944,979 cost / 12.18 QALY / 16.21 discounted LY
- Pembro: 744,153 cost / 7.55 QALY / 10.61 discounted LY
- Incremental: 200,826 cost / 4.62 QALY / 5.60 discounted LY / 9.68 undiscounted LY
- ICER: 43,466 per QALY
- NMB @100K: 261,200   @150K: 492,213
- RF/LR/DM QALY (combo): 8.98 / 1.66 / 1.54
- Incremental RF/LR/DM QALY: 4.41 / 0.36 / -0.15
- Administration cost per arm: 13,338

## Analytic threshold prices
- 100,000/QALY -> 461,199.83 per course
- 150,000/QALY -> 692,212.77 per course

## Scenario ICERs
- Base case: 43,466 (dCost 200,826, dQALY 4.62)
- Weibull OS: 57,641 (dCost 370,135, dQALY 6.42)
- Rank #1 joint model: 61,864 (dCost 395,410, dQALY 6.39)
- Unconstrained (no GP): 52,372 (dCost 377,957, dQALY 7.22)
- GP floor 48 mo: 43,452 (dCost 200,323, dQALY 4.61)
- GP floor 60 mo (base): 43,466 (dCost 200,826, dQALY 4.62)
- GP floor 72 mo: 43,513 (dCost 202,199, dQALY 4.65)
- GP floor 96 mo: 43,709 (dCost 207,040, dQALY 4.74)
- 10-year horizon: 77,468 (dCost 101,209, dQALY 1.31)
- 20-year horizon: 42,072 (dCost 133,564, dQALY 3.17)
- 30-year horizon: 42,319 (dCost 179,359, dQALY 4.24)
- 40-year horizon (base): 43,466 (dCost 200,826, dQALY 4.62)
- Waning V1 (OS-survival convergence): 21,025 (dCost 41,567, dQALY 1.98)
- Waning V2 (hazard convergence): 106,692 (dCost 243,869, dQALY 2.29)
- No direct OS benefit: -14,460 (dCost -10,498, dQALY 0.73)
- 0% discount: 38,177 (dCost 296,608, dQALY 7.77)
- 5% discount: 50,493 (dCost 171,515, dQALY 3.40)
