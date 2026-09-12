# Canonical prose values (do not hand-type; regenerate with results_pipeline.py)

- Combo: 940,510 cost / 12.18 QALY / 16.21 discounted LY
- Pembro: 744,153 cost / 7.55 QALY / 10.61 discounted LY
- Incremental: 196,357 cost / 4.62 QALY / 5.60 discounted LY / 9.68 undiscounted LY
- ICER: 42,499 per QALY
- NMB @100K: 265,669   @150K: 496,682
- RF/LR/DM QALY (combo): 8.98 / 1.66 / 1.54
- Incremental RF/LR/DM QALY: 4.41 / 0.36 / -0.15
- Administration cost per arm: 8,869

## Analytic threshold prices
- 100,000/QALY -> 465,668.83 per course
- 150,000/QALY -> 696,681.77 per course

## Scenario ICERs
- Base case: 42,499 (dCost 196,357, dQALY 4.62)
- Weibull OS: 56,945 (dCost 365,666, dQALY 6.42)
- Rank #1 joint model: 61,165 (dCost 390,941, dQALY 6.39)
- Unconstrained (no GP): 51,753 (dCost 373,488, dQALY 7.22)
- GP floor 48 mo: 42,483 (dCost 195,854, dQALY 4.61)
- GP floor 60 mo (base): 42,499 (dCost 196,357, dQALY 4.62)
- GP floor 72 mo: 42,551 (dCost 197,730, dQALY 4.65)
- GP floor 96 mo: 42,765 (dCost 202,571, dQALY 4.74)
- 10-year horizon: 74,047 (dCost 96,740, dQALY 1.31)
- 20-year horizon: 40,664 (dCost 129,095, dQALY 3.17)
- 30-year horizon: 41,265 (dCost 174,890, dQALY 4.24)
- 40-year horizon (base): 42,499 (dCost 196,357, dQALY 4.62)
- Waning V1 (OS-survival convergence): 18,765 (dCost 37,098, dQALY 1.98)
- Waning V2 (hazard convergence): 104,737 (dCost 239,400, dQALY 2.29)
- No direct OS benefit: -20,616 (dCost -14,967, dQALY 0.73)
- 0% discount: 37,597 (dCost 292,108, dQALY 7.77)
- 5% discount: 49,183 (dCost 167,066, dQALY 3.40)
