# Canonical prose values (do not hand-type; regenerate with results_pipeline.py)

- Combo: 881,300 cost / 11.46 QALY / 15.20 discounted LY
- Pembro: 719,353 cost / 7.41 QALY / 10.38 discounted LY
- Incremental: 161,947 cost / 4.05 QALY / 4.82 discounted LY / 7.79 undiscounted LY
- ICER: 39,976 per QALY
- NMB @100K: 243,167   @150K: 445,724
- RF/LR/DM QALY (combo): 8.55 / 1.53 / 1.37
- Incremental RF/LR/DM QALY: 4.03 / 0.29 / -0.27
- Administration cost per arm: 2,303

## Analytic threshold prices
- 100,000/QALY -> 443,167.15 per course
- 150,000/QALY -> 645,724.06 per course

## Scenario ICERs
- Base case: 39,976 (dCost 161,947, dQALY 4.05)
- Weibull OS: 54,745 (dCost 308,820, dQALY 5.64)
- Rank #1 joint model: 59,624 (dCost 334,543, dQALY 5.61)
- Unconstrained (no GP): 51,741 (dCost 376,791, dQALY 7.28)
- GP floor 48 mo: 39,937 (dCost 161,018, dQALY 4.03)
- GP floor 60 mo (base): 39,976 (dCost 161,947, dQALY 4.05)
- GP floor 72 mo: 40,055 (dCost 163,614, dQALY 4.08)
- GP floor 96 mo: 40,334 (dCost 168,709, dQALY 4.18)
- 10-year horizon: 73,937 (dCost 95,047, dQALY 1.29)
- 20-year horizon: 39,848 (dCost 122,804, dQALY 3.08)
- 30-year horizon: 39,587 (dCost 154,949, dQALY 3.91)
- 40-year horizon (base): 39,976 (dCost 161,947, dQALY 4.05)
- Waning V1 (OS-survival convergence): 20,398 (dCost 39,319, dQALY 1.93)
- Waning V2 (hazard convergence): 104,971 (dCost 226,421, dQALY 2.16)
- No additional OS hazard benefit beyond 5 years: 21,491 (dCost 37,783, dQALY 1.76)
- 0% discount: 33,050 (dCost 211,207, dQALY 6.39)
- 5% discount: 47,535 (dCost 147,152, dQALY 3.10)
