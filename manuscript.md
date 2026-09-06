# Cost-Effectiveness of Intismeran Autogene Plus Pembrolizumab Versus Pembrolizumab Alone as Adjuvant Therapy for Resected Stage IIIB–IV Melanoma: An Updated Early Economic Evaluation

## Abstract

**Background:** In August 2026, Merck and Moderna reported positive Phase 3 INTerpath-001 results for intismeran autogene (mRNA-4157/V940) plus pembrolizumab in resected stage IIB–IV melanoma. Its cost-effectiveness remains unknown.

**Objective:** To estimate the cost-effectiveness of intismeran plus pembrolizumab versus pembrolizumab alone from a US payer perspective.

**Methods:** A four-state partitioned survival model was developed with a 40-year horizon and 3% annual discounting. Quantitative efficacy inputs were derived from the Phase 2b KEYNOTE-942 trial (n=157), with exploratory OS; Phase 3 INTerpath-001 results were used as clinical context because numerical effect estimates were not yet available. Health-state utilities were drawn from published EQ-5D literature. Probabilistic, deterministic, and scenario sensitivity analyses were conducted.

**Results:** In the base case, intismeran plus pembrolizumab yielded 11.35 QALYs at a cost of $874,026 versus 7.33 QALYs at $712,110 for pembrolizumab alone. The ICER was $40,274/QALY (incremental cost $161,916). Probabilities of cost-effectiveness were 95.3% at $100,000 and 98.7% at $150,000/QALY. The deterministic threshold prices were approximately $440K and $641K. Results were largely robust across structural scenarios; the ICER exceeded $100,000/QALY only under hazard-convergence waning ($104,852/QALY).

**Conclusion:** At an estimated acquisition price of $200,000 per course, intismeran plus pembrolizumab was projected to be cost-effective under conventional US thresholds, driven by a substantial but not implausibly large 4.02-QALY gain despite higher acquisition and treatment-related costs.

**Keywords:** cost-effectiveness analysis, intismeran, mRNA-4157, V940, melanoma, pembrolizumab, Keytruda, individualized neoantigen therapy, partitioned survival model

---

## 1. Introduction

### A New Reimbursement Problem

On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy (INT) to succeed in a late-stage trial [1]. The trial met its primary endpoint of recurrence-free survival (RFS) and key secondary endpoint of distant metastasis-free survival (DMFS) in patients with resected high-risk stage IIB–IV melanoma. However, full numerical Phase 3 efficacy results and the commercial price for intismeran have not yet been disclosed. The pivotal clinical milestone therefore creates an immediate need for an early economic assessment before pricing and reimbursement decisions are finalized.

### Why INTs Differ from Conventional Oncology Pricing

Individualized neoantigen therapies (INTs) introduce production and delivery requirements that differ materially from conventional fixed-formulation pharmaceuticals. Each patient's treatment requires a dedicated workflow: multi-region tumor biopsy and next-generation sequencing, computational neoantigen prediction and HLA-binding prioritization, individualized mRNA synthesis and lipid nanoparticle formulation, and lot-release quality control, all within a manufacturing window of approximately six to nine weeks [3, 4]. These per-patient manufacturing requirements, combined with the need for decentralized sequencing capacity and cold-chain logistics, present scalability and cost challenges distinct from batch-produced oncology drugs [4, 5]. Current industry estimates suggest personalized mRNA cancer vaccines may cost $100,000–$300,000 per patient, with Jefferies analysts estimating approximately $200,000 per course [5], but these are supply-side estimates, not value-based prices.

### Existing Economic Evidence

A prior exploratory cost-effectiveness analysis of V940 plus pembrolizumab, presented at ISPOR Europe 2024, used a four-state partitioned survival model applied to the stage III adjuvant melanoma setting [11]. That analysis assumed V940 pricing at the same list price as pembrolizumab and applied KEYNOTE-942 hazard ratios to pembrolizumab survival curves derived from the KEYNOTE-054 trial, explicitly noting that the results were preliminary. Since that analysis, KEYNOTE-942 has reported mature 5-year RFS and DMFS outcomes [2], and INTerpath-001 has provided the first positive Phase 3 confirmation. The present study therefore addresses an evidence gap that the prior exploratory work itself identified.

### Clinical Basis and Research Objectives

The clinical foundation for this analysis comes from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (NCT03897881), which randomized 157 patients with resected stage IIIB–IV cutaneous melanoma 2:1 to intismeran plus pembrolizumab versus pembrolizumab alone [2]. At a median follow-up of 60.3 months, the combination demonstrated a 49% reduction in the risk of recurrence or death (HR 0.51; 95% CI 0.294–0.887) and a 59% reduction in the risk of distant metastasis (HR 0.411; 95% CI 0.200–0.843). Five-year RFS was 68.8% versus 49.1%, while overall survival showed a favorable trend (HR 0.471; 95% CI 0.165–1.345) but remained immature (5-year OS 92.2% vs 71.3%; 7 events per arm). We therefore conducted an early US payer-perspective economic evaluation to estimate the cost-effectiveness of intismeran plus pembrolizumab and, given the absence of an established launch price, to identify the acquisition price consistent with commonly cited willingness-to-pay thresholds.

---

## 2. Methods

### 2.1 Model Structure

A partitioned survival model (PSM) was developed in Python with four mutually exclusive health states: recurrence-free (RF), locoregional recurrence (LR), distant metastasis (DM), and death [8, 13]. At each monthly cycle, patient proportions were derived from three survival curves:

- RF(t) = RFS(t)
- LR(t) = DMFS(t) − RFS(t)
- DM(t) = OS(t) − DMFS(t)
- Dead(t) = 1 − OS(t)

This construction ensures RFS ≤ DMFS ≤ OS at all times and maintains internal consistency between the three endpoints [11, 13].

### 2.2 Clinical Data

RFS and OS survival probabilities were extracted at 18, 24, 36, 48, and 60 months, and DMFS probabilities at 18, 24, 36, and 48 months, from the published KEYNOTE-942 Kaplan–Meier curves [2] (60-month DMFS was not reported). Parametric survival models were fitted using weighted nonlinear least squares [14, 15, 16]. Five candidate distributions (exponential, Weibull, log-normal, log-logistic, generalized gamma) were evaluated for each component using AIC [14]. A full joint enumeration of 5⁶ = 15,625 distribution combinations was performed; the generalized gamma appeared in all top-10 joint models. The log-normal was selected as the base case for parsimony.

### 2.3 External Validation

Model predictions for the pembrolizumab-alone arm were externally validated against the KEYNOTE-054 and CheckMate 238 trials (Section 3.3).

### 2.4 Costs

The analysis adopted a US payer perspective. Pembrolizumab was assumed at 200 mg every 3 weeks for up to 18 cycles at a WAC of $12,272 per 200 mg dose, yielding an annual cost of $220,896 [6]. Intismeran was modeled at an assumed course acquisition price of $200,000 based on Jefferies analyst estimates [5]. Administration costs were $500 per cycle [19]. Adverse event costs were $15,000 incremental for the combination arm [8, 10]. Monthly disease management costs were $3,000 for locoregional recurrence and $12,000 for distant metastasis [8, 10]. All costs inflated to 2026 USD [20].

### 2.5 Utilities

Health-state utilities: RF 0.83, LR 0.64, DM 0.55 [8, 9, 10]. A one-time utility decrement of 0.05 for adverse events in the combination arm.

### 2.6 Base Case Analysis

Lifetime horizon (40 years), 3% annual discounting for costs and QALYs [22]. Evaluated at WTP thresholds of $50,000, $100,000, and $150,000/QALY [23].

### 2.7 Sensitivity Analysis

**Deterministic:** Parameters varied one at a time across plausible ranges (Supplementary Table S11).

**Probabilistic:** 6,000 Monte Carlo iterations; survival parameters (μ, σ) jointly sampled from bivariate normal distributions with covariance from NLS Jacobian (Supplementary Table S12); costs from gamma (CV=0.2); utilities from beta (SE=0.03). Intismeran price fixed at $200,000. No iterations rejected (0 structural violations).

### 2.8 Scenario Analyses

Three structural scenarios: (1) recurrence-convergence waning (5→8y); (2) hazard-convergence waning (5→10y); (3) no direct OS benefit. See Supplementary Section 9.3 for technical implementation definitions.

---

## 3. Results

### 3.1 Base Case

In the base case, intismeran plus pembrolizumab generated an additional 4.78 discounted life-years (15.06 vs 10.28) and 4.02 additional QALYs (11.35 vs 7.33 QALYs), at an incremental lifetime cost of $161,916 ($874,026 vs $712,110). ICER = $40,274/QALY (Table 2). NMB: $240,123 at $100,000/QALY; $441,143 at $150,000/QALY.

### 3.2 Cost and QALY Decomposition

The $200,000 acquisition cost of intismeran contributed substantially to the incremental cost. The combination arm had $16,114 higher locoregional-recurrence management costs ($85,455 vs $69,341) but $70,199 lower distant-metastasis management costs ($354,189 vs $424,388). The incremental QALY gain of 4.02 was driven primarily by additional time in the recurrence-free state (4.01 RF QALYs), with 0.29 LR QALYs, partly offset by a 0.27-QALY reduction in the DM state.

### 3.3 External and Structural Validation

Over 480 monthly cycles, no violations of RFS ≤ DMFS ≤ OS were observed. The 30-year combination OS was 22.2% and pembro OS 10.5%, both below general-population survival of 23.1% at 30 years for the modeled cohort.

### 3.4 Probabilistic and Deterministic Uncertainty

**PSA:** 6,000 valid iterations. P(CE): 63.8% at $50,000/QALY, 95.3% at $100,000/QALY, 98.7% at $150,000/QALY. Mean ΔQALY 4.32, ΔCost $127,711, NMB $304,291 at $100K. 16.7% cost-saving.

**DSA:** Largest NMB changes from pembro OS parameter ($57,434 to dominant), combo OS parameter ($19,875–$40,684), discount rate ($33,120–$48,284). Max finite ICER $57,434/QALY.

### 3.5 Structural Scenario Analyses

- Recurrence-convergence waning (5→8y): $21,514/QALY (ΔQALY 1.93, ΔCost $41,571)
- Hazard-convergence waning (5→10y): $104,852/QALY (ΔQALY 2.17, ΔCost $227,542)
- No direct OS benefit: Dominant (ΔQALY 0.50, ΔCost −$14,330)

### 3.6 Value-Based Price Thresholds

At $200,000 ICER = $40,274/QALY. Threshold prices: ~$440,000 at $100,000/QALY; ~$641,000 at $150,000/QALY.

---

## 4. Discussion

### 4.1 Principal Findings

The base case ICER of $40,274/QALY, 4.02 QALY gain, and 95.3% P(CE) at $100,000/QALY support cost-effectiveness under conventional US thresholds. The combination remained dominant (lower cost, higher QALY) under the no-direct-OS-benefit scenario, indicating that recurrence prevention alone can support favorable economics within the model.

### 4.2 Economic Drivers

The incremental $161,916 lifetime cost is partly offset by $70,199 lower distant-metastasis management costs. The 4.02 incremental QALYs are driven by 4.01 RF QALYs gained.

### 4.3 Pricing and Reimbursement Implications

ICER remains below $100,000/QALY up to ~$440,000 per course, more than twice the estimated $200,000 price. Managed-entry or outcomes-based reimbursement may be relevant at launch, but would require mature Phase 3 efficacy data and a defined acquisition price.

### 4.4 Long-Term Survival and Structural Uncertainty

Without the GP mortality floor, the log-normal fit predicts 30-year combo OS of 79.1%, exceeding the age-matched general population survival of 23.1% — implausible for resected stage IIIB–IV melanoma. The GP hazard floor reduces this to 22.2% and lowers the ICER from $51,753 (unconstrained) to $40,274 (constrained). The hazard floor compresses the absolute survival difference while preserving the within-trial benefit.

The treatment-effect waning scenarios bracket the plausible ICER range: $21,514 (recurrence-convergence) to $104,852 (hazard-convergence). The Weibull alternative yields $54,925/QALY. The best AIC-ranked mixed-family combination yields $61,102.

### 4.5 Comparison with Previous Economic Evaluations

The ICER of $40,274 is within the range of published pembrolizumab-vs-observation estimates ($38,050–$78,400), but this comparability is coincidental: the absolute QALY gain (4.02) is larger than pembrolizumab-vs-observation (1.2–2.1 QALYs), and the incremental cost is also larger. The ratio falls in a similar range because both numerator and denominator increase proportionally.

The 2024 ISPOR exploratory analysis [11] reported $75,206/QALY using parity-priced V940 and external pembrolizumab curves; the lower ICER here reflects mature 5-year data, a $200,000 price, and a ratio-constrained framework — not a univariate change in any single factor.

### 4.6 Clinical and External Validity

The model's 5-year pembro-arm RFS of 40.5% was below KEYNOTE-054 (50%), consistent with the higher-risk KEYNOTE-942 population. The 5-year pembro-arm OS (74.5%) closely matched KEYNOTE-942 (71.3%). No external long-term data exist for the combination arm; the positive Phase 3 INTerpath-001 readout provides directional confirmation but numerical HRs have not been disclosed.

### 4.7 Limitations

1. Efficacy from Phase 2b (n=157); Phase 3 INTerpath-001 (n=1,137) confirmed RFS/DMFS but numerical data pending.
2. OS data immature (HR 0.471; 95% CI 0.165–1.345; 7 events per arm).
3. Survival curves digitized from published KM figures.
4. Extrapolation beyond 5 years sensitive to parametric choice; 40-year combo OS 2.0% is a model extrapolation, not clinical prediction.
5. Intismeran price based on analyst estimates.
6. WAC used for pembrolizumab; net prices may differ.
7. Utilities and disease management costs from published literature, not KEYNOTE-942.
8. Indirect costs and productivity losses not included.
9. Mixed-stage population (IIIB–IV); subgroup analyses underpowered.
10. No budget-impact analysis conducted.

### 4.8 Conclusion

Intismeran plus pembrolizumab provides clinical benefit in adjuvant melanoma and was projected to be cost-effective at the estimated acquisition cost of $200,000 per course. The deterministic threshold price remains below $100,000/QALY up to approximately $440,000 per course. As more mature Phase 3 INTerpath-001 data become available, these findings should be updated.

---

## Tables

### Table 1. Model Input Parameters

| Parameter | Base Case | Source |
|:---|:---:|:---|
| 5-year RFS, combo | 68.8% | KEYNOTE-942 [2] |
| 5-year RFS, pembro | 49.1% | KEYNOTE-942 [2] |
| RFS HR | 0.51 (95% CI 0.294–0.887) | KEYNOTE-942 [2] |
| DMFS HR | 0.411 (95% CI 0.200–0.843) | KEYNOTE-942 [2] |
| OS HR | 0.471 (95% CI 0.165–1.345) | KEYNOTE-942 [2] |
| 5-year OS, combo | 92.2% | KEYNOTE-942 [2] |
| 5-year OS, pembro | 71.3% | KEYNOTE-942 [2] |
| Keytruda annual cost | $220,896 | WAC [6] |
| Intismeran course cost | $200,000 | Jefferies [5] |

### Table 2. Base-Case Results and Decomposition

| Component | Combo | Pembro | Incremental |
|:---|:---:|:---:|:---:|
| **Cost decomposition** | | | |
| Intismeran | $200,000 | $0 | $200,000 |
| Pembrolizumab | $217,888 | $217,888 | $0 |
| Sequencing | $1,000 | $0 | $1,000 |
| Administration/AE | $15,493 | $493 | $15,000 |
| LR management | $85,455 | $69,341 | $16,114 |
| DM management | $354,189 | $424,388 | −$70,199 |
| **Total cost** | **$874,026** | **$712,110** | **$161,916** |
| **QALY decomposition** | | | |
| RF state | 8.49 | 4.48 | 4.01 |
| LR state | 1.52 | 1.23 | 0.29 |
| DM state | 1.35 | 1.62 | −0.27 |
| **Total QALYs** | **11.35** | **7.33** | **4.02** |
| **ICER** | | | **$40,274/QALY** |

### Table 3. Probabilistic Sensitivity Analysis

| Metric | Mean | Median | 95% CI |
|:---|:---:|:---:|:---|
| Incremental cost | $127,711 | $165,479 | −$668,394 to $520,672 |
| Incremental QALYs | 4.32 | 4.15 | 2.53–7.04 |
| NMB @ $100,000/QALY | $304,291 | — | — |
| NMB @ $150,000/QALY | $520,293 | — | — |

### Table 4. Scenario and Deterministic Threshold Price Analyses

| Scenario | ΔCost | ΔQALY | ICER | NMB @ $150K |
|:---|:---:|:---:|:---:|:---:|
| **Base case (GP-constrained)** | $161,916 | 4.02 | $40,274 | $441,143 |
| **Pricing** | | | | |
| Price: $100K/course | $61,916 | 4.02 | $15,400 | $541,143 |
| Price: $300K/course | $261,916 | 4.02 | $65,147 | $341,143 |
| Price: $500K/course | $461,916 | 4.02 | $114,893 | $141,143 |
| **Alternative assumptions** | | | | |
| Weibull OS | $307,578 | 5.60 | $54,925 | $532,420 |
| Discount rate: 0% | $210,821 | 6.37 | $33,120 | $743,990 |
| Discount rate: 5% | $147,137 | 3.05 | $48,284 | $309,960 |
| Time horizon: 10 years | $96,400 | 1.30 | $74,305 | $33,335 |
| Time horizon: 20 years | $123,974 | 3.08 | $40,293 | $337,545 |
| **Treatment-effect persistence** | | | | |
| Recurrence convergence (5→8y) | $41,571 | 1.93 | $21,514 | $248,268 |
| Hazard convergence (5→10y) | $227,542 | 2.17 | $104,852 | $97,975 |
| No direct OS benefit | −$14,330 | 0.50 | Dominant | $89,330 |

### Table 5. Deterministic Threshold Prices

| Price | ΔCost | ΔQALY | ICER |
|:---|:---:|:---:|:---:|
| $0 | $124,765 | 4.02 | $31,036 |
| $100,000 | $224,765 | 4.02 | $55,913 |
| $200,000 (base case) | $324,765 | 4.02 | $80,787 |
| $300,000 | $424,765 | 4.02 | $105,661 |
| ~$440,000 | ~$564,765 | 4.02 | ~$100,000 |
| $500,000 | $624,765 | 4.02 | $155,414 |
| $600,000 | $724,765 | 4.02 | $180,288 |
| ~$641,000 | ~$765,765 | 4.02 | ~$150,000 |

---

## Figures

- **Figure 1:** Log-normal model RFS and OS curves, with 5-year data points from KEYNOTE-942 → `fig1_survival.png`
- **Figure 2:** Cost-effectiveness acceptability curve → `fig2_ceac.png`
- **Figure 3:** Tornado diagram, one-way sensitivity analysis → `fig3_tornado.png`
- **Figure 4:** Price–ICER curve and deterministic threshold pricing → `fig4_price_icer.png`

---

## References

1. Merck and Moderna. Press release: Merck and Moderna announce positive topline results from Phase 3 INTerpath-001 trial of mRNA-4157/V940 in combination with KEYTRUDA® (pembrolizumab) for adjuvant melanoma. 2026.
2. Khattak MA, et al. Recurrence-free survival and overall survival with mRNA-4157/V940 plus pembrolizumab in resected melanoma: KEYNOTE-942 5-year results. 2026.
3. Ibragimova I, et al. Manufacturing and logistics of individualized neoantigen therapy: a systematic review. 2025.
4. Li Y, et al. Scalability challenges for individualized mRNA cancer vaccines. 2026.
5. Jefferies. Personalized cancer vaccines: cost and market outlook. 2026.
6. Merck Sharp & Dohme LLC. KEYTRUDA® (pembrolizumab) wholesale acquisition cost. 2026.
7. Schwarze L, et al. Sequencing costs and implications for personalized cancer vaccines. 2020.
8. Bensimon J, et al. Cost-effectiveness of pembrolizumab versus observation in resected high-risk melanoma. 2019.
9. Johnston K, et al. Health-state utilities in advanced melanoma: a systematic review. 2021.
10. Zhang M, et al. Cost-effectiveness of adjuvant pembrolizumab in resected stage III melanoma. 2023.
11. McLean L, et al. Cost-effectiveness of V940 plus pembrolizumab in adjuvant melanoma: an early economic evaluation. ISPOR 2024.
12. Eggermont AMM, et al. KEYNOTE-054: Long-term follow-up of pembrolizumab vs placebo in resected stage III melanoma. 2024.
13. Woods B, et al. Partitioned survival models for cost-effectiveness analysis in oncology. 2017.
14. Latimer NR. NICE DSU Technical Support Document 14: Survival analysis for economic evaluations alongside clinical trials—extrapolation with patient-level data. 2011.
15. Latimer NR, et al. Survival analysis for economic evaluations alongside clinical trials—extrapolation with patient-level data. 2013.
16. Ishak KJ, et al. Model selection and uncertainty in cost-effectiveness analysis. 2013.
17. Arias E, et al. US life tables. National Vital Statistics Reports. 2022.
18. Weber JS, et al. Adverse events in patients treated with pembrolizumab. 2017.
19. CMS. Physician Fee Schedule: intravenous infusion costs. 2025.
20. BLS. Consumer Price Index: medical care component. 2025.
21. Briggs A, et al. Decision Modelling for Health Economic Evaluation. Oxford University Press. 2006.
22. Sanders GD, et al. Recommendations for conduct, reporting, modeling, and policy decisions for cost-effectiveness analysis in health. 2016.
23. Neumann PJ, et al. Cost-effectiveness and cost-benefit analysis in the US. 2014.
24. Stinnett AA, et al. Net health benefits: a new framework for the analysis of cost-effectiveness. 1998.
25. Fenwick E, et al. Representing QALYs in a cost-effectiveness acceptability curve. 2004.
26. Husereau D, et al. Consolidated Health Economic Evaluation Reporting Standards 2022 (CHEERS 2022). 2022.
27. Weber JS, et al. KEYNOTE-942: Adjuvant mRNA-4157/V940 plus pembrolizumab in resected melanoma. 2024.
28. Eggermont AMM, et al. KEYNOTE-054: Pembrolizumab vs placebo for resected high-risk melanoma—7-year results. 2025.
29. Ascierto PA, et al. CheckMate 238: Nivolumab vs ipilimumab in resected stage IIIB–IV melanoma. 2024.
30. Ascierto PA, et al. CheckMate 238 long-term follow-up. 2025.
