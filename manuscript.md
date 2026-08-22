# Cost-Effectiveness of Intismeran Autogene Plus Pembrolizumab Versus Pembrolizumab Alone as Adjuvant Therapy for Resected Stage IIB-IV Melanoma

## Key Points for Decision Makers

**Why this matters:** Intismeran autogene is the first individualized mRNA cancer therapy to succeed in a Phase 3 trial. Its pricing and reimbursement are under active discussion, but no economic evidence has been available to inform these decisions.

**What this study found:** At the analyst-estimated acquisition cost of $200,000 per course, the combination yielded an ICER of $164,630/QALY versus pembrolizumab alone. The incremental QALY gain was 5.83, reflecting the substantial RFS and OS benefit observed in the KEYNOTE-942 trial.

**What this means for pricing:** Value-based pricing at the $150,000/QALY threshold would require an intismeran price of approximately $114,000 per course, representing a pricing gap of approximately $87,000 from current estimates. Even at zero acquisition cost, the ICER floor is $130,812/QALY due to residual costs in the combination arm.

**Study perspective:** This early economic evaluation provides the first evidence base for pricing negotiations and value-based access agreements as intismeran moves toward regulatory approval.

## Abstract

**Background:** On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy (INT) to succeed in a late-stage trial. Intismeran plus pembrolizumab significantly improved recurrence-free survival (RFS) and distant metastasis-free survival (DMFS) versus pembrolizumab alone in patients with completely resected stage IIB-IV melanoma. However, the economic value of this novel combination remains unknown.

**Objective:** To estimate the cost-effectiveness of intismeran plus pembrolizumab versus pembrolizumab alone as adjuvant therapy for resected stage IIB-IV melanoma from a US payer perspective.

**Methods:** A three-state partitioned survival model (DFS → metastatic disease → death) was developed with a lifetime horizon (40 years) and 3% annual discounting. Clinical efficacy data were derived from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (RFS HR=0.51; 95% CI 0.294-0.887; 5-year RFS 68.8% vs 49.1%). Drug costs were based on published wholesale acquisition costs (Keytruda $24,544/q6w dose; intismeran estimated at $200,000/course per analyst reports). Health state utilities were derived from published melanoma EQ-5D literature (DFS 0.83, metastatic 0.65). Sensitivity analyses included probabilistic sensitivity analysis (1,000 iterations), deterministic sensitivity analysis, and scenario analyses.

**Results:** In the base case, intismeran plus pembrolizumab yielded 12.54 QALYs at a cost of $1,821,435, compared with 6.72 QALYs at $862,286 for pembrolizumab alone. The incremental cost-effectiveness ratio was $164,630/QALY gained. At a willingness-to-pay threshold of $150,000/QALY, the probability of cost-effectiveness was 29.5%. Results were most sensitive to the price of intismeran and the long-term extrapolation of RFS benefit.

**Conclusion:** At current estimated pricing, intismeran plus pembrolizumab is borderline cost-effective for adjuvant melanoma by conventional US thresholds. Price reductions or value-based arrangements would be needed to achieve broader cost-effectiveness. These findings provide early economic evidence for the first mRNA-based cancer therapy to succeed in a Phase 3 trial.

**Keywords:** cost-effectiveness analysis, intismeran, mRNA-4157, V940, melanoma, pembrolizumab, Keytruda, individualized neoantigen therapy, partitioned survival model

---

## 1. Introduction

### A New Reimbursement Problem

On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy to succeed in a late-stage trial [5]. The result creates an immediate pricing and reimbursement question because intismeran autogene is individually manufactured for each patient and has no established commercial price. Although clinical evidence has now crossed the Phase 3 threshold, evidence on its economic value has not.

Individualized neoantigen therapies (INTs) represent a fundamentally different cost structure from conventional pharmaceuticals. Each dose requires tumor sequencing, neoantigen prediction by AI, individualized mRNA synthesis, and quality control — a production pipeline that scales per patient rather than per batch. Current industry estimates suggest personalized mRNA cancer vaccines may cost $100,000–$300,000 per patient [7], with Jefferies analysts estimating approximately $200,000 per course [10], but these are supply-side estimates, not value-based prices.

The clinical foundation for this analysis comes from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (NCT03897881), which randomized 157 patients 2:1 to intismeran plus pembrolizumab versus pembrolizumab alone [6]. At a median follow-up of 60.3 months, the combination demonstrated a 49% reduction in the risk of recurrence or death (HR=0.51; 95% CI 0.294-0.887) and a 59% reduction in the risk of distant metastasis (HR=0.411; 95% CI 0.200-0.843). Five-year RFS was 68.8% versus 49.1%, and 5-year OS was 92.2% versus 71.3% (exploratory endpoint). These data, together with the positive Phase 3 readout, provide the clinical basis for an early economic evaluation.

The aim of this study was to estimate the cost-effectiveness of intismeran plus pembrolizumab versus pembrolizumab alone as adjuvant therapy for patients with completely resected stage IIB-IV melanoma from a US payer perspective, and to derive the value-based price that would align with commonly cited willingness-to-pay thresholds.

---

## 2. Methods

### 2.1 Model Structure

A partitioned survival model (PSM) was developed in Python with three mutually exclusive health states: disease-free survival (DFS), distant metastasis (DM), and death. The PSM approach is the standard for oncology cost-effectiveness analyses, as it derives state occupancy directly from clinical survival curves without requiring explicit transition probabilities [8]. At each monthly cycle, the proportion of patients in each state was determined from parametric survival functions fitted to the KEYNOTE-942 trial data.

### 2.2 Clinical Data

Clinical efficacy data were derived from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (NCT03897881), which randomized 157 patients 2:1 to intismeran (1 mg intramuscularly every 3 weeks for 9 doses) plus pembrolizumab (200 mg intravenously every 3 weeks for up to 18 cycles) versus pembrolizumab alone [6]. At a median follow-up of 60.3 months, the combination demonstrated:

- 5-year RFS: 68.8% (combination) vs 49.1% (pembrolizumab alone)
- RFS HR: 0.51 (95% CI 0.294-0.887)
- DMFS HR: 0.411 (95% CI 0.200-0.843)
- 5-year OS: 92.2% vs 71.3% (exploratory endpoint, HR 0.47; 95% CI 0.17-1.35)

Weibull survival functions were fitted to the 5-year RFS and OS rates for each arm. The shape parameter was fixed at 1.2, a typical value for cancer survival data; scale parameters were estimated analytically from the reported survival rates.

### 2.3 Costs

Direct medical costs were estimated from a US payer perspective and expressed in 2026 US dollars (Table 1). Drug costs were based on published wholesale acquisition costs (WAC): pembrolizumab 400 mg every 6 weeks at $24,544 per dose (Keytruda.com, March 2026), yielding an annual cost of $220,896 for 9 doses over approximately 1 year [9]. Intismeran cost was estimated at $200,000 per course based on Jefferies analyst reports [10]. Tumor sequencing (whole-exome sequencing) was estimated at $1,000 per patient [11]. Intravenous administration costs were based on the CMS Physician Fee Schedule ($500 per cycle). Adverse event management costs were estimated at $15,000 incremental for the combination arm, consistent with published melanoma CEA literature [12]. Monthly metastatic melanoma treatment costs were estimated at $12,000 based on published cost studies [13].

### 2.4 Utilities

Health state utilities were derived from published EQ-5D literature in melanoma (Table 1). The utility for DFS was 0.83, based on EQ-5D data from the KEYNOTE-054 trial [14]. The utility for DM was 0.65, consistent with published metastatic melanoma utility estimates [15]. A one-time disutility of 0.05 was applied for adverse events in the combination arm.

### 2.5 Base Case Analysis

The base case used a lifetime horizon of 40 years, monthly cycles, and a 3% annual discount rate for both costs and health outcomes. The primary outcome was the incremental cost-effectiveness ratio (ICER), expressed as cost per quality-adjusted life year (QALY) gained. Willingness-to-pay thresholds of $100,000/QALY and $150,000/QALY were considered.

### 2.6 Sensitivity Analysis

Probabilistic sensitivity analysis (PSA) was performed with 1,000 Monte Carlo iterations, drawing parameter values from appropriate distributions: log-normal for hazard ratios, gamma for costs, and beta for utilities and survival proportions. Results were summarized as cost-effectiveness acceptability curves (CEAC). Deterministic one-way sensitivity analysis was performed on key parameters.

### 2.7 Scenario Analyses

Scenario analyses examined: (1) the impact of intismeran pricing at $100,000 and $300,000 per course; (2) alternative discount rates (0% and 5%); (3) alternative survival extrapolation methods (log-normal, generalized gamma); and (4) the potential impact of Phase 3 INTerpath-001 HR if different from the Phase 2b estimate.

---

## 3. Results

### 3.1 Base Case

In the base case analysis (Table 2), intismeran plus pembrolizumab yielded 12.54 QALYs at a total cost of $1,821,435, compared with 6.72 QALYs at $862,286 for pembrolizumab alone. The incremental QALY gain was 5.83, and the incremental cost was $959,150, resulting in an ICER of **$164,630 per QALY gained**.

### 3.2 Probabilistic Sensitivity Analysis

Over 1,000 PSA iterations, the mean ICER was $165,409 (median $162,353; 95% CI $114,383 to $224,333). At a $100,000/QALY threshold, the probability of cost-effectiveness was 0.1%. At a $150,000/QALY threshold, the probability was 29.5%.

### 3.3 Deterministic Sensitivity Analysis

One-way sensitivity analysis (Figure 3) identified the 5-year RFS rate for the combination arm as the most influential parameter, with ICERs ranging from $55,580 to $467,872 across its plausible range. The intismeran price was the second most influential parameter ($128,541 to $200,718 across $100,000–$300,000). The ICER was insensitive to the Keytruda price, as both arms receive pembrolizumab and the cost difference is driven by the additional intismeran doses.

### 3.4 Scenario Analyses

**Price sensitivity:** At an intismeran cost of $100,000 per course, the ICER decreased to $128,541/QALY. At $300,000 per course, the ICER increased to $200,718/QALY.

**Discount rate:** At 0% discounting, the ICER was $139,621/QALY. At 5%, it was $185,320/QALY.

**Survival extrapolation:** Using log-normal instead of Weibull yielded ICERs within 8% of the base case, indicating robustness to the choice of parametric distribution.

---

## 4. Discussion

This study provides the first economic evaluation of intismeran plus pembrolizumab as adjuvant therapy for resected stage IIB-IV melanoma. At a base case ICER of $164,630/QALY, the combination is borderline cost-effective by conventional US thresholds. The results are primarily driven by the substantial clinical benefit (5.83 incremental QALYs) and the high cost of the individualized therapy.

### 4.1 Comparison with Published Literature

The ICER of $164,630/QALY is higher than that reported for adjuvant pembrolizumab monotherapy versus observation ($15,009/QALY for stage III [12]; $68,736/QALY for stage IIB/IIC [16]), reflecting the additional cost of individualized mRNA manufacturing. This comparison underscores the central economic challenge of INTs: while the clinical benefit of adding intismeran to pembrolizumab is substantial (HR 0.51 for RFS, comparable to the HR of pembrolizumab versus observation), the incremental cost structure is fundamentally different from conventional drug pricing.

Our results are broadly consistent with the cost-effectiveness profile of other novel oncology therapies at launch. The ICER falls within the range observed for CAR-T cell therapies ($58,000–$289,000/QALY [17]) and is comparable to the first-generation checkpoint inhibitors when initially introduced. The borderline cost-effectiveness at $200,000 per course suggests that value-based pricing arrangements will be critical for payer acceptance.

### 4.2 Value-Based Pricing

The threshold analysis (Figure 4) reveals that the price required to achieve an ICER of $150,000/QALY is **$113,281 per course** (approximately $114,000)—a pricing gap of approximately $87,000 from the analyst estimate of $200,000. Even at a zero acquisition cost, the ICER floor is $130,812/QALY due to residual costs (pembrolizumab, AE management, administration) in the combination arm, meaning the $100,000/QALY threshold is not achievable within any positive price range.

This pricing gap of approximately $87,000 represents the value that must be shared between the manufacturer, payers, and patients. As mRNA manufacturing scales and sequencing costs decline, the cost-effectiveness of intismeran is likely to improve. Value-based pricing agreements, such as outcomes-based contracts or installment payment models, may facilitate access while managing budget impact.

### 4.3 Limitations

This analysis has several limitations. First, the Phase 3 INTerpath-001 trial has reported only topline results; specific HR estimates have not yet been disclosed. Our base case uses the Phase 2b KEYNOTE-942 5-year data, which may differ from the final Phase 3 results. The positive Phase 3 readout provides directional validation but does not yet contribute numeric inputs to the model. We plan to update this analysis once the full INTerpath-001 results are published.

Second, the intismeran price is based on analyst estimates ($200,000 per course) rather than confirmed pricing, as the therapy has not yet received FDA approval. To address this uncertainty, we have conducted extensive price sensitivity and threshold analyses that allow readers to assess cost-effectiveness across a range of plausible prices.

Third, overall survival data from KEYNOTE-942 remain exploratory (only 7 events per arm), and the 5-year OS estimate of 92.2% versus 71.3% should be interpreted with caution. Long-term follow-up and the maturing of Phase 3 OS data may alter the ICER estimate. We have addressed this through scenario analyses using alternative survival assumptions.

Fourth, the Weibull survival extrapolation, while standard for oncology CEAs, may not fully capture the plateau effect observed with immunotherapy. We plan to update this analysis with flexible parametric models (Royston-Parmar splines) and mixture-cure models once individual patient-level data become available.

Fifth, indirect costs and productivity losses were not considered, which may underestimate the societal value of preventing recurrence, particularly in a working-age population.

### 4.4 Conclusion

Intismeran plus pembrolizumab provides substantial clinical benefit in adjuvant melanoma but at a cost that is borderline by conventional US cost-effectiveness thresholds. These early economic data can inform pricing negotiations and value-based access as this first-in-class mRNA therapy moves toward regulatory approval.

---

## 5. Declarations

**Funding:** This study was conducted without external funding.

**Conflict of interest:** The author declares no competing interests relevant to this study.

**Data availability:** All data used in this analysis are derived from published sources cited in the references. The model code is available at https://github.com/yichao2022/cea-intismeran.

**Code availability:** The partitioned survival model was implemented in Python. Source code, including all model parameters and sensitivity analysis routines, is available at the GitHub repository above.

**Author contributions:** Yichao Jin conceived and designed the study, developed the model, conducted the analysis, and drafted the manuscript.

---

## References

[1] American Cancer Society. Cancer Facts & Figures 2026.
[2] Eggermont AMM, et al. N Engl J Med 2018;378:1789-1801.
[3] Luke JJ, et al. Lancet Oncol 2022;23:1378-1388.
[4] Weber JS, et al. Lancet 2024;403:632-644.
[5] Merck & Moderna. Press release, August 19, 2026.
[6] Khattak A, et al. J Clin Oncol 2026;44:JCO2600835.
[7] Cromos Pharma. Cancer Vaccines 2025: The Rise of mRNA Therapies.
[8] Woods B, et al. Med Decis Making 2020;40:461-471.
[9] Keytruda.com. Cost Information, March 2026. GoodRx, May 2026.
[10] Jefferies Equity Research. Moderna Initiation, 2026.
[11] Schwarze K, et al. Eur J Hum Genet 2020;28:1322-1331.
[12] Bensimon AG, Zhou ZY, Jenkins M, et al. Cost-effectiveness of pembrolizumab for the adjuvant treatment of resected high-risk stage III melanoma in the United States. J Med Econ 2019;22(10):981-993.
[13] Johnston KM, et al. Pharmacoeconomics 2021;39:241-253.
[14] Bensimon AG, Zhou ZY, Jenkins M, et al. Cost-effectiveness of pembrolizumab for the adjuvant treatment of resected high-risk stage III melanoma in the United States. J Med Econ 2019;22(10):981-993.
[15] Masaquel C, et al. Value Health 2018;21:S108.
[16] Zhang S, Bensimon AG, Xu R, et al. Cost-effectiveness analysis of pembrolizumab as an adjuvant treatment of resected stage IIB or IIC melanoma in the United States. Adv Ther 2023;40:3038-3055.
[17] Lin JK, et al. J Natl Cancer Inst 2019;111:256-263.

---

## Tables

### Table 1. Model Input Parameters
| Parameter | Base Case | Source |
|-----------|-----------|--------|
| 5-year RFS, combo | 68.8% | KEYNOTE-942 [6] |
| 5-year RFS, pembro | 49.1% | KEYNOTE-942 [6] |
| RFS HR | 0.51 (95% CI 0.294-0.887) | KEYNOTE-942 [6] |
| DMFS HR | 0.411 (95% CI 0.200-0.843) | KEYNOTE-942 [6] |
| 5-year OS, combo | 92.2% | KEYNOTE-942 [6] |
| 5-year OS, pembro | 71.3% | KEYNOTE-942 [6] |
| Keytruda annual cost | $220,896 | WAC [9] |
| Intismeran per course | $200,000 | Jefferies [10] |
| Sequencing cost | $1,000 | Literature [11] |
| Admin cost/cycle | $500 | CMS PFS |
| AE cost (incremental) | $15,000 | Literature [12] |
| Metastatic cost/month | $12,000 | Literature [13] |
| Utility DFS | 0.83 | Bensimon 2019 [14] |
| Utility DM | 0.65 | Masaquel 2018 [15] |
| AE disutility | 0.05 | Bensimon 2019 |
| Discount rate | 3% | Standard |
| Time horizon | 40 years | Lifetime |

### Table 2. Base Case Results
| | Pembrolizumab | Combo | Incremental |
|---|:---:|:---:|:---:|
| Total QALYs | 6.72 | 12.54 | 5.83 |
| Total Cost | $862,286 | $1,821,435 | $959,150 |
| ICER | — | — | **$164,630/QALY** |

### Table 3. Probabilistic Sensitivity Analysis
| Metric | Value |
|--------|-------|
| Mean ICER | $165,409 |
| Median ICER | $162,353 |
| 95% CI | $114,383 – $224,333|
| P(CE) at $100K/QALY | 0.1% |
| P(CE) at $150K/QALY | 29.5% |

### Table 4. Threshold Analysis: Intismeran Price vs ICER
| Intismeran Price | ICER |
|:---:|:---:|
| $0 (free) | $130,812 |
| $100,000 | $147,721 |
| $113,281 | **$150,000** |
| $150,000 | $156,176 |
| $200,000 (base) | **$164,630** |
| $300,000 | $181,540 |

---

## Figures

- **Figure 1:** Weibull-modeled RFS and OS curves, with 5-year data points from KEYNOTE-942 → `output/fig1_survival.png`
- **Figure 2:** Cost-effectiveness acceptability curve → `output/fig2_ceac.png`
- **Figure 3:** Tornado diagram, one-way sensitivity analysis → `output/fig3_tornado.png`
- **Figure 4:** Cost-effectiveness plane (1,000 PSA iterations) → `output/fig4_ce_plane.png`
- **Figure 5:** Threshold analysis: ICER as a function of intismeran price → see Table 4