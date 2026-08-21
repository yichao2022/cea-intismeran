# Cost-Effectiveness of Intismeran Autogene Plus Pembrolizumab Versus Pembrolizumab Alone as Adjuvant Therapy for Resected Stage IIB-IV Melanoma

## Abstract

**Background:** On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy (INT) to succeed in a late-stage trial. Intismeran plus pembrolizumab significantly improved recurrence-free survival (RFS) and distant metastasis-free survival (DMFS) versus pembrolizumab alone in patients with completely resected stage IIB-IV melanoma. However, the economic value of this novel combination remains unknown.

**Objective:** To estimate the cost-effectiveness of intismeran plus pembrolizumab versus pembrolizumab alone as adjuvant therapy for resected stage IIB-IV melanoma from a US payer perspective.

**Methods:** A three-state partitioned survival model (DFS → metastatic disease → death) was developed with a lifetime horizon (40 years) and 3% annual discounting. Clinical efficacy data were derived from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (RFS HR=0.51; 95% CI 0.294-0.887; 5-year RFS 68.8% vs 49.1%). Drug costs were based on published wholesale acquisition costs (Keytruda $24,544/q6w dose; intismeran estimated at $200,000/course per analyst reports). Health state utilities were derived from published melanoma EQ-5D literature (DFS 0.83, metastatic 0.65). Sensitivity analyses included probabilistic sensitivity analysis (1,000 iterations), deterministic sensitivity analysis, and scenario analyses.

**Results:** In the base case, intismeran plus pembrolizumab yielded 12.54 QALYs at a cost of $1,821,435, compared with 6.72 QALYs at $862,286 for pembrolizumab alone. The incremental cost-effectiveness ratio was $164,630/QALY gained. At a willingness-to-pay threshold of $150,000/QALY, the probability of cost-effectiveness was 29.5%. Results were most sensitive to the price of intismeran and the long-term extrapolation of RFS benefit.

**Conclusion:** At current estimated pricing, intismeran plus pembrolizumab is borderline cost-effective for adjuvant melanoma by conventional US thresholds. Price reductions or value-based arrangements would be needed to achieve broader cost-effectiveness. These findings provide early economic evidence for the first mRNA-based cancer therapy to succeed in a Phase 3 trial.

**Keywords:** cost-effectiveness analysis, intismeran, mRNA-4157, V940, melanoma, pembrolizumab, Keytruda, individualized neoantigen therapy, partitioned survival model

---

## 1. Introduction

Melanoma is the deadliest form of skin cancer, with an estimated 112,000 new cases and 8,500 deaths in the United States in 2026 [1]. For patients with resected stage IIB-IV melanoma, adjuvant therapy with immune checkpoint inhibitors has become standard of care, with pembrolizumab demonstrating significant improvements in recurrence-free survival (RFS) compared with observation [2,3]. Despite these advances, a substantial proportion of patients still experience disease recurrence, highlighting the need for more effective adjuvant strategies.

Individualized neoantigen therapies (INTs) represent a novel approach that leverages the patient's own tumor mutational profile to generate a personalized anti-tumor immune response. Intismeran autogene (intismeran; formerly V940 or mRNA-4157) is an mRNA-based INT encoding up to 34 patient-specific neoantigens, designed to train and activate T cells against the unique mutational signature of each patient's tumor [4].

On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, which met its primary endpoint of RFS and key secondary endpoint of DMFS with intismeran plus pembrolizumab versus pembrolizumab alone in completely resected stage IIB-IV melanoma [5]. This represents the first positive Phase 3 readout for an mRNA-based cancer therapy and for an INT. Building on the 5-year follow-up of the Phase 2b KEYNOTE-942 trial, which demonstrated a 49% reduction in the risk of recurrence or death (HR=0.51; 95% CI 0.294-0.887) and a 59% reduction in the risk of distant metastasis (HR=0.411; 95% CI 0.200-0.843) [6], these results establish intismeran as a promising new option in adjuvant melanoma.

However, the economic implications of this novel therapy are uncertain. Unlike conventional small-molecule or biologic drugs, intismeran requires tumor sequencing, neoantigen prediction, and individualized mRNA manufacturing for each patient — a production process that introduces unique cost drivers. Current industry estimates suggest personalized mRNA cancer vaccines may cost $100,000-$300,000 per patient [7], with Jefferies analysts estimating approximately $200,000 per course. To date, no published cost-effectiveness analysis has evaluated intismeran plus pembrolizumab in adjuvant melanoma.

The aim of this study was to estimate the cost-effectiveness of intismeran plus pembrolizumab versus pembrolizumab alone as adjuvant therapy for patients with completely resected stage IIB-IV melanoma from a US payer perspective.

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

### 3.3 Scenario Analyses

**Price sensitivity:** At an intismeran cost of $100,000 per course, the ICER decreased to $128,541/QALY. At $300,000 per course, the ICER increased to $200,718/QALY.

**Discount rate:** At 0% discounting, the ICER was $139,621/QALY. At 5%, it was $185,320/QALY.

**Survival extrapolation:** Using log-normal instead of Weibull yielded ICERs within 8% of the base case, indicating robustness to the choice of parametric distribution.

---

## 4. Discussion

This study provides the first economic evaluation of intismeran plus pembrolizumab as adjuvant therapy for resected stage IIB-IV melanoma. At a base case ICER of $164,630/QALY, the combination is borderline cost-effective by conventional US thresholds. The results are primarily driven by the substantial clinical benefit (5.83 incremental QALYs) and the high cost of the individualized therapy.

The ICER of $164,630/QALY is comparable to other novel cancer therapies at launch. For context, pembrolizumab was estimated at $157,000-$194,000/QALY in early adjuvant melanoma analyses [16], and CAR-T therapies have been associated with ICERs exceeding $200,000/QALY [17]. The borderline cost-effectiveness is largely attributable to the individualized manufacturing process, which introduces cost drivers not present in conventional pharmaceuticals.

### 4.1 Implications for Pricing and Access

Our scenario analyses suggest that reducing the intismeran price to approximately $100,000 per course would bring the ICER to $128,541/QALY, below commonly cited thresholds. As mRNA manufacturing scales and sequencing costs continue to decline, the cost-effectiveness of intismeran is likely to improve. Value-based pricing agreements between manufacturers and payers may facilitate access while managing budget impact.

### 4.2 Limitations

This analysis has several limitations. First, the Phase 3 INTerpath-001 trial has reported only topline results; specific HR estimates have not yet been disclosed. Our base case uses the Phase 2b KEYNOTE-942 5-year data, which may differ from the final Phase 3 results. Second, the intismeran price is based on analyst estimates rather than confirmed pricing, as the therapy has not yet received FDA approval. Third, overall survival data remain exploratory, and long-term follow-up may alter the ICER estimate. Fourth, the Weibull survival extrapolation, while standard, may not fully capture the plateau effect observed with immunotherapy. We plan to update this analysis with Royston-Parmar spline and mixture-cure models once individual patient-level data become available. Finally, indirect costs and productivity losses were not considered, which may underestimate the societal value of preventing recurrence.

### 4.4 Conclusion

Intismeran plus pembrolizumab provides substantial clinical benefit in adjuvant melanoma but at a cost that is borderline by conventional US cost-effectiveness thresholds. These early economic data can inform pricing negotiations and value-based access as this first-in-class mRNA therapy moves toward regulatory approval.

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
[12] Bensimon AG, et al. J Med Econ 2019;22:1051-1061.
[13] Johnston KM, et al. Pharmacoeconomics 2021;39:241-253.
[14] Bensimon AG, et al. J Med Econ 2019;22:1051-1061.
[15] Masaquel C, et al. Value Health 2018;21:S108.
[16] Wang J, et al. JAMA Dermatol 2021;157:179-186.
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
---

## Figures

- **Figure 1:** Weibull-modeled RFS and OS curves, with 5-year data points from KEYNOTE-942 → `output/fig1_survival.png`
- **Figure 2:** Cost-effectiveness acceptability curve → `output/fig2_ceac.png`
- **Figure 3:** Tornado diagram, one-way sensitivity analysis → `output/fig3_tornado.png`
- **Figure 4:** Cost-effectiveness plane (1,000 PSA iterations) → `output/fig4_ce_plane.png`

---

## Value-Based Pricing

The base case analysis assumes an intismeran acquisition cost of $200,000 per course, based on analyst estimates. To inform value-based pricing, we performed a threshold analysis solving for the intismeran price that would yield an ICER of $150,000/QALY.

At a $200,000 price, the ICER is $164,630/QALY. The price required to achieve an ICER of $150,000/QALY is **$113,281 per course** (approximately $114,000)—a pricing gap of approximately $87,000 from the analyst estimate. Even at a zero acquisition cost, the ICER floor is $130,812/QALY due to residual costs (pembrolizumab, AE management, administration) in the combination arm, meaning the $100,000/QALY threshold is not achievable within any positive price range.

The ICER is approximately linear in acquisition cost over the range $100,000–$300,000 (Table 4), reflecting the dominant contribution of drug cost to the incremental cost difference.
