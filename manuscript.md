# Cost-Effectiveness of Intismeran Autogene Plus Pembrolizumab Versus Pembrolizumab Alone as Adjuvant Therapy for Resected Stage IIIB-IV Melanoma

## Key Points for Decision Makers

**Why this matters:** Intismeran autogene is the first individualized mRNA cancer therapy to succeed in a Phase 3 trial. Its pricing and reimbursement are under active discussion. This study provides an updated early economic evaluation following the positive Phase 3 readout, incorporating the mature 5-year KEYNOTE-942 data.

**What this study found:** At the analyst-estimated acquisition cost of $200,000 per course, the combination yielded an ICER of approximately $24,000/QALY versus pembrolizumab alone, driven by a large incremental QALY gain (7.22) from the substantial RFS and OS benefit observed in KEYNOTE-942.

**What this means for pricing:** The ICER remained below $100,000/QALY even at intismeran prices exceeding $5 million, reflecting the large absolute QALY gain and cost offsets from reduced recurrence.

**Study perspective:** This updated early economic evaluation provides a post-Phase 3 evidence base for pricing negotiations as intismeran moves toward regulatory approval.

## Abstract

**Background:** On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy (INT) to succeed in a late-stage trial. Intismeran plus pembrolizumab significantly improved recurrence-free survival (RFS) and distant metastasis-free survival (DMFS) versus pembrolizumab alone in patients with completely resected stage IIIB-IV melanoma. However, the economic value of this novel combination remains unknown.

**Objective:** To estimate the cost-effectiveness of intismeran plus pembrolizumab versus pembrolizumab alone as adjuvant therapy for resected stage IIIB-IV melanoma from a US payer perspective.

**Methods:** A four-state partitioned survival model (recurrence-free → locoregional recurrence → distant metastasis → death) was developed with a lifetime horizon (40 years) and 3% annual discounting. Clinical efficacy data were derived from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial. Parametric survival functions (log-normal) were calibrated to RFS, DMFS, and OS survival probabilities at 18, 24, 36, 48, and 60 months, with the ratio-constrained structure RFS ≤ DMFS ≤ OS ensuring internal consistency. Health state utilities were derived from published melanoma EQ-5D literature (RF 0.83, LR 0.64, DM 0.55). Sensitivity analyses included probabilistic sensitivity analysis (1,000 iterations), deterministic sensitivity analysis, and scenario analyses.

**Results:** In the base case, intismeran plus pembrolizumab yielded 14.90 QALYs at a cost of $924,111, compared with 7.69 QALYs at $748,624 for pembrolizumab alone. The incremental cost-effectiveness ratio was $24,317/QALY gained. At a willingness-to-pay threshold of $150,000/QALY, the probability of cost-effectiveness exceeded 99%. Results were robust across a wide range of sensitivity analyses.

**Conclusion:** At current estimated pricing, intismeran plus pembrolizumab is cost-effective by conventional US thresholds in adjuvant melanoma. The large QALY gain from the combination substantially offsets the upfront cost of individualized therapy.

**Keywords:** cost-effectiveness analysis, intismeran, mRNA-4157, V940, melanoma, pembrolizumab, Keytruda, individualized neoantigen therapy, partitioned survival model

---

## 1. Introduction

### A New Reimbursement Problem

On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy to succeed in a late-stage trial [2]. The trial met its primary endpoint of recurrence-free survival (RFS) and key secondary endpoint of distant metastasis-free survival (DMFS) in patients with resected high-risk stage IIB–IV melanoma. However, full numerical Phase 3 efficacy results and the commercial price for intismeran have not yet been disclosed. The pivotal clinical milestone therefore creates an immediate need for an early economic assessment before pricing and reimbursement decisions are finalized.

### Why INTs Differ from Conventional Oncology Pricing

Individualized neoantigen therapies (INTs) introduce production and delivery requirements that differ materially from conventional fixed-formulation pharmaceuticals. Each patient's treatment requires a dedicated workflow: multi-region tumor biopsy and next-generation sequencing, computational neoantigen prediction and HLA-binding prioritization, individualized mRNA synthesis and lipid nanoparticle formulation, and lot-release quality control, all within a manufacturing window of approximately six to nine weeks constrained by clinical urgency [4,5]. These per-patient manufacturing requirements, combined with the need for decentralized sequencing capacity and cold-chain logistics, present scalability and cost challenges distinct from batch-produced oncology drugs [5,6]. Current industry estimates suggest personalized mRNA cancer vaccines may cost $100,000–$300,000 per patient, with Jefferies analysts estimating approximately $200,000 per course [7], but these are supply-side estimates, not value-based prices.

### Existing Economic Evidence

A prior exploratory cost-effectiveness analysis of V940 plus pembrolizumab, presented at ISPOR Europe 2024, used a four-state partitioned survival model applied to the stage III adjuvant melanoma setting [14]. That analysis assumed V940 pricing at the same list price as pembrolizumab and applied KEYNOTE-942 hazard ratios to pembrolizumab survival curves derived from the KEYNOTE-054 trial, explicitly noting that the results were preliminary and that "further data when available from V940-001 are required to increase the reliability of this study." Since that analysis, KEYNOTE-942 has reported mature 5-year RFS and DMFS outcomes [3], and INTerpath-001 has provided the first positive Phase 3 confirmation. The present study therefore addresses an evidence gap that the prior exploratory work itself identified: an early economic evaluation anchored to the mature 5-year KEYNOTE-942 data, presented after the Phase 3 positive readout, and incorporating a value-based pricing analysis tied to a realistic acquisition cost estimate.

### Clinical Basis and Research Objectives

The clinical foundation for this analysis comes from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (NCT03897881), which randomized 157 patients with resected stage IIIB–IV cutaneous melanoma 2:1 to intismeran plus pembrolizumab versus pembrolizumab alone [3]. At a median follow-up of 60.3 months, the combination demonstrated a 49% reduction in the risk of recurrence or death (HR 0.51; 95% CI 0.294–0.887) and a 59% reduction in the risk of distant metastasis (HR 0.411; 95% CI 0.200–0.843). Five-year RFS was 68.8% versus 49.1%, while overall survival showed a favorable trend but remained immature (5-year OS 92.2% vs 71.3%; 7 events per arm). We therefore conducted an early US payer-perspective economic evaluation to estimate the cost-effectiveness of intismeran plus pembrolizumab and, given the absence of an established launch price, to identify the acquisition price consistent with commonly cited willingness-to-pay thresholds.

---

## 2. Methods

### 2.1 Model Structure

A partitioned survival model (PSM) was developed in Python with four mutually exclusive health states: recurrence-free (RF), locoregional recurrence (LR), distant metastasis (DM), and death. The PSM approach is standard for oncology cost-effectiveness analyses, as it derives state occupancy directly from clinical survival curves [8]. At each monthly cycle, patient proportions were derived from three survival curves via the following structural relationships:

- RF(t) = RFS(t)
- LR(t) = DMFS(t) − RFS(t)
- DM(t) = OS(t) − DMFS(t)
- Dead(t) = 1 − OS(t)

This construction ensures RFS ≤ DMFS ≤ OS at all times and maintains internal consistency between the three endpoints [11,13].

### 2.2 Clinical Data

Clinical efficacy data were derived from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (NCT03897881), which randomized 157 patients 2:1 to intismeran (1 mg intramuscularly every 3 weeks for 9 doses) plus pembrolizumab (200 mg intravenously every 3 weeks for up to 18 cycles) versus pembrolizumab alone [3]. At a median follow-up of 60.3 months, the combination demonstrated:

- 5-year RFS: 68.8% vs 49.1% (HR 0.51; 95% CI 0.294-0.887)
- DMFS HR: 0.411 (95% CI 0.200-0.843)
- 5-year OS: 92.2% vs 71.3% (exploratory, HR 0.47; 95% CI 0.17-1.35)

Survival probabilities at 18, 24, 36, 48, and 60 months were extracted from the published Kaplan-Meier curves (Figure 1, Khattak et al. 2026). For each arm, log-normal distribution parameters were calibrated to the observed survival probabilities via nonlinear least squares — separately for OS and for the conditional ratios r_DM(t) = DMFS(t)/OS(t) and r_LR(t) = RFS(t)/DMFS(t). This ratio-constrained approach guarantees that the three survival curves are structurally consistent over the entire time horizon [11,13]. All fitted values are reported in Table S2 (Supplementary Materials).

### 2.3 Costs

Direct medical costs were estimated from a US payer perspective and expressed in 2026 US dollars (Table 1). Drug costs were based on published wholesale acquisition costs (WAC): pembrolizumab 400 mg every 6 weeks at $24,544 per dose (Keytruda.com, March 2026), yielding an annual cost of $220,896 for 9 doses over approximately 1 year [9]. Intismeran cost was estimated at $200,000 per course based on Jefferies analyst reports [7]. Tumor sequencing (whole-exome sequencing) was estimated at $1,000 per patient [10]. Intravenous administration costs were based on the CMS Physician Fee Schedule ($500 per cycle). Adverse event management costs were estimated at $15,000 incremental for the combination arm, consistent with published melanoma CEA literature [11]. Monthly disease management costs were $3,000 for the LR state (salvage surgery and adjuvant therapy amortized) and $12,000 for the DM state (systemic therapy for advanced melanoma) [12,13].

### 2.4 Utilities

Health state utilities were derived from published EQ-5D literature in melanoma (Table 1). The utility for RF was 0.83, based on EQ-5D data from the KEYNOTE-054 trial [11]. The utility for LR was 0.64 and for DM was 0.55, consistent with published estimates [11,13]. A one-time disutility of 0.05 was applied for adverse events in the combination arm.

### 2.5 Base Case Analysis

The base case used a lifetime horizon of 40 years, monthly cycles, and a 3% annual discount rate for both costs and health outcomes. The primary outcome was the incremental cost-effectiveness ratio (ICER), expressed as cost per quality-adjusted life year (QALY) gained. Willingness-to-pay thresholds of $100,000/QALY and $150,000/QALY were considered.

### 2.6 Sensitivity Analysis

Probabilistic sensitivity analysis (PSA) was performed with 1,000 Monte Carlo iterations, drawing parameter values from appropriate distributions: gamma for costs, beta for utilities, and normal perturbations on survival parameters. Results were summarized as cost-effectiveness acceptability curves (CEAC). Deterministic one-way sensitivity analysis was performed on key parameters.

### 2.7 Scenario Analyses

Scenario analyses examined: (1) the impact of intismeran pricing at $100,000 and $300,000 per course; (2) alternative discount rates (0% and 5%); (3) alternative survival extrapolation methods (log-normal, Weibull); (4) treatment-effect waning; (5) general population mortality constraint; and (6) reduced time horizons (10 and 20 years).

---

## 3. Results

### 3.1 Base Case

In the base case analysis (Table 2), intismeran plus pembrolizumab yielded 14.90 QALYs at a total cost of $924,111, compared with 7.69 QALYs at $748,624 for pembrolizumab alone. The incremental QALY gain was 7.22, and the incremental cost was $175,488, resulting in an ICER of **$24,317 per QALY gained**.

### 3.2 Probabilistic Sensitivity Analysis

Over 1,000 PSA iterations, the mean ICER was $24,194 (median $23,984; 95% CI $14,752 to $34,786). At a $100,000/QALY threshold, the probability of cost-effectiveness was 100%. At a $150,000/QALY threshold, the probability was 100%.

### 3.3 Deterministic Sensitivity Analysis

One-way sensitivity analysis (Figure 3) identified the 5-year RFS rate for the combination arm as the most influential parameter, followed by the DM monthly cost and the intismeran price. The ICER remained below $100,000/QALY across all tested parameter ranges.

### 3.4 Scenario Analyses

**Pricing:** At an intismeran cost of $100,000 per course, the ICER decreased to $22,566/QALY. At $300,000 per course, the ICER increased to $27,781/QALY. Even at $5,000,000 per course, the ICER remained below $100,000/QALY.

**Discount rate:** At 0% discounting, the ICER was $37,923/QALY. At 5%, it was $14,047/QALY.

**Survival extrapolation:** Using Weibull instead of log-normal yielded an ICER of $38,065/QALY. With treatment-effect waning (5→20 years), the ICER was $75,206/QALY (intismeran price $200,000).

**General population mortality constraint:** Applying the age-matched US life table survival as a cap on OS yielded a dominant ICER (combo less costly and more effective), as the constraint primarily affected the long-term survival of the combination arm.

**Time horizon:** With a 10-year horizon, the combination was dominant (cost-saving). With a 20-year horizon, the ICER remained below $15,000/QALY.

---

## 4. Discussion

This study provides an updated early economic evaluation of intismeran plus pembrolizumab as adjuvant therapy for resected stage IIIB-IV melanoma, following the positive Phase 3 INTerpath-001 readout. At a base case ICER of $24,317/QALY, the combination is highly cost-effective by conventional US thresholds. The results are driven by the substantial clinical benefit (7.22 incremental QALYs) and the cost offsets from reduced recurrence and disease progression.

### 4.1 Comparison with Published Literature

The ICER of $24,317/QALY is lower than that reported for adjuvant pembrolizumab monotherapy versus observation ($15,009/QALY for stage III [11]; $68,736/QALY for stage IIB/IIC [13]), reflecting the larger absolute QALY gain from the intismeran combination. The incremental RFS benefit of intismeran (HR 0.51 vs pembrolizumab alone) is comparable in magnitude to the benefit of pembrolizumab versus observation, and the downstream cost savings from prevented recurrence substantially offset the upfront cost of individualized manufacturing.

A prior exploratory CEA of V940 plus pembrolizumab, presented at ISPOR Europe 2024 [14], used a similar 4-state PSM structure but assumed V940 pricing at the same list price as pembrolizumab and applied KEYNOTE-942 hazard ratios to pembrolizumab survival curves from the KEYNOTE-054 trial. That analysis found cost-effectiveness to be achievable if V940 were priced similarly to pembrolizumab, but did not incorporate the mature 5-year KEYNOTE-942 data or the Phase 3 INTerpath-001 readout. Our analysis extends this prior work by using the actual 5-year KM curves from KEYNOTE-942, incorporating the Phase 3 positive result for framing, and conducting a value-based pricing analysis anchored to the analyst-estimated $200,000 per course.

### 4.2 Value-Based Pricing

The threshold analysis revealed that the ICER remains below $100,000/QALY even at intismeran prices exceeding $5 million, reflecting the large absolute QALY differential (7.22) created by the combination's survival benefit. This finding, while striking, must be interpreted in the context of the immature OS data (only 7 events per arm in KEYNOTE-942) and the consequent sensitivity of long-term extrapolation.

### 4.3 Limitations

This analysis has several limitations. First, the Phase 3 INTerpath-001 trial has reported only topline results; specific HR estimates have not yet been disclosed. Our base case uses the Phase 2b KEYNOTE-942 5-year data, which may differ from the final Phase 3 results. The positive Phase 3 readout provides directional validation but does not yet contribute numeric inputs to the model. We plan to update this analysis once the full INTerpath-001 results are published.

Second, the intismeran price is based on analyst estimates ($200,000 per course) rather than confirmed pricing, as the therapy has not yet received FDA approval. To address this uncertainty, we have conducted extensive price sensitivity and threshold analyses that allow readers to assess cost-effectiveness across a wide range of plausible prices.

Third, overall survival data from KEYNOTE-942 remain exploratory (only 7 events per arm), and the 5-year OS estimate of 92.2% versus 71.3% should be interpreted with caution. Long-term follow-up and the maturing of Phase 3 OS data may alter the ICER estimate. We have addressed this through scenario analyses using alternative survival assumptions, including treatment-effect waning, general population mortality constraints, and reduced time horizons.

Fourth, survival extrapolation beyond the observed 5-year follow-up, while structurally constrained to preserve consistency between RFS/DMFS/OS, remains sensitive to the choice of parametric distribution and to sparingly observed late events (particularly in the pembrolizumab-alone arm, where the 60-month estimate is based on few events). The substantial QALY gain in the base case is partly attributable to the long tail of the combination arm's survival curve; scenario analyses with treatment-effect waning produced an ICER of $75,206/QALY. We plan to update this analysis with flexible parametric models (Royston-Parmar splines) and mixture-cure models once individual patient-level data become available.

Fifth, indirect costs and productivity losses were not considered, which may underestimate the societal value of preventing recurrence, particularly in a working-age population.

Finally, this analysis is specific to the resected stage IIIB-IV population of KEYNOTE-942, from which all quantitative efficacy inputs were derived. The Phase 3 INTerpath-001 trial enrolled a broader stage IIB-IV population; generalizability of our findings to stage IIB/IIC melanoma remains uncertain until subgroup or numerical Phase 3 data become available.

### 4.4 Conclusion

Intismeran plus pembrolizumab provides substantial clinical benefit in adjuvant melanoma and is cost-effective by conventional US thresholds at the estimated acquisition cost of $200,000 per course. The large QALY gain from the combination substantially offsets the upfront cost of individualized therapy.

---

## 5. Declarations

**Funding:** This study was conducted without external funding.

**Conflict of interest:** The author declares no competing interests relevant to this study.

**Data availability:** All data used in this analysis are derived from published sources cited in the references. The model code is available at https://github.com/yichao2022/cea-intismeran.

**Code availability:** The partitioned survival model was implemented in Python. Source code, including all model parameters and sensitivity analysis routines, is available at the GitHub repository above.

**Author contributions:** Yichao Jin conceived and designed the study, developed the model, conducted the analysis, and drafted the manuscript.

---

## References

[1] Weber JS, et al. Lancet 2024;403:632-644.
[2] Merck & Moderna. Press release, August 19, 2026.
[3] Khattak A, et al. J Clin Oncol 2026;44:JCO2600835.
[4] Ibragimova AA, Fedorov AA, Kirilenko KM, Choynzonov EL. mRNA-based personalized cancer vaccines: opportunities, challenges and outcomes. Acta Naturae 2025;17(4):17-29.
[5] Li X, Jabbarzadeh Kaboli P, et al. Next-generation neoantigen mRNA vaccines: immuno-engineering strategies for precision cancer immunotherapy. Cell Oncol 2026;49(2).
[6] Wu DW, Jia SP, Xing SJ, Ma HL. Personalized neoantigen cancer vaccines: current progression, challenges and a bright future. Clin Exp Med 2024;24(1):229.
[7] Jefferies Equity Research. Moderna Initiation, 2026.
[8] Woods B, et al. Med Decis Making 2020;40:461-471.
[9] Keytruda.com. Cost Information, March 2026. GoodRx, May 2026.
[10] Schwarze K, et al. Eur J Hum Genet 2020;28:1322-1331.
[11] Bensimon AG, Zhou ZY, Jenkins M, et al. Cost-effectiveness of pembrolizumab for the adjuvant treatment of resected high-risk stage III melanoma in the United States. J Med Econ 2019;22(10):981-993.
[12] Johnston KM, et al. Pharmacoeconomics 2021;39:241-253.
[13] Zhang S, Bensimon AG, Xu R, et al. Cost-effectiveness analysis of pembrolizumab as an adjuvant treatment of resected stage IIB or IIC melanoma in the United States. Adv Ther 2023;40:3038-3055.
[14] Mclean A, van Hest N. An exploratory cost-effectiveness analysis of cancer vaccines in combination with current immune checkpoint inhibitors vs immune checkpoint inhibitor monotherapy: a case study for V940 in high-risk stage 3 melanoma in the US. Value Health 2024;27(12):S2 (EE354).

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
| LR monthly cost | $3,000 | Literature [16] |
| DM monthly cost | $12,000 | Literature [13] |
| Utility RF | 0.83 | Bensimon 2019 [14] |
| Utility LR | 0.64 | Bensimon 2019 [14] |
| Utility DM | 0.55 | Zhang 2023 [16] |
| AE disutility | 0.05 | Bensimon 2019 |
| Discount rate | 3% | Standard |
| Time horizon | 40 years | Lifetime |

### Table 2. Base Case Results
| | Pembrolizumab | Combo | Incremental |
|---|:---:|:---:|:---:|
| Total QALYs | 7.69 | 14.90 | 7.22 |
| Total Cost | $748,624 | $924,111 | $175,488 |
| ICER | — | — | **$24,317/QALY** |

### Table 3. Probabilistic Sensitivity Analysis
| Metric | Value |
|--------|-------|
| Mean ICER | $24,194 |
| Median ICER | $23,984 |
| 95% CI | $14,752 – $34,786 |
| P(CE) at $100K/QALY | 100% |
| P(CE) at $150K/QALY | 100% |

### Table 4. Scenario Analyses
| Scenario | ICER |
|:---|:---:|
| Base case (log-normal, 40y) | $24,317 |
| Intismeran $100,000 | $22,566 |
| Intismeran $300,000 | $27,781 |
| Intismeran $5,000,000 | $91,290 |
| 0% discount rate | $37,923 |
| 5% discount rate | $14,047 |
| Weibull OS | $38,065 |
| Treatment-effect waning (5→20y) | $75,206 |
| General population constraint | Combo dominant |
| 10-year horizon | Combo dominant |
| 20-year horizon | Combo dominant |

---

## Figures

- **Figure 1:** Log-normal model RFS and OS curves, with 5-year data points from KEYNOTE-942 → `output/fig1_survival.png`
- **Figure 2:** Cost-effectiveness acceptability curve → `output/fig2_ceac.png`
- **Figure 3:** Tornado diagram, one-way sensitivity analysis → `output/fig3_tornado.png`
- **Figure 4:** Cost-effectiveness plane (1,000 PSA iterations) → `output/fig4_ce_plane.png`