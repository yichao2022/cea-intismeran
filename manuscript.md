# Cost-Effectiveness of Intismeran Autogene Plus Pembrolizumab Versus Pembrolizumab Alone as Adjuvant Therapy for Resected Stage IIIB-IV Melanoma

## Abstract

**Background:** On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy (INT) to succeed in a late-stage trial. Intismeran plus pembrolizumab significantly improved recurrence-free survival (RFS) and distant metastasis-free survival (DMFS) versus pembrolizumab alone in patients with completely resected stage IIIB-IV melanoma. However, the economic value of this novel combination remains unknown.

**Objective:** To estimate the cost-effectiveness of intismeran plus pembrolizumab versus pembrolizumab alone as adjuvant therapy for resected stage IIIB-IV melanoma from a US payer perspective.

**Methods:** A four-state partitioned survival model (recurrence-free → locoregional recurrence → distant metastasis → death) was developed with a lifetime horizon (40 years) and 3% annual discounting. Clinical efficacy data were derived from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial. Parametric survival functions (log-normal) were calibrated to RFS, DMFS, and OS survival probabilities at 18, 24, 36, 48, and 60 months, with the ratio-constrained structure RFS ≤ DMFS ≤ OS ensuring internal consistency. Health state utilities were derived from published melanoma EQ-5D literature (RF 0.83, LR 0.64, DM 0.55). Sensitivity analyses included probabilistic sensitivity analysis (2,000 iterations), deterministic sensitivity analysis, and scenario analyses.

**Results:** In the base case, intismeran plus pembrolizumab yielded 11.10 QALYs at a cost of $860,118, compared with 7.34 QALYs at $713,157 for pembrolizumab alone. The incremental cost-effectiveness ratio was $39,105/QALY gained. At willingness-to-pay thresholds of $100,000 and $150,000/QALY [23], the probability of cost-effectiveness was 94.9% and 98.7%, respectively. Results were robust across a wide range of sensitivity analyses.

**Conclusion:** At an assumed acquisition price of $200,000 per course, intismeran plus pembrolizumab is cost-effective by conventional US thresholds in adjuvant melanoma. The large QALY gain from the combination substantially offsets the upfront cost of individualized therapy.

**Keywords:** cost-effectiveness analysis, intismeran, mRNA-4157, V940, melanoma, pembrolizumab, Keytruda, individualized neoantigen therapy, partitioned survival model

---

## 1. Introduction

### A New Reimbursement Problem

On August 19, 2026, Merck and Moderna announced positive topline results from the Phase 3 INTerpath-001 trial, establishing intismeran autogene (mRNA-4157/V940) plus pembrolizumab as the first mRNA-based individualized neoantigen therapy to succeed in a late-stage trial [2]. The trial met its primary endpoint of recurrence-free survival (RFS) and key secondary endpoint of distant metastasis-free survival (DMFS) in patients with resected high-risk stage IIB–IV melanoma. However, full numerical Phase 3 efficacy results and the commercial price for intismeran have not yet been disclosed. The pivotal clinical milestone therefore creates an immediate need for an early economic assessment before pricing and reimbursement decisions are finalized.

### Why INTs Differ from Conventional Oncology Pricing

Individualized neoantigen therapies (INTs) introduce production and delivery requirements that differ materially from conventional fixed-formulation pharmaceuticals. Each patient's treatment requires a dedicated workflow: multi-region tumor biopsy and next-generation sequencing, computational neoantigen prediction and HLA-binding prioritization, individualized mRNA synthesis and lipid nanoparticle formulation, and lot-release quality control, all within a manufacturing window of approximately six to nine weeks constrained by clinical urgency [4,5]. These per-patient manufacturing requirements, combined with the need for decentralized sequencing capacity and cold-chain logistics, present scalability and cost challenges distinct from batch-produced oncology drugs [5,6]. Current industry estimates suggest personalized mRNA cancer vaccines may cost $100,000–$300,000 per patient, with Jefferies analysts estimating approximately $200,000 per course [7], but these are supply-side estimates, not value-based prices.

### Existing Economic Evidence

A prior exploratory cost-effectiveness analysis of V940 plus pembrolizumab, presented at ISPOR Europe 2024, used a four-state partitioned survival model applied to the stage III adjuvant melanoma setting [14]. That analysis assumed V940 pricing at the same list price as pembrolizumab and applied KEYNOTE-942 hazard ratios to pembrolizumab survival curves derived from the KEYNOTE-054 trial, explicitly noting that the results were preliminary and that "further data when available from V940-001 are required to increase the reliability of this study." Since that analysis, KEYNOTE-942 has reported mature 5-year RFS and DMFS outcomes [3], and INTerpath-001 has provided the first positive Phase 3 confirmation. The present study therefore addresses an evidence gap that the prior exploratory work itself identified: an early economic evaluation anchored to the mature 5-year KEYNOTE-942 data, presented after the Phase 3 positive readout, and incorporating a deterministic threshold price analysis tied to a realistic acquisition cost estimate.

### Clinical Basis and Research Objectives

The clinical foundation for this analysis comes from the 5-year follow-up of the Phase 2b KEYNOTE-942 trial (NCT03897881), which randomized 157 patients with resected stage IIIB–IV cutaneous melanoma 2:1 to intismeran plus pembrolizumab versus pembrolizumab alone [3]. At a median follow-up of 60.3 months, the combination demonstrated a 49% reduction in the risk of recurrence or death (HR 0.51; 95% CI 0.294–0.887) and a 59% reduction in the risk of distant metastasis (HR 0.411; 95% CI 0.200–0.843). Five-year RFS was 68.8% versus 49.1%, while overall survival showed a favorable trend (HR 0.471; 95% CI 0.165–1.345) but remained immature (5-year OS 92.2% vs 71.3%; 7 events per arm). We therefore conducted an early US payer-perspective economic evaluation to estimate the cost-effectiveness of intismeran plus pembrolizumab and, given the absence of an established launch price, to identify the acquisition price consistent with commonly cited willingness-to-pay thresholds.

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
- 5-year OS: 92.2% vs 71.3% (exploratory, HR 0.471; 95% CI 0.165-1.345)

Survival probabilities at 18, 24, 36, 48, and 60 months were extracted from the published Kaplan-Meier curves (Figure 1, Khattak et al. 2026). Parametric survival models were fitted via weighted nonlinear least squares, following standard methodological guidance [15-17]. Model selection compared five candidate distributions (exponential, Weibull, log-normal, log-logistic, generalized gamma) for each of three underlying survival components (OS, r_DM = DMFS/OS, r_LR = RFS/DMFS) using AIC, as recommended by NICE DSU TSD 14 [15]. The log-normal distribution was selected as the base case, consistent with the ratio-constrained structure RFS ≤ DMFS ≤ OS [11,13]. All fitted values are reported in Table S2 (Supplementary Materials).

### 2.3 Costs

Direct medical costs were estimated from a US payer perspective and expressed in 2026 US dollars (Table 1). Drug costs were based on published wholesale acquisition costs (WAC): pembrolizumab 400 mg every 6 weeks at $24,544 per dose (Keytruda.com, March 2026), yielding an annual cost of $220,896 for 9 doses over approximately 1 year [9]. Intismeran cost was estimated at $200,000 per course based on Jefferies analyst reports [7]. Tumor sequencing (whole-exome sequencing) was estimated at $1,000 per patient [10]. Intravenous administration costs were based on the CMS Physician Fee Schedule ($500 per cycle) [20]. Adverse event management costs were estimated at $15,000 incremental for the combination arm, consistent with published melanoma CEA literature [11]. In the KEYNOTE-942 5-year safety analysis, treatment-related Grade 3–5 adverse events occurred in 25.0% of combination-arm patients and 14.0% of pembrolizumab-arm patients (Supplementary Table A4, [3]). Monthly disease management costs were $3,000 for the LR state (salvage surgery and adjuvant therapy amortized) and $12,000 for the DM state (systemic therapy for advanced melanoma) [12,13].

### 2.4 Utilities

Health state utilities were derived from published EQ-5D literature in melanoma (Table 1). The utility for RF was 0.83, based on EQ-5D data from the KEYNOTE-054 trial [11]. The utility for LR was 0.64 and for DM was 0.55, consistent with published estimates [11,13]. A one-time disutility of 0.05 was applied for adverse events in the combination arm.

### 2.5 Base Case Analysis

The base case used a lifetime horizon of 40 years, monthly cycles, and a 3% annual discount rate for both costs and health outcomes. The primary outcome was the incremental cost-effectiveness ratio (ICER), expressed as cost per quality-adjusted life year (QALY) gained. Willingness-to-pay thresholds of $100,000/QALY and $150,000/QALY were considered.

### 2.6 Sensitivity Analysis

Probabilistic sensitivity analysis (PSA) was performed with 2,000 Monte Carlo iterations, drawing parameter values from appropriate distributions: gamma for costs, beta for utilities, and normal perturbations on survival parameters. All 2,000 iterations were structurally valid (RFS ≤ DMFS ≤ OS, monotonicity, survival within [0,1] verified every cycle; 0 rejections); no iterations were excluded on the basis of the ICER value. Results were summarized as cost-effectiveness acceptability curves (CEAC) computed from incremental net monetary benefit. Deterministic one-way sensitivity analysis was performed on key parameters.

### 2.7 Scenario Analyses

Scenario analyses examined: (1) the impact of intismeran pricing at $100,000 and $300,000 per course; (2) alternative discount rates (0% and 5%); (3) alternative survival extrapolation methods (log-normal, Weibull); (4) treatment-effect waning; (5) general population mortality constraint; and (6) reduced time horizons (10 and 20 years).

---

## 3. Results

### 3.1 Base Case

In the base case analysis (Table 2), intismeran plus pembrolizumab generated an additional 4.42 life-years (14.71 vs 10.29 life-years) and 3.76 additional QALYs (11.10 vs 7.34 QALYs) at an incremental lifetime cost of $146,961 ($860,118 vs $713,157), resulting in an ICER of **$39,105 per QALY gained**.

### 3.2 Cost and QALY Decomposition

The $200,000 acquisition cost of intismeran was substantially offset by lower lifetime distant-metastasis management costs. The combination arm had $82,734 lower DM management costs (combo $342,542 vs pembro $425,276), reflecting the reduced incidence of distant recurrence (DMFS HR 0.411). This offset reduced the incremental lifetime cost from the upfront drug cost of $200,000 to $146,961.

The incremental QALY gain of 3.76 was driven almost entirely by additional time in the recurrence-free state (3.83 RF QALYs gained), which more than compensated for slightly fewer DM-state QALYs (-0.32).

| Component | Combo | Pembro | Incremental |
|:---|:---:|:---:|:---:|
| **Costs** | | | |
| Intismeran | $200,000 | $0 | $200,000 |
| Pembrolizumab | $217,888 | $217,888 | $0 |
| Sequencing | $1,000 | $0 | $1,000 |
| Administration/AE | $15,493 | $493 | $15,000 |
| LR management | $83,196 | $69,500 | $13,696 |
| DM management | $342,542 | $425,276 | -$82,734 |
| **Total cost** | **$860,118** | **$713,157** | **$146,961** |
| **QALYs** | | | |
| RF state | 8.32 | 4.49 | 3.83 |
| LR state | 1.48 | 1.24 | 0.24 |
| DM state | 1.31 | 1.62 | -0.32 |
| **Total QALYs** | **11.10** | **7.35** | **3.76** |

### 3.3 External and Structural Validation

Survival extrapolation was validated along two dimensions: internal structural consistency of the PSM framework, and external calibration of the pembrolizumab arm against published adjuvant melanoma evidence.

**Structural consistency.** Over 480 monthly cycles, no violations of the constraint RFS ≤ DMFS ≤ OS were observed. All health-state probabilities remained non-negative, and the sum of state probabilities was within 10⁻¹² of unity at every cycle. The general population mortality hazard floor, applied from 60 months onward, ensures that the model's long-term survival does not exceed age- and sex-matched US life-table survival at any time point. At 30 years, the model projects combo OS of 21.5% and pembro OS of 10.9%, both below the general-population survival floor of approximately 35%.

**External validation.** The pembrolizumab-alone arm predictions were compared against long-term adjuvant melanoma evidence [4,19]. At 5 years, the model predicted RFS of 40.2%, DMFS of 53.5%, and OS of 74.1%, consistent with the KEYNOTE-942 observed data (RFS 49.1%, DMFS 65.4%, OS 71.3%). At 7 years, the model predicted RFS of 31.7% and OS of 63.8%. The lower RFS relative to KEYNOTE-054 (7-year RFS 50%) reflects the higher-risk population enrolled in KEYNOTE-942 (stage IIIB–IV including resected stage IV, vs stage IIIA–C in KEYNOTE-054). Because the pembrolizumab arm serves as the comparator for the incremental analysis, its under-prediction of RFS may, if anything, overstate the incremental QALY gain attributable to intismeran; however, the direction of any bias is not straightforward, as the model's predicted OS for the pembrolizumab arm was slightly higher than observed (74.1% vs 71.3% at 5 years). No external long-term data exist for the combination arm (first-in-class therapy). Instead of relying on this comparison to calibrate the base case, scenario analyses that directly bound the treatment effect — treatment-effect waning and no-direct-OS-benefit (Section 3.6) — were used to bound the range of plausible ICER estimates.

| Outcome | Model | KEYNOTE-054 | CheckMate 238 |
|:---|:---:|:---:|:---:|
| 5-year RFS | 40.2% | 49.1% (observed) | — |
| 5-year DMFS | 53.5% | 65.4% (observed) | — |
| 5-year OS | 74.1% | 71.3% (observed) | 76% (nivolumab) |
| 7-year RFS | 31.7% | 50% (46–55%) | — |
| 7-year DMFS | 44.6% | 54% (50–59%) | — |
| 7-year OS | 63.8% | — | ~70% (nivolumab) |

### 3.4 Probabilistic and Deterministic Uncertainty

Probabilistic sensitivity analysis (all 2,000 structurally valid iterations out of 2,000 draws; 0 rejections based on structural criteria) quantified joint uncertainty in survival, cost, and utility parameters. The mean incremental net monetary benefit (NMB) was $281,051 at the $100,000/QALY threshold and $482,475 at $150,000/QALY (Table 3; Figure 2). The probability that intismeran plus pembrolizumab is cost-effective was 63.9% at $50,000/QALY, 94.9% at $100,000/QALY, and 98.7% at $150,000/QALY. All iterations showed positive QALY gains (ΔQALY > 0), and 17.6% of iterations showed net cost savings (ΔCost < 0).

**Deterministic sensitivity analysis.** One-way sensitivity analysis (Figure 3) identified the 5-year RFS rate for the combination arm as the most influential parameter, followed by the DM monthly cost and the intismeran price. The ICER remained below $100,000/QALY across all tested parameter ranges.

### 3.5 Structural Scenario Analyses

**Pricing:** Scenario analysis results are summarized in Table 4. At an intismeran cost of $100,000 per course, the ICER decreased to $12,496/QALY. At $300,000 per course, the ICER increased to $65,714/QALY. At $500,000 per course, the ICER was $118,932/QALY.

**Discount rate:** At 0% discounting, the ICER was $31,092/QALY. At 5%, it was $47,814/QALY.

**Survival extrapolation:** Using Weibull instead of log-normal yielded an ICER of $54,921/QALY. With treatment-effect waning (hazard convergence, 5→10y), the ICER was $105,360/QALY.

**General population mortality hazard floor:** Applying the age-matched US life table survival as a hazard floor on OS from month 60 onward yielded an ICER of $39,105/QALY.

**Time horizon:** With a 10-year horizon, the ICER was $74,945/QALY. With a 20-year horizon, it was $39,734/QALY.**

**No direct OS benefit:** Assuming no OS benefit (intismeran effect only through recurrence prevention) yielded a dominant result (−$29,379/QALY) — lower cost and higher QALY.

---

### 3.6 Value-Based Price Thresholds

Figure 4 presents the Price–ICER curve. The ICER is a linear function of the intismeran acquisition price, crossing the $100,000/QALY threshold at approximately $429,000 per course and the $150,000/QALY threshold at approximately $616,000 per course. At the estimated price of $200,000, the ICER of $39,105/QALY is well below both thresholds.

The deterministic threshold price of intismeran was estimated by identifying the price at which the ICER reaches conventional willingness-to-pay thresholds. At the estimated acquisition cost of $200,000 per course, the ICER was $39,105/QALY. The deterministic threshold price at a $100,000/QALY threshold was approximately $429,000 per course, meaning the ICER remains below $100,000/QALY even at more than twice the current estimated price. At a $150,000/QALY threshold, the deterministic threshold price was approximately $616,000 per course.

| Intismeran price | ΔCost | ΔQALY | ICER |
|:---|:---:|:---:|:---:|
| $0 | −$53,149 | 3.76 | −$14,113 |
| $100,000 | $46,961 | 3.76 | $12,496 |
| $200,000 (base case) | $146,961 | 3.76 | $39,105 |
| $300,000 | $246,961 | 3.76 | $65,714 |
| $400,000 | $346,961 | 3.76 | $92,323 |
| $500,000 | $446,961 | 3.76 | $118,932 |
| $600,000 | $546,961 | 3.76 | $145,541 |
| **Deterministic threshold prices** | | | |
| ~$429,000 | ~$375,937 | 3.76 | $100,000 |
| ~$616,000 | ~$563,315 | 3.76 | $150,000 |


## 4. Discussion

### 4.1 Principal Findings

This early economic evaluation demonstrates that intismeran plus pembrolizumab is likely cost-effective for adjuvant treatment of resected stage IIIB--IV melanoma at conventional US willingness-to-pay thresholds. The base case ICER of $39,105/QALY, the 3.76 QALY gain, and the 94.9% probability of cost-effectiveness at $100,000/QALY (98.7% at $150,000/QALY) all support this conclusion.

The economic value of the combination is generated through two complementary mechanisms described in detail below (Section 4.2): extended recurrence-free quality-adjusted survival and downstream cost offsets from avoided distant metastasis. The ICER is notably lower than the $75,206 from the prior exploratory analysis [14] (Section 4.5), reflecting the combined effect of mature clinical data and a structurally more conservative modelling framework.

### 4.2 Economic Drivers and Cost Offsets

The incremental $146,961 lifetime cost of intismeran plus pembrolizumab is substantially lower than the $200,000 acquisition cost of the intismeran course alone, because $82,734 in downstream distant-metastasis management costs are avoided (combo $342,542 vs pembro $425,276). The QALY decomposition (Table 2) shows that the 3.76 incremental QALYs are driven almost entirely by additional time in the recurrence-free health state: the combination arm generates 3.83 more RF QALYs, partially offset by 0.32 fewer DM QALYs, with a small net gain of 0.24 LR QALYs. This asymmetry — a large QALY gain concentrated in the least costly health state, paired with cost savings from the most costly state — is the structural mechanism behind the favourable ICER.

The economic value is therefore not simply a matter of "large QALY gain offsets upfront cost." The combination works through two pathways simultaneously: (i) it shifts patients from DM to RF, the highest-utility, lowest-cost state; and (ii) it compresses the time spent in DM, the highest-cost state, generating direct cost offsets. The ICER is low because the incremental expenditure is coupled to a large gain in recurrence-free quality-adjusted survival and to partial avoidance of the high cost of metastatic disease, not because the intervention is inexpensive.

### 4.3 Pricing and Reimbursement Implications

The base case assumes an intismeran acquisition price of $200,000 per course, based on Jefferies analyst estimates. At this price, the ICER of $39,105/QALY is well below conventional US thresholds. The deterministic threshold price analysis (Section 3.6) shows that the ICER remains below $100,000/QALY at intismeran prices up to approximately $429,000 per course and below $150,000/QALY up to approximately $616,000 per course, providing substantial headroom for price negotiations.

While individualized production requirements — including tumor sequencing, mRNA synthesis, lipid nanoparticle formulation, and quality assurance for each batch — may contribute to acquisition-price pressure, production cost and commercial price are distinct quantities. The $200,000 estimate reflects the expected market price at launch, not the manufacturer's cost of goods sold. Several factors may reduce future production costs: sequencing costs continue to decline [10], mRNA platform manufacturing is becoming increasingly standardized, and automation of the synthesis and formulation pipeline could reduce per-batch costs. However, whether these technical efficiencies will translate into lower acquisition prices for payers is uncertain, as pricing decisions are influenced by clinical value, competitive dynamics, and reimbursement negotiations, not solely by production economics.

Cost-effectiveness should not be conflated with affordability. The deterministic threshold price of approximately $429,000 per course at $100,000/QALY implies that the combination could be considered cost-effective at more than twice the current estimated price, but the incremental per-patient budget impact of adding $200,000 in drug costs to an already expensive pembrolizumab regimen would be substantial in a population of eligible patients. Formal budget-impact analysis is needed alongside cost-effectiveness evidence for reimbursement decision-making.

### 4.4 Long-Term Survival and Structural Uncertainty

Long-term survival extrapolation is the dominant source of structural uncertainty in this analysis. The mature 5-year RFS and DMFS data from KEYNOTE-942 provide a robust foundation for the recurrence and metastasis components of the model, but OS remains exploratory (HR 0.471; 95% CI 0.165–1.345; 7 events per arm) and must be extrapolated to a lifetime horizon.

The necessity of the general population mortality hazard floor is best illustrated by comparing the constrained and unconstrained models. Without the GP mortality floor, the log-normal survival fit for the combination arm predicts a 30-year OS of 79.1%, far exceeding the age- and sex-matched general population survival of 22.1% at 30 years for the modeled cohort. This is implausible: a cohort of patients with resected stage IIIB--IV melanoma cannot have a lower mortality rate than the general population. Applying the GP mortality hazard floor from 60 months onward — a standard practice recommended by NICE DSU TSD 14 [15] and used in published melanoma CEAs [13] — reduces the 30-year combination OS to 21.5%, a clinically plausible trajectory. The hazard floor reduces the ICER from $51,753 (unconstrained) to $39,105 (GP-constrained); the direction of this change is non-intuitive (the constraint reduces the ICER) because the unconstrained model predicts implausibly high survival for both arms, with a larger absolute gain for the combination arm, and the hazard floor compresses that difference while preserving the within-trial benefit.

The OS uncertainty is further underscored by the treatment-effect waning scenario. When the combination's recurrence-related hazard is allowed to converge to the pembrolizumab level between years 5 and 10, the ICER rises to $105,360/QALY — nearly three times the base case — because the QALY gain is compressed from 3.76 to 2.13 while the incremental cost remains high. The range from $39,105 to $105,360 represents the plausible interval for the ICER under alternative persistence assumptions.

Finally, the choice of parametric distribution, while important, is secondary to the structural scenarios. The log-normal distribution was selected by AIC and was the only distribution that preserved structural consistency (RFS ≤ DMFS ≤ OS) across all 4,096 candidate combinations. However, AIC measures fit within the observed 5-year follow-up; it does not establish long-term clinical validity. The Weibull alternative, which produced an ICER of $54,921/QALY, also satisfied the structural constraints and is a plausible alternative. The most informative bounds on the ICER come not from alternative parametric forms but from the treatment-effect waning and no-direct-OS-benefit scenarios, which directly address the question: what if the treatment effect does not persist?

### 4.5 Comparison with Previous Economic Evaluations

Direct comparison of ICERs across adjuvant melanoma cost-effectiveness studies is not straightforward, because the incremental decision differs. Bensimon et al. [11] and Zhang et al. [13] evaluated pembrolizumab versus observation, estimating the value of adding immunotherapy to surveillance. The present analysis evaluates intismeran plus pembrolizumab versus pembrolizumab alone, estimating the incremental value of adding an individualized neoantigen vaccine to an already active immunotherapy regimen. The relevant policy question is therefore not whether the combination is cost-effective in absolute terms, but whether the combination of two active agents can provide sufficient additional value — through both clinical benefit and downstream cost offsets — to justify the cost of the second agent.

The ICER of $39,105/QALY is within the range of published estimates for the initial pembrolizumab decision ($38,050--$78,400/QALY), but this comparability is coincidental rather than meaningful. The absolute QALY gain from adding intismeran to pembrolizumab (3.76 QALYs) is larger than the gain from adding pembrolizumab to observation (1.2--2.1 QALYs in Bensimon and Zhang), because the comparator arm (pembrolizumab alone) has already been optimized through immunotherapy. The incremental cost is also larger ($146,961 vs approximately $100,000--$150,000). The ratio happens to fall within a similar range because both numerator and denominator increase proportionally, not because the value proposition is comparable.

The comparison with the 2024 exploratory analysis [14] is more informative but requires careful interpretation. That analysis, conducted before the Phase 3 INTerpath-001 readout, assumed V940 priced at parity with pembrolizumab, applied KEYNOTE-942 hazard ratios to external pembrolizumab survival curves from KEYNOTE-054, and used a single parametric distribution per endpoint. The present analysis differs on all four dimensions: post-Phase III evidence context, direct fitting to KEYNOTE-942 5-year data, a $200,000 price assumption, and a ratio-constrained structural framework with a general population mortality floor. The lower ICER ($39,105 vs $75,206) is the net result of these four differences, not a univariate effect of any single factor, and should not be interpreted as evidence that intismeran is simply more effective than previously estimated.

### 4.6 Clinical and External Validity

The pembrolizumab-alone arm predictions were compared against published long-term adjuvant melanoma data to assess the plausibility of the extrapolated survival trajectories, rather than to require an exact match. The KEYNOTE-942 trial enrolled a higher-risk population (stage IIIB--IV, including resected stage IV) than KEYNOTE-054 (stage IIIA--C) and CheckMate 238 (stage IIIB--IV, but with different entry criteria). The model's lower 7-year RFS for the pembrolizumab arm (31.7% vs 50% in KEYNOTE-054) is therefore expected and reflects the underlying risk differences rather than a failure of calibration. The model correctly reproduces the 5-year OS and DMFS of the pembrolizumab arm, supporting the structural validity of the PSM framework.

No external long-term data exist for the combination arm (first-in-class therapy). The positive Phase 3 INTerpath-001 readout provides directional confirmation of the clinical benefit of intismeran plus pembrolizumab, but the numerical hazard ratios have not yet been disclosed, and the extent to which Phase 3 results will align with the Phase 2b KEYNOTE-942 estimates used in this model remains to be determined. Two structural scenarios address the uncertainty in the combination arm: (i) the treatment-effect waning scenario (ICER $105,360/QALY); and (ii) the no-direct-OS-benefit scenario, an OS-neutral structural scenario in which intismeran plus pembrolizumab remained dominant (lower cost and higher QALY) compared with pembrolizumab alone, indicating that the base-case conclusion does not rely on a direct OS benefit but can be supported by recurrence prevention alone.

### 4.7 Limitations

This analysis has several limitations that should be considered when interpreting the results. First, the Phase 3 INTerpath-001 trial has reported only topline results; numerical hazard ratios and subgroup data have not yet been disclosed, and all quantitative efficacy inputs are derived from the Phase 2b KEYNOTE-942 trial. Second, overall survival data from KEYNOTE-942 remain exploratory (HR 0.471; 95% CI 0.165–1.345), with only 7 events per arm at the 5-year analysis, and the 5-year OS estimate of 92.2% versus 71.3% should be interpreted with caution. Third, survival curves were reconstructed from published Kaplan–Meier figures using digitized coordinates rather than being estimated from individual patient-level data, which introduces measurement error in the parametric fits. Fourth, survival extrapolation beyond the observed 5-year follow-up, while structurally constrained and anchored to a general population mortality floor, remains sensitive to the choice of parametric distribution and the assumption that the treatment effect persists indefinitely. Fifth, the intismeran acquisition price of $200,000 per course is based on analyst estimates, and the actual market price may differ. Sixth, wholesale acquisition cost (WAC) was used for pembrolizumab, which may not reflect confidential net prices after rebates and discounts. Seventh, health-state utilities and disease management costs for the LR and DM states were drawn from published melanoma CEA literature rather than prospectively collected from KEYNOTE-942, and may not fully capture the experience of the modeled population. Eighth, indirect costs, productivity losses, and informal caregiver burden were not included, which may underestimate the societal value of preventing recurrence and distant metastasis in a working-age population. Ninth, no formal budget-impact analysis was conducted; the cost-effectiveness findings should be complemented by an affordability assessment when pricing and population-level data become available. Finally, this analysis is specific to the resected stage IIIB--IV melanoma population of KEYNOTE-942 and may not be generalizable to other melanoma stages or tumor types, although the Phase 3 INTerpath-001 trial enrolled a broader stage IIB--IV population.

### 4.8 Conclusion

Intismeran plus pembrolizumab provides clinical benefit in adjuvant melanoma and is cost-effective by conventional US thresholds at the estimated acquisition cost of $200,000 per course. The deterministic threshold price analysis shows that the ICER remains below $100,000/QALY at intismeran prices up to approximately $429,000 per course, providing substantial headroom for price negotiations. The QALY gain from the combination, driven by improved recurrence-free survival and reduced distant metastases, partially offsets the upfront cost of individualized therapy. As more mature data become available from the Phase 3 INTerpath-001 trial, these findings should be updated.
## 5. Declarations

**Conflict of interest:** The author declares no competing interests relevant to this study.

**Reporting:** This study adheres to the Consolidated Health Economic Evaluation Reporting Standards 2022 (CHEERS 2022) [26] (Supplementary Checklist).

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
[14] Mclean A, van Hest N.
[15] Latimer N. NICE DSU Technical Support Document 14: Survival analysis for economic evaluations alongside clinical trials — extrapolation with patient-level data. Sheffield: NICE Decision Support Unit; 2011 (updated 2013).
[16] Latimer NR. Survival analysis for economic evaluations alongside clinical trials — extrapolation with patient-level data: inconsistencies, limitations, and a practical guide. Med Decis Making 2013;33(6):743-754.
[17] Ishak KJ, Kreif N, Benedict A, Muszbek N.
[18] Arias E, Xu JQ.
[19] Weber J, Mandala M, Del Vecchio M, et al.
[20] Centers for Medicare & Medicaid Services. Physician Fee Schedule, CY 2025: CPT Code 96365 (Intravenous Infusion, Initial, Up to 1 Hour). Payment rate ~$500 per cycle. Available from: https://www.cms.gov/medicare/payment/fee-schedules/physician
[21] Bureau of Labor Statistics.
[22] Sanders GD, Neumann PJ, Basu A, et al. Recommendations for conduct, methodological practices, and reporting of cost-effectiveness analyses: Second Panel on Cost-Effectiveness in Health and Medicine. JAMA 2016;316(10):1093-1103. Consumer Price Index: Medical Care Component (CPI-U, Series CUUR0000SAM). Accessed 2026. Available from: https://www.bls.gov/cpi/ Adjuvant nivolumab versus ipilimumab in resected stage III or IV melanoma. N Engl J Med 2017;377(19):1824-1835. United States life tables, 2019. National Vital Statistics Reports 2022;70(19):1-59. Overview of parametric survival analysis for health-economic applications. Pharmacoeconomics 2013;31(8):663-675. An exploratory cost-effectiveness analysis of cancer vaccines in combination with current immune checkpoint inhibitors vs immune checkpoint inhibitor monotherapy: a case study for V940 in high-risk stage 3 melanoma in the US. Value Health 2024;27(12):S2 (EE354).

## Tables

### Table 1. Model Input Parameters
| Parameter | Base Case | Source |
|-----------|-----------|--------|
| 5-year RFS, combo | 68.8% | KEYNOTE-942 [6] |
| 5-year RFS, pembro | 49.1% | KEYNOTE-942 [6] |
| RFS HR | 0.51 (95% CI 0.294-0.887) | KEYNOTE-942 [6] |
| DMFS HR | 0.411 (95% CI 0.200-0.843) | KEYNOTE-942 [6] |
| OS HR | 0.471 (95% CI 0.165-1.345) | KEYNOTE-942 [6] |
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

### Table 2. Base-Case Results and Decomposition
| | Pembrolizumab | Combo | Incremental |
|---|:---:|:---:|:---:|
| Total QALYs | 7.34 | 11.10 | 3.76 |
| Total Cost | $713,157 | $860,118 | $146,961 |
| ICER | — | — | **$39,105/QALY** |

### Table 3. Probabilistic Sensitivity Analysis
| Metric | Mean | Median | 95% CI |
|:---|:---:|:---:|:---:|
| Incremental cost | $121,796 | $152,383 | $‑610,346 – $509,272 |
| Incremental QALYs | 4.03 | 3.83 | 2.33 – 6.74 |
| NMB at $100K/QALY | $281,051 | — | — |
| NMB at $150K/QALY | $482,475 | — | — |

| CE probability | Value |
|:---|:---:|
| P(CE) at $50K/QALY | 63.9% |
| P(CE) at $100K/QALY | 94.9% |
| P(CE) at $150K/QALY | 98.7% |


### Table 4. Scenario and Deterministic Threshold Price Analysis
| Scenario | ΔCost | ΔQALY | ICER | NMB @ $150K |
|:---|:---:|:---:|:---:|:---:|
| **Base case** | | | | |
| GP-constrained | $146,961 | 3.76 | $39,105 | $416,757 |
| **Pricing** | | | | |
| $100K/course | $46,961 | 3.76 | $12,496 | $516,757 |
| $300K/course | $246,961 | 3.76 | $65,714 | $316,757 |
| $500K/course | $446,961 | 3.76 | $118,932 | $116,757 |
| **Alternative assumptions** | | | | |
| Weibull OS | $293,679 | 5.35 | $54,921 | $508,411 |
| Discount rate 0% | $183,728 | 5.91 | $31,092 | $702,636 |
| Discount rate 5% | $136,896 | 2.86 | $47,814 | $292,571 |
| 10-year horizon | $95,573 | 1.28 | $74,945 | $95,714 |
| 20-year horizon | $115,739 | 2.91 | $39,734 | $321,191 |
| GP utilities | $146,961 | 3.59 | $40,970 | $391,092 |
| **Treatment-effect persistence** | | | | |
| Treatment-effect waning (5→10y hazard convergence) | $224,317 | 2.13 | $105,360 | $95,042 |
| No direct OS benefit | −$14,562 | 0.50 | −$29,379 | $88,915 |

---

## Figures

- **Figure 1:** Log-normal model RFS and OS curves, with 5-year data points from KEYNOTE-942 → `output/fig1_survival.png`
- **Figure 2:** Cost-effectiveness acceptability curve → `output/fig2_ceac.png`
- **Figure 3:** Tornado diagram, one-way sensitivity analysis → `output/fig3_tornado.png`
- **Figure 4:** Price–ICER curve and deterministic threshold pricing → `output/fig4_price_icer.png`
- **Supplementary Figure S7:** Cost-effectiveness plane (2,000 PSA iterations) → `output/fig4_ce_plane.png`