# Peer Review Report — Cost-Effectiveness of Intismeran Autogene Plus Pembrolizumab in Adjuvant Melanoma

**Target Journal:** PharmacoEconomics  
**Manuscript ID:** CEINTIS-2026  

---

## Phase 0: Reviewer Configuration Card

| Role | Persona | Expertise | Focus |
|------|---------|-----------|-------|
| EIC | PharmacoEconomics Editor, prof. in HTA methodology | Health economics publishing, CHEERS reporting, journal scope fit | Novelty, journal fit, editorial framing |
| R1 (Methodology) | Health economist, PSM/survival modeling specialist | Partitioned survival models, PSA/DSA, parametric extrapolation | Model structure, fitting, PSA design, numerical consistency |
| R2 (Domain) | Melanoma oncologist-immunotherapy trialist | KEYNOTE-942/054, CheckMate 238, INTerpath-001 | Clinical data accuracy, external validity, trial representation |
| R3 (Perspective) | HTA/payer policy researcher | US payer economics, budget impact, coverage decisions | Policy relevance, affordability, price assumptions |
| DA (Devil's Advocate) | Adversarial methodology critic | Logical fallacy detection, cherry-picking identification | Core argument attack, strongest counter-arguments |

---

## EIC Report (Editor-in-Chief)

**总体评价：**  
本文在 INTerpath-001 三期读出后数周、intismeran 定价尚未公布的时间窗口,对其开展了早期成本效果评价,选题时点与阈值价格分析($527K/$852K)对 PharmacoEconomics 的 HTA 与支付方读者高度契合。模型结构常规但报告透明度高:GitHub 代码、4096 组合生存模型选择、外部验证、CHEERS 2022 清单齐备。讨论对不确定性（OS 不成熟、waning 情景、预算影响）的态度总体审慎。

**评分：**
- 期刊契合度: **8/10** — 早期经济学评价 + Phase 3 背景 = 本刊核心定位
- 原创性: **8/10** — 首篇 intismeran+ pembrolizumab CEA，阈值价格分析有新意
- 重要性: **7/10** — INT 个体化疫苗的定价问题具有高度时效性

**Issues:**

| # | 严重度 | 位置 | 问题 | 建议 |
|---|--------|------|------|------|
| E1 | MAJOR | 标题 | 标题写 "Stage IIIB–IV"，但摘要背景写 "Stage IIB–IV"，人群分期不一致 | 统一为 Stage IIIB–IV（KEYNOTE-942 纳入标准） |
| E2 | MAJOR | 摘要 Results | "Results were robust across scenarios" 但 hazard-convergence waning ICER = $104,247 超 $100K 阈值，措辞过度 | 改为 "Results were largely robust; ICER exceeded $100,000/QALY only under hazard-convergence waning" |
| E3 | MAJOR | Discussion §4.2 与 Results §3.1-3.2 | 成本/QALY 分解讨论（§4.2 第 2 段）与 Results §3.2 内容高度重复 | 删除 §4.2 第 2 段重复内容，或仅引用 "as described in Section 3.2" |
| E4 | MINOR | 标题 | "An Updated Early Economic Evaluation" 中 "Updated" 不够精确（此前只有 ISPOR 摘要，非正式发表论文） | 改为 "An Early Economic Evaluation" 或保留但加注说明 |
| E5 | MINOR | 摘要 Methods | 未提及 Phase 3 INTerpath-001 在方法学部分的衔接 | 在 Methods 末尾加一句说明 Phase 3 数据将更新模型 |

**初步建议：Minor Revision**  
选题时效性好、报告透明度高、模型方法学扎实。主要问题在措辞一致性（标题/摘要分期矛盾、"robust" 过度概括）和结果/讨论重复，属于可快速修正的表述问题。

---

## Reviewer 1 Report — Methodology

**总体评价：**  
该文采用 ratio-constrained partitioned survival model 结合 GP mortality hazard floor 方法,方法学选择合理。4096 种分布组合的系统筛选、基于 NLS Jacobian 的 PSA 协方差设计、0 rejections 的结构验证,均体现了较高的方法学严谨性。代码开源可复现。主要方法学问题在于:GP hazard floor 从 month 60 开始切换,但缺乏对切换点敏感性的正式分析;PSA 中 intismeran 价格固定($200,000)而非作为不确定性参数处理,这与 DSA 中价格是最敏感参数的事实不一致;Table S8 中 OS mu combo DSA 范围(4.228–4.628)与 base case(8.268)差距极大,疑似 typo。

**评分：**
- 方法严谨性: **8/10**
- 可复现性: **9/10** (GitHub 代码 + CSV 结果)
- 内部一致性: **7/10** (Table S8 参数范围疑似错误)

**Issues:**

| # | 严重度 | 位置 | 问题 | 建议 |
|---|--------|------|------|------|
| M1 | CRITICAL | Table S8, line 354 | OS mu combo base = 8.268, DSA low = 4.228, DSA high = 4.628 — low/high 范围仅 0.4 宽且均远低于 base,疑似复制粘贴错误。DSA 结果中 OS mu combo 的 lo/hi ICER 应基于 base±0.4 而非 base=8.396,与正文不一致 | 核实 DSA range 是否应为 8.068–8.528 (±0.2) 或其他合理范围;更新 Table S8 |
| M2 | MAJOR | Methods §2.6, PSA | Intismeran 价格在 PSA 中固定为 $200,000,但 DSA 显示价格是最敏感参数(范围 $43,702–$55,979)。PSA 不纳入最敏感参数的不确定性会低估 ICER 的总体变异 | 将 intismeran 价格纳入 PSA,从合理分布(如 Gamma 或 Uniform $160K–$240K)抽样;或明确说明为何固定(如"价格尚未确定,PSA 聚焦于临床/模型参数不确定性") |
| M3 | MAJOR | Methods §2.4, GP floor | GP mortality hazard floor 从 month 60 开始,但未正式分析切换点(60 个月)选择对结果的影响。如果从 month 36 或 month 72 开始,ICER 会怎样变化? | 添加 GP floor start month(如 36/60/84)的敏感性分析,或引用文献说明 60 个月是标准选择 |
| M4 | MINOR | 摘要 vs §3.1 | 摘要写 combo QALY 14.16,正文 §3.1 Table 2 也写 14.16——数值一致 ✓。但折扣 life-years 摘要写 19.01 vs 10.75,正文 §3.1 写 19.01 vs 10.75 ——一致 ✓。Table 2 中 Combo DM cost $524,121 − Pembro $453,082 = $71,039 ✓,与正文一致 | 无需修改(验证通过) |
| M5 | MINOR | Supplementary §S13, PSA covariance | 正文说 PSA 协方差从 NLS Jacobian 推导,但未说明交叉治疗组协方差为何独立抽样(arm-specific parameters sampled independently)。这是常见做法但应明确说明局限 | 在 Methods 中添加一句:"Cross-arm covariance could not be estimated from aggregate KM data; arm-specific survival parameters were therefore sampled independently" — 已有但可更显著 |
| M6 | MINOR | Table S4, 4096 combinations | 所有 top 10 组合的 pembro OS40 = 0.0%,但 combo OS40 = 71–76%。这表明 unconstrained 下 pembro arm 已被 GP floor 完全约束,而 combo arm 尚未——为什么 combo 的 unconstrained OS 在 40 年仍高达 76%? | 在补充材料中添加解释:combo arm 的 log-normal 参数 mu=8.396 远高于 pembro 的 mu=4.841,导致 40 年 extrapolation 仍在 GP floor 之上,因此 hazard floor 从 month 60 开始后会显著约束 combo 的长期生存 |

---

## Reviewer 2 Report — Domain (Clinical)

**总体评价：**  
作为 melanoma oncologist,我认为该文对 KEYNOTE-942 5 年数据的使用总体准确,HR 和生存率引用与 Khattak et al. 2026 报告一致。外部验证(pembro arm vs KEYNOTE-054/CheckMate 238)的框架合理,对人群差异的解释令人信服。但 30–40 年 OS 外推( combo arm 40y OS = 58.7%)在 stage IIIB–IV melanoma 队列中临床上难以置信——这是在 KEYNOTE-942 中位随访仅 5 年、OS 事件仅 7/104 的基础上外推 35 年。Phase 3 INTerpath-001 已宣布阳性结果但 HRs 未公开,论文在标题中用 "Updated" 但正文仍依赖 Phase 2b 数据,应更审慎措辞。

**评分：**
- 临床数据准确性: **8/10** — 关键 HR 和 5 年率引用准确
- 外部效度: **6/10** — 40 年 OS 外推临床上缺乏支撑
- 文献覆盖: **7/10** — 缺少对 KEYNOTE-054 和 CheckMate 238 最新更新的完整引用

**Issues:**

| # | 严重度 | 位置 | 问题 | 建议 |
|---|--------|------|------|------|
| D1 | CRITICAL | §3.3, Table S5 | Combo arm 40 年 OS = 58.7%,即 stage IIIB–IV melanoma 患者到约 101 岁仍有 58.7% 存活。即使有 hazard floor 约束,这个数字在临床上极不现实——KEYNOTE-942 的 5 年 OS 92.2% 基于仅 7 个事件(n=104 combo), extrapolation 至 40 年的不确定性被严重低估 | (1) 在 Limitations 中明确说明 40 年 OS 外推基于模型假设而非临床证据;(2) 添加 20 年 horizon 的主要结果作为保守替代(表中已有:ICER $41,785);(3) 讨论中承认 combo 40y OS 58.7% 是模型产物而非临床预测 |
| D2 | MAJOR | 摘要 Background | "Stage IIB–IV melanoma" 但 KEYNOTE-942 纳入的是 completely resected stage IIIB–IV,INTerpath-001 扩展到 IIB–IV。应区分 | 标题/摘要统一为 "Stage IIIB–IV"(与 KEYNOTE-942 一致);如引用 INTerpath-001 范围,明确标注 |
| D3 | MAJOR | §1.4, §4.6 | 论文标题说 "Updated" 但模型输入完全来自 Phase 2b KEYNOTE-942。Phase 3 INTerpath-001 已宣布阳性但 HRs 未公开。读者可能误以为已整合 Phase 3 数据 | 在摘要 Methods 和 §1.4 明确说明:模型输入基于 Phase 2b KEYNOTE-942;Phase 3 数据将在 HRs 公布后更新 |
| D4 | MAJOR | §3.3, External validation | 外部验证仅用 published landmark values(5y/7y 的 RFS/DMFS/OS),未做 formal statistical calibration(如 chi-squared or visual inspection with confidence bands).模型 pembro 5y RFS 40.5% vs CheckMate 238 51%,差异 -10.5pp 仅以 "higher-risk population" 解释,未量化 | 添加 formal goodness-of-fit 统计量或至少在 figure 中显示 95% CI;讨论中量化人群差异对结果的影响 |
| D5 | MINOR | §2.2, Clinical data | 60-month DMFS 未报告(Supplementary Table S2 显示 "---"),这是已知的 KEYNOTE-942 局限。但论文未充分讨论这个缺失对 r_DM 参数拟合和长期 DMFS 外推的影响 | 在 Methods 中说明 DMFS 在 60 个月缺失,因此 r_DM 拟合仅基于 18/24/36/48 个月数据点;这增加了长期 DMFS 外推不确定性 |
| D6 | MINOR | References | CheckMate 238 的引用(ascierto2023)是 2017 年 NEJM 原始发表,未引用后续 5 年/7 年更新数据(ascierto2025)。Table S7 引用了 ascierto2025 但 references.bib 中需确认是否包含 | 确保 ascierto2025(7 年更新)在 bib 文件中并正确引用 |

---

## Reviewer 3 Report — Perspective (Policy/HTA)

**总体评价：**  
该文对 US payer 视角下的 intismeran 定价问题提供了有价值的早期证据。阈值价格分析($527K/$852K)直接回答了"多贵仍具成本效果"的支付方核心问题。但论文将 cost-effectiveness 与 affordability 的区分停留在一句话(Level §4.3),未提供 budget impact 分析,这对实际覆盖决策是关键缺失。$200,000 价格基于 Jefferies 分析师估计,缺乏对不同定价策略(如基于疗效的定价)的讨论。"no-direct-OS-benefit = dominant"场景虽有趣,但 ΔQALY 仅 0.52 且 dominance 依赖于模型假设(复发预防驱动 QALY 增益),对政策读者来说可能过度乐观。

**评分：**
- 政策相关性: **7/10** — 阈值价格分析对 payer 有直接价值
- 实践价值: **6/10** — 缺少 budget impact 分析
- 假设合理性: **7/10** — $200K 价格假设合理但缺乏对替代定价策略的讨论

**Issues:**

| # | 严重度 | 位置 | 问题 | 建议 |
|---|--------|------|------|------|
| P1 | MAJOR | 全文 | 缺少 budget impact analysis (BIA)。论文多次提到"cost-effectiveness ≠ affordability"(§4.3)但未提供任何 budget impact 量化。对于实际 coverage decision 这是关键缺失 | 添加一个补充表格:不同市场规模假设下(如每年 5,000/10,000/20,000 患者)的 5 年 budget impact |
| P2 | MAJOR | §4.3 | $200,000 价格仅引用 Jefferies 分析师估计,未讨论替代定价策略(如基于疗效的 outcomes-based contract, risk-sharing agreement)。INT 个体化生产特性可能支持价值-based pricing | 添加对替代定价机制的讨论,特别是 outcomes-based pricing 在 INT 类产品中的适用性 |
| P3 | MAJOR | §3.5, Table 4 | "No direct OS benefit"场景报告为 "Dominant"(ΔCost = −$22,543, ΔQALY = 0.52)。但 dominance 仅由 0.52 QALY 差异驱动,且取决于复发预防的 QALY 贡献(RF utility 0.83)和 LR/DM 成本假设。若 RF utility 降至 0.80(DSA 低值),dominance 可能消失 | 在讨论该场景时明确说明 dominance 对 utility 和 cost 假设的敏感性;建议添加该场景的 PSA |
| P4 | MINOR | §4.3 | 提到"formal budget-impact analysis is needed"但未说明 BIA 何时/如何进行 | 建议在 Limitations 中明确承诺将在 Phase 3 数据可用后进行 BIA |
| P5 | MINOR | §3.6, Figure 4 | 阈值价格分析假设 ICER 与价格呈线性关系(图中为直线)。实际中,如果价格变化影响患者选择/市场渗透,关系可能非线性 | 注明线性假设的局限;如可能,讨论价格弹性对阈值的影响 |

---

## Devil's Advocate Report

**最强反方论证：**  
本文的核心论证——intismeran + pembrolizumab 在 $200,000 定价下具成本效果——几乎完全依赖于 GP mortality hazard floor 下 log-normal OS 外推的 combo arm 长期生存轨迹。40 年 OS 58.7% 对 stage IIIB–IV melanoma 患者而言,意味着到 101 岁仍有近六成人存活,这在临床上令人震惊。论文巧妙地使用 hazard floor(而非 survival cap)框架:floor 仅约束死亡率不低于 GP,但不约束生存曲线——因此 combo arm 在早期 5 年 RFS 优势(68.8% vs 49.1%)产生的"超额存活"可沿 log-normal 尾部延伸 35 年而不被 GP 曲线拉回。这是一种结构性乐观(structural optimism):模型的数学构造本身就偏向于将短期 RFS 收益转化为巨大的终身 QALY 增益。OS HR 0.471 的 95% CI 上限达 1.345(即可能无 OS 获益),但论文的 base case 和所有场景均假设 combo 有 OS 优势——"no direct OS benefit"场景虽存在,但被标记为"Dominant"(成本更低),暗示无论 OS 是否获益,intismeran 都值得用,这对一个仅基于 7 个 OS 事件的 Phase 2b 试验来说,结论过强。

**Issue List:**

| # | 严重度 | 位置 | 攻击点 | 作者如何可能回应 |
|---|--------|------|--------|----------------|
| DA1 | CRITICAL | §3.3, Table S5 | Combo 40y OS = 58.7% 在 stage IIIB–IV melanoma 中临床上不可能。Hazard floor 框架是否被有意选择以产生比 survival cap 更乐观的 extrapolation? | 作者可能回应:hazard floor 是 NICE DSU TSD 14 推荐的标准方法;survival cap 会在 OS > GP 时错误地强制两曲线重合。但问题是:floor 从 month 60 开始,而 combo arm 在 60 个月时 OS 已达 92.2%,远高于 GP 的 95.5%——这意味着 floor 在早期并不激活,combo 的 log-normal 尾部可以在少约束的情况下延伸 |
| DA2 | CRITICAL | 摘要 | "92.2% vs 71.3% 5-year OS" 仅基于 7 个事件/arm(n=104),95% CI 极宽(0.165–1.345)。将此点估计放在摘要中作为"Results"呈现,可能误导读者认为这是一个精确估计 | 作者可能回应:这是已发表的 KEYNOTE-942 5 年数据的忠实引用。但问题在于:摘要中将其作为确定性结果呈现,而非作为高度不确定的探索性终点 |
| DA3 | MAJOR | §3.5, Table 4 | "No direct OS benefit"场景 ΔQALY = 0.52 的 dominance 极其脆弱——仅依赖 RF utility(0.83)和 LR/DM 成本假设。若 RF utility 降至 0.80,dominance 消失 | 作者可能回应:已在 DSA 中测试了 RF utility 的敏感性。但 DSA 仅在 base case 下测试,未在"no direct OS"场景下重复 |
| DA4 | MAJOR | §1.4, §4.6 | Phase 3 INTerpath-001 已宣布阳性但论文仍完全基于 Phase 2b 数据。在 Phase 3 HRs 未公开的情况下,用"Updated"标题是否在利用 Phase 3 的光环效应( halo effect)来增强 Phase 2b 分析的可信度? | 作者可能回应:"Updated"指更新了 5 年随访数据(从 ISPOR 2024 的早期数据)。但标题确实可能被解读为"已整合 Phase 3 数据" |
| DA5 | MINOR | Table S3 | Combo OS 选 log-normal(AIC -50.08),pembro OS 选 Weibull(AIC -32.00)。但论文对 combo arm 统一使用 log-normal(所有 6 条曲线),pembro arm 也用 log-normal(parsimony)。这意味着 pembro OS 不是其最优拟合(Weibull 更好)——这是否为了保持模型一致性而牺牲了 pembro arm 的拟合质量? | 作者可能回应:parsimony 和跨组件一致性是标准做法。但应量化 Weibull pembro OS 对 base case ICER 的影响(已在 Weibull 场景中测试:ICER $59,959) |

**被忽略的替代解释：**
论文未讨论的一种可能性:combo arm 的 RFS 优势可能部分来自 lead-time bias(更频繁的监测导致更早检测复发),而非真实的治疗效果。KEYNOTE-942 是 open-label 试验,复发检测可能受组间监测差异影响。

**"So What?" Test 结论：**  
论文的核心价值(阈值价格分析)不依赖于 OS extrapolation 的精确性——即使 OS 外推有偏,价格-ICER 线性关系仍然成立。但 base case ICER $49,841 的精确性确实依赖于 extrapolation 质量。建议将 20 年 horizon 结果($41,785)作为更保守的主报告,40 年结果作为补充。

---

## Editorial Decision (Phase 2 Synthesis)

### Consensus Issues (Multiple Reviewers Agree)

| Issue | Reviewers | Severity | Action Required |
|-------|-----------|----------|-----------------|
| 标题/摘要分期不一致(IIIB–IV vs IIB–IV) | EIC, R2 | MAJOR | 统一为 Stage IIIB–IV |
| Combo 40y OS = 58.7% 临床不合理 | R2, DA | CRITICAL | 添加明确 disclaimers;考虑将 20-year horizon 作为 co-primary |
| "Robust across scenarios" 措辞过度 | EIC | MAJOR | 修改为 "largely robust, with ICER exceeding $100K only under hazard-convergence waning" |
| 缺少 budget impact analysis | R3 | MAJOR | 添加 BIA 表格或在 Limitations 中明确承诺 |
| Intismeran 价格未纳入 PSA | R1 | MAJOR | 纳入 PSA 或明确说明为何固定 |
| Phase 2b vs Phase 3 数据混淆 | R2, DA | MAJOR | 在摘要和 §1.4 明确说明模型输入基于 Phase 2b |
| 结果/讨论重复 | EIC | MAJOR | 精简 §4.2 重复内容 |
| Table S8 DSA 参数范围疑似错误 | R1 | CRITICAL | 核实 OS mu combo DSA 范围 |

### Devil's Advocate CRITICAL Adjudication

| DA Issue | Verdict | Rationale |
|----------|---------|-----------|
| DA1: Combo 40y OS 58.7% | **Validated** — 确实是模型产物而非临床预测。建议在 Limitations 和 Results 中明确标注 | Hazard floor 是标准方法但 40 年 extrapolation 仍是假设驱动 |
| DA2: 5y OS 92.2% with 7 events | **Partially validated** — 点估计引用准确但应在摘要中标注不确定性 | 建议在摘要 Results 中加注 "(HR 0.471; 95% CI 0.165–1.345; 7 events)" |

### Editorial Decision Letter

**Decision: Minor Revision**

Dear Dr. Jin,

Thank you for submitting your manuscript "Cost-Effectiveness of Intismeran Autogene Plus Pembrolizumab Versus Pembrolizumab Alone as Adjuvant Therapy for Resected Stage IIIB–IV Melanoma" to PharmacoEconomics.

Your manuscript has been evaluated by the Editor-in-Chief and three independent peer reviewers with expertise in health economic modeling, melanoma oncology, and HTA policy. The reviewers agreed that your paper addresses an important and timely topic — the cost-effectiveness of the first mRNA-based individualized neoantigen therapy — and commended the transparency of your modeling approach, including the open-source code, systematic distribution selection, and comprehensive sensitivity analyses.

However, several issues require attention before the manuscript can be accepted:

**Required revisions:**

1. **Title/Abstract consistency (EIC, R2):** Unify the population staging to "Stage IIIB–IV" throughout (currently "Stage IIB–IV" appears in the abstract background).

2. **Long-term OS extrapolation disclaimers (R2, DA — CRITICAL):** The 40-year OS of 58.7% for the combination arm in a stage IIIB–IV melanoma cohort is a model output, not a clinical prediction. Add explicit disclaimers in the Results, Discussion, and Limitations sections. Consider presenting the 20-year horizon result ($41,785/QALY) as a co-primary alongside the 40-year base case.

3. **"Robust"措辞 (EIC):** Replace "Results were robust across scenarios" with language that acknowledges the hazard-convergence waning scenario exceeds $100,000/QALY.

4. **PSA price uncertainty (R1):** Either include the intismeran acquisition price as a probabilistic parameter in the PSA, or explicitly justify why it was held fixed despite being the most sensitive parameter in the DSA.

5. **Table S8 DSA parameter ranges (R1 — CRITICAL):** Verify and correct the OS μ combo DSA range in Table S8 (currently showing 4.228–4.628 against a base case of 8.268).

6. **Phase 2b vs Phase 3 clarity (R2, DA):** Explicitly state in the abstract Methods and Section 1.4 that model inputs are based on Phase 2b KEYNOTE-942 data, with Phase 3 data to be incorporated when hazard ratios become available.

7. **Budget impact (R3):** Add a budget impact analysis table in the supplement or explicitly commit to one in the Limitations section.

8. **Redundancy (EIC):** Reduce repetition between Results §3.2 and Discussion §4.2 regarding cost/QALY decomposition.

The revised manuscript should be returned within 30 days.

Sincerely,  
Editor-in-Chief, PharmacoEconomics

---

### Revision Roadmap (Prioritized)

| Priority | Issue | Reviewer | Effort | Description |
|----------|-------|----------|--------|-------------|
| P0 | DA1: 40y OS disclaimers | R2, DA | Low | Add text disclaimers + consider 20-year co-primary |
| P0 | M1: Table S8 DSA range | R1 | Low | Verify/correct OS μ combo DSA range |
| P1 | E1: Title/abstract staging | EIC, R2 | Trivial | Unify to Stage IIIB–IV |
| P1 | E2: "Robust"措辞 | EIC | Trivial | Modify to "largely robust" |
| P1 | M2: PSA price uncertainty | R1 | Medium | Add price to PSA or justify fixed |
| P1 | D3: Phase 2b vs 3 clarity | R2, DA | Low | Add clarifying text |
| P2 | P1: Budget impact | R3 | Medium | Add BIA table |
| P2 | E3: Discussion redundancy | EIC | Low | Deduplicate §4.2 |
| P3 | D4: External validation stats | R2 | Medium | Add formal calibration |
| P3 | M3: GP floor start sensitivity | R1 | Medium | Add sensitivity analysis |

---

## Addendum: EIC 完整报告补充发现（GLM-5.3 EIC，22:40）

| # | 严重度 | 位置 | 问题 | 建议 |
|---|--------|------|------|------|
| A1 | 中度 | Table S8 (supplementary.tex L354-356) | OS μ 数值口径矛盾：S8 用旧 V1 值 combo=8.268 / pembro=3.725，但正文 DSA 与 S12 用 V2 canonical 值 8.396/4.841。同稿内两套参数并存 | 从 CEAModelV2 重生成 Table S8 全部 fitted 参数行 |
| A2 | 中度 | 摘要 Methods | 未点明疗效输入证据等级（Phase 2b, n=157），读者仅凭摘要会高估结论强度 | 摘要 Methods 加 "based on Phase 2b KEYNOTE-942 5-year efficacy data" |
| A3 | 轻度 | manuscript.tex L371 | Conclusion 段末与 \section*{Declarations} 粘连（缺空行），编译后段落不换行 | Declarations 前补空行 |
| A4 | 建议 | 全文 | "Phase 3 读出后首个 CEA" 这一核心新颖性卖点未显式声明 | 在 Intro/Discussion 明确一句 |
