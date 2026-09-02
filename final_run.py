"""
Re-run primary with corrected life table (age 61, 64% male).
Then run PSA, DSA, scenarios, price threshold.
No manuscript update until complete.
"""
import sys, os, json
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from scipy.stats import norm, gamma, beta
from model import ModelParams, lognorm_surv, weibull_surv

# ══════════════════════════════════════════
# 1. Life table (updated: age 61, 64% male)
# ══════════════════════════════════════════
with open("/tmp/life_table_surv.json") as f:
    LT = json.load(f)
LT_surv_year = np.array([LT["survival_by_year"].get(str(y), 0.0) for y in range(0, 43)])
monthly_log_surv = np.interp(np.arange(480) / 12.0, np.arange(0, 43), np.log(np.maximum(LT_surv_year, 1e-12)))
LT_surv_monthly = np.exp(monthly_log_surv)
LT_hazard_monthly = np.zeros(480)
for mo in range(479):
    LT_hazard_monthly[mo] = -np.log(max(LT_surv_monthly[mo+1] / max(LT_surv_monthly[mo], 1e-12), 1e-12))
LT_hazard_monthly[-1] = LT_hazard_monthly[-2]

# ══════════════════════════════════════════
# 2. Model (V2 with hazard-based GP constraint)
# ══════════════════════════════════════════
from rerun_primary import CEAModelV2, run_diagnostics

def run_base(p: ModelParams) -> dict:
    """Run base case, return full results dict."""
    m = CEAModelV2(p)
    t = m.t
    sur = m._survival()
    states = m._state_proportions(sur)
    res = m.run()
    disc = np.exp(-p.discount_rate * t)
    for arm in ["combo", "pembro"]:
        ly_disc = np.sum(sur[f"os_{arm}"] * disc) * p.cycle_length
        res[arm]['ly_disc'] = ly_disc
    return {"m": m, "t": t, "sur": sur, "states": states, "res": res, "p": p}

# ── Primary: GP-constrained ──
p = ModelParams(constraint_general_pop=True)
base = run_base(p)
r = base["res"]
icer = base["m"].icer(r)
print(f"{'='*60}")
print("PRIMARY BASE CASE (GP mortality constrained, age 61 pooled)")
print(f"{'='*60}")
print(f"Combo:  LY_disc={r['combo']['ly_disc']:.2f}  QALY={r['combo']['qaly']:.2f}  Cost=${r['combo']['cost']:,.0f}")
print(f"Pembro: LY_disc={r['pembro']['ly_disc']:.2f}  QALY={r['pembro']['qaly']:.2f}  Cost=${r['pembro']['cost']:,.0f}")
print(f"ΔLY={r['combo']['ly_disc']-r['pembro']['ly_disc']:.2f}  ΔQALY={r['combo']['qaly']-r['pembro']['qaly']:.2f}  ΔCost=${r['combo']['cost']-r['pembro']['cost']:,.0f}")
print(f"ICER = ${icer:,.0f}/QALY")

# ══════════════════════════════════════════
# 3. PSA (1,000 iterations)
# ══════════════════════════════════════════
print(f"\n{'='*60}")
print("PSA (1,000 iterations)")
print(f"{'='*60}")

rng = np.random.default_rng(42)
icers_psa = []
costs_diff = []
qalys_diff = []
n_ok = 0

for _ in range(1000):
    pp = ModelParams(constraint_general_pop=True)
    # Perturb costs (gamma, CV=0.2)
    for attr, val in [("cost_keytruda_annual", 220_896), ("cost_intismeran", 200_000),
                      ("cost_lr_monthly", 3_000), ("cost_dm_monthly", 12_000)]:
        cv = 0.2
        setattr(pp, attr, rng.gamma(1/cv**2, val*cv**2))
    # Perturb utilities (beta, SE=0.03)
    for attr, val in [("util_rf", 0.83), ("util_lr", 0.64), ("util_dm", 0.55)]:
        se = 0.03
        alpha = val * (val*(1-val)/se**2 - 1)
        beta_p = (1-val) * (val*(1-val)/se**2 - 1)
        setattr(pp, attr, rng.beta(max(alpha,1), max(beta_p,1)))
    # Perturb OS params (mu ~ N(0, 0.1))
    pp.os_mu_combo += rng.normal(0, 0.1)
    pp.os_mu_pembro += rng.normal(0, 0.1)

    m = CEAModelV2(pp)
    try:
        r_psa = m.run()
        ic = m.icer(r_psa)
        if 0 < ic < 1e6:
            icers_psa.append(ic)
            costs_diff.append(r_psa["combo"]["cost"] - r_psa["pembro"]["cost"])
            qalys_diff.append(r_psa["combo"]["qaly"] - r_psa["pembro"]["qaly"])
            n_ok += 1
    except:
        continue

icers_arr = np.array(icers_psa)
print(f"Valid iterations: {n_ok}")
print(f"Mean ICER: ${np.mean(icers_arr):,.0f}")
print(f"Median ICER: ${np.median(icers_arr):,.0f}")
print(f"95% CI: ${np.percentile(icers_arr, 2.5):,.0f} – ${np.percentile(icers_arr, 97.5):,.0f}")
for thresh, label in [(50_000, "$50K"), (100_000, "$100K"), (150_000, "$150K")]:
    pce = np.mean(icers_arr <= thresh)
    print(f"P(CE at {label}): {pce:.4f}")

# ══════════════════════════════════════════
# 4. DSA (one-way)
# ══════════════════════════════════════════
print(f"\n{'='*60}")
print("DSA (One-Way Sensitivity Analysis)")
print(f"{'='*60}")

dsa_params = {
    "Utility RF": ("util_rf", 0.83, 0.80, 0.86),
    "Utility LR": ("util_lr", 0.64, 0.61, 0.67),
    "Utility DM": ("util_dm", 0.55, 0.52, 0.58),
    "Keytruda cost": ("cost_keytruda_annual", 220_896, 176_717, 265_075),
    "Intismeran cost": ("cost_intismeran", 200_000, 160_000, 240_000),
    "LR monthly cost": ("cost_lr_monthly", 3_000, 2_400, 3_600),
    "DM monthly cost": ("cost_dm_monthly", 12_000, 9_600, 14_400),
    "Discount rate": ("discount_rate", 0.03, 0.0, 0.05),
    "OS mu combo": ("os_mu_combo", p.os_mu_combo, p.os_mu_combo-0.2, p.os_mu_combo+0.2),
    "OS mu pembro": ("os_mu_pembro", p.os_mu_pembro, p.os_mu_pembro-0.2, p.os_mu_pembro+0.2),
}

for label, (attr, base_val, lo, hi) in sorted(dsa_params.items()):
    icers_dsa = []
    for val, tag in [(lo, "low"), (hi, "high")]:
        pp = ModelParams(constraint_general_pop=True)
        setattr(pp, attr, val)
        m = CEAModelV2(pp)
        r_dsa = m.run()
        icers_dsa.append(m.icer(r_dsa))
    print(f"  {label:<25s}  low=${icers_dsa[0]:>,.0f}  high=${icers_dsa[1]:>,.0f}  range=${abs(icers_dsa[1]-icers_dsa[0]):>,.0f}")

# ══════════════════════════════════════════
# 5. Scenarios
# ══════════════════════════════════════════
print(f"\n{'='*60}")
print("SCENARIO ANALYSES")
print(f"{'='*60}")

scenarios = [
    ("GP Mortality Constrained (PRIMARY)", {"constraint_general_pop": True}),
    ("Unconstrained log-normal", {"constraint_general_pop": False}),
    ("Treatment-effect waning", {"constraint_general_pop": True, "treatment_waning": True}),
    ("No direct OS benefit", {"constraint_general_pop": True, "no_direct_os_benefit": True}),
    ("0% discount rate", {"constraint_general_pop": True, "discount_rate": 0.0}),
    ("5% discount rate", {"constraint_general_pop": True, "discount_rate": 0.05}),
    ("Weibull OS", {"constraint_general_pop": True, "os_distribution": "weibull"}),
    ("10-year horizon", {"constraint_general_pop": True, "time_horizon_years": 10}),
    ("20-year horizon", {"constraint_general_pop": True, "time_horizon_years": 20}),
    ("Waning + GP + 20y (most conservative)", {"constraint_general_pop": True, "treatment_waning": True, "time_horizon_years": 20}),
]

print(f"{'Scenario':<45s} {'ΔQALY':>8s} {'ΔCost':>10s} {'ICER':>10s}")
for name, kwargs in scenarios:
    extra = {k: v for k, v in kwargs.items() if k == "no_direct_os_benefit"}
    pp = ModelParams(**{k: v for k, v in kwargs.items() if k != "no_direct_os_benefit"})
    for k, v in extra.items():
        setattr(pp, k, v)
    m = CEAModelV2(pp)
    r_sc = m.run()
    dq = r_sc["combo"]["qaly"] - r_sc["pembro"]["qaly"]
    dc = r_sc["combo"]["cost"] - r_sc["pembro"]["cost"]
    ic = dc / dq if dq > 0 else float('inf')
    print(f"{name:<45s} {dq:>8.2f} ${dc:>8,.0f} ${ic:>8,.0f}")

# ══════════════════════════════════════════
# 6. Price threshold analysis
# ══════════════════════════════════════════
print(f"\n{'='*60}")
print("PRICE THRESHOLD ANALYSIS")
print(f"{'='*60}")

for thresh in [100_000, 150_000, 200_000]:
    # Binary search for max intismeran price at each WTP
    lo, hi = 50_000, 2_000_000
    for _ in range(30):
        mid = (lo + hi) / 2
        pp = ModelParams(constraint_general_pop=True, cost_intismeran=mid)
        m = CEAModelV2(pp)
        r_pt = m.run()
        ic = m.icer(r_pt)
        if ic <= thresh:
            lo = mid
        else:
            hi = mid
    print(f"  At WTP=${thresh:,}/QALY: max intismeran price = ${lo:,.0f}")

print(f"\n{'='*60}")
print("COMPLETE - all results ready for manuscript update")
print(f"{'='*60}")