"""
Comprehensive re-run of cea-intismeran with:
1. GP mortality constraint: hazard floor from 60 months (primary)
2. Unconstrained log-normal (methodological comparison)
3. Treatment-effect waning (conservative)
4. No additional OS hazard benefit beyond 5 years (strong conservative)

Structural integrity: RFS ≤ DMFS ≤ OS at every month, state sum = 1.0.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import norm
from model import ModelParams, CEAModel, lognorm_surv, weibull_surv, general_pop_surv
from model import discount_factors  # noqa: E402

# ── Life table monthly hazards (proper interpolation) ──
# Repo copy first (data/life_table_surv.json); /tmp kept as a legacy fallback.
_LT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "data", "life_table_surv.json")
if not os.path.exists(_LT_PATH):
    _LT_PATH = "/tmp/life_table_surv.json"

with open(_LT_PATH) as f:
    LT = json.load(f)
LT_surv_year = np.array([LT["survival_by_year"].get(str(y), 0.0) for y in range(0, 43)])
# Monthly GP survival: linearly interpolate survival (or use hazard interpolation)
# Build monthly survival via annual hazard compounding
annual_hazard = np.zeros(42)
for y in range(1, 42):
    annual_hazard[y] = -np.log(max(LT_surv_year[y] / max(LT_surv_year[y-1], 1e-12), 1e-12))
# Monthly survival: interpolate the annual log-survival linearly
monthly_log_surv = np.interp(np.arange(480) / 12.0, np.arange(0, 43), np.log(np.maximum(LT_surv_year, 1e-12)))
LT_surv_monthly = np.exp(monthly_log_surv)
# Monthly hazard: -ln(S(m+1)/S(m))
LT_hazard_monthly = np.zeros(480)
for mo in range(479):
    LT_hazard_monthly[mo] = -np.log(max(LT_surv_monthly[mo+1] / max(LT_surv_monthly[mo], 1e-12), 1e-12))
LT_hazard_monthly[-1] = LT_hazard_monthly[-2]

# ── Modified model with hazard-based GP constraint ──
class CEAModelV2(CEAModel):
    """V2: hazard-based GP mortality floor from 60 months, both arms."""

    def _os_dist(self, arm: str) -> str:
        """Effective OS distribution for an arm: per-arm override, else the shared field."""
        prefix = "combo" if arm == "combo" else "pembro"
        return getattr(self.p, f"os_distribution_{prefix}", "") or self.p.os_distribution

    def _os_base(self, mo: np.ndarray, arm: str) -> np.ndarray:
        """Fitted OS curve for one arm, before the GP floor / waning / no-direct-OS steps."""
        prefix = "combo" if arm == "combo" else "pembro"
        if self._os_dist(arm) == "weibull":
            scale = getattr(self.p, f"os_weibull_scale_{prefix}")
            shape = getattr(self.p, f"os_weibull_shape_{prefix}")
            return weibull_surv(mo, scale, shape)
        mu = getattr(self.p, f"os_mu_{prefix}")
        sigma = getattr(self.p, f"os_sigma_{prefix}")
        return lognorm_surv(mo, mu, sigma)

    def _pembro_os(self, mo: np.ndarray) -> np.ndarray:
        """Pembro OS with the GP floor applied (the comparison curve for waning / no-direct-OS)."""
        p_os = self._os_base(mo, "pembro")
        if self.p.constraint_general_pop:
            p_os = self._apply_gp_hazard_floor(p_os, mo)
        return p_os

    def _os_fn(self, mo: np.ndarray, arm: str) -> np.ndarray:
        prefix = "combo" if arm == "combo" else "pembro"
        # Base OS -- honours per-arm distribution overrides (os_distribution_combo /
        # os_distribution_pembro); reading only p.os_distribution silently ignored them.
        os_val = self._os_base(mo, arm)

        # GP mortality hazard floor from 60 months onward
        if self.p.constraint_general_pop:
            os_val = self._apply_gp_hazard_floor(os_val, mo)

        # Treatment-effect waning
        if self.p.treatment_waning and arm == "combo":
            p_os = self._pembro_os(mo)
            # Wane from full benefit at 5y to zero at 20y
            t_y = self.t[:len(mo)]
            wane = np.clip((t_y - 5) / (20 - 5), 0, 1)
            benefit = os_val - p_os
            os_val = p_os + benefit * (1 - wane)

        # No additional OS hazard benefit beyond the observed window: from month 60
        # onward the combination arm takes the pembrolizumab mortality hazard, but the
        # curve is integrated forward from the combination arm's own survival at month
        # 60. Setting the curves equal outright (as before) would have dropped survival
        # from ~92% to ~71% within one month -- an implausible discontinuity.
        if getattr(self.p, 'no_os_hazard_benefit', False) and arm == "combo":
            mask_post = mo > 60
            idx = np.nonzero(mask_post)[0]
            if len(idx):
                p_os = self._pembro_os(mo)
                h_p = -np.log(np.clip(p_os[1:] / np.clip(p_os[:-1], 1e-12, None),
                                      1e-12, None))
                s = os_val[idx[0] - 1]
                for i in idx:
                    s *= float(np.exp(-h_p[i - 1]))
                    os_val[i] = s

        return os_val

    def _apply_gp_hazard_floor(self, os_val: np.ndarray, mo: np.ndarray) -> np.ndarray:
        """Apply GP mortality hazard floor from 60 months onward."""
        n = len(os_val)
        # Convert to hazard
        h_model = np.zeros(n)
        for i in range(n - 1):
            h_model[i] = -np.log(max(os_val[i + 1] / max(os_val[i], 1e-12), 1e-12))
        # Apply floor from gp_floor_start_months (default 60)
        mask = mo >= getattr(self.p, 'gp_floor_start_months', 60)
        gp_h = np.interp(np.arange(n), np.arange(480), LT_hazard_monthly)
        h_final = np.where(mask, np.maximum(h_model, gp_h), h_model)
        h_final[-1] = h_final[-2]  # last point
        # Reconstruct survival S(i) = exp(-sum_{j<i} h_j): index i is month i.
        cum_h = np.concatenate([[0.0], np.cumsum(h_final[:-1])])
        new_os = np.exp(-cum_h)
        # Ensure monotonic
        for i in range(1, n):
            new_os[i] = min(new_os[i], new_os[i-1])
        return new_os


def run_diagnostics(scenario_name: str, p: ModelParams, extra_attrs: dict = None):
    """Run model, print full diagnostics, return results."""
    if extra_attrs:
        for k, v in extra_attrs.items():
            setattr(p, k, v)
    m = CEAModelV2(p)
    t = m.t
    mo = t * 12
    sur = m._survival()
    states = m._state_proportions(sur)
    res = m.run()

    print(f"\n{'='*80}")
    print(f"SCENARIO: {scenario_name}")
    print(f"GP constraint={p.constraint_general_pop}, "
          f"Waning={p.treatment_waning}, "
          f"No additional OS hazard benefit={getattr(p, 'no_os_hazard_benefit', False)}")
    print(f"{'='*80}")

    # ── Verifications ──
    print("\n--- Structural Integrity (RFS ≤ DMFS ≤ OS, state sum = 1) ---")
    violations = {"rfs_dmfs": 0, "dmfs_os": 0, "sum": 0.0}
    for i in range(len(t)):
        for arm in ["combo", "pembro"]:
            rfs = sur[f"rfs_{arm}"][i]
            dmfs = sur[f"dmfs_{arm}"][i]
            osv = sur[f"os_{arm}"][i]
            if rfs > dmfs + 1e-10:
                violations["rfs_dmfs"] += 1
            if dmfs > osv + 1e-10:
                violations["dmfs_os"] += 1
            st = states[arm]
            ssum = st["rf"][i] + st["lr"][i] + st["dm"][i] + st["dead"][i]
            violations["sum"] = max(violations["sum"], abs(ssum - 1.0))
    print(f"  RFS ≤ DMFS violations: {violations['rfs_dmfs']} / {len(t)} months")
    print(f"  DMFS ≤ OS violations:  {violations['dmfs_os']} / {len(t)} months")
    print(f"  Max state sum deviation: {violations['sum']:.2e}")
    ok = violations["rfs_dmfs"] == 0 and violations["dmfs_os"] == 0 and violations["sum"] < 1e-10
    print(f"  OVERALL: {'✓ PASS' if ok else '✗ FAIL'}")

    # ── OS at key years ──
    print(f"\n--- OS at Key Time Points ---")
    print(f"{'Year':>5s} {'Combo OS':>10s} {'Pembro OS':>10s} {'Diff':>10s}")
    for yr in [5, 10, 20, 30, 40]:
        idx = int(yr / p.cycle_length)
        if idx < len(t):
            cos = sur[f"os_combo"][idx]
            pos = sur[f"os_pembro"][idx]
            print(f"{yr:5d} {cos:>10.4f} {pos:>10.4f} {cos-pos:>10.4f}")

    # ── Discounted LY ──
    disc = discount_factors(p.discount_rate, t)
    for arm, label in [("combo", "Combo"), ("pembro", "Pembro")]:
        ly_disc = np.sum(sur[f"os_{arm}"] * disc) * p.cycle_length
        res[arm]['ly_disc'] = ly_disc

    # ── Results summary ──
    print(f"\n--- Results ---")
    print(f"{'':>12s} {'Combo':>12s} {'Pembro':>12s} {'Δ':>12s}")
    c, p_n = res["combo"], res["pembro"]
    dly = c["ly_disc"] - p_n["ly_disc"]
    dq = c["qaly"] - p_n["qaly"]
    dc = c["cost"] - p_n["cost"]
    icer = dc / dq if dq > 0 else float('inf')
    print(f"Disc LY   {c['ly_disc']:>12.2f} {p_n['ly_disc']:>12.2f} {dly:>12.2f}")
    print(f"Disc QALY {c['qaly']:>12.2f} {p_n['qaly']:>12.2f} {dq:>12.2f}")
    print(f"Cost      ${c['cost']:>9,.0f} ${p_n['cost']:>9,.0f} ${dc:>9,.0f}")
    print(f"ICER      {'':>12s} {'':>12s} ${icer:>8,.0f}/QALY")

    # ── State occupancy at key years ──
    print(f"\n--- State Occupancy (Combo) ---")
    print(f"{'Year':>5s} {'RF':>8s} {'LR':>8s} {'DM':>8s} {'Dead':>8s} {'OS':>8s}")
    for yr in [5, 10, 20, 30, 40]:
        idx = int(yr / p.cycle_length)
        if idx < len(t):
            st = states["combo"]
            print(f"{yr:5d} {st['rf'][idx]:>8.4f} {st['lr'][idx]:>8.4f} {st['dm'][idx]:>8.4f} {st['dead'][idx]:>8.4f} {sur['os_combo'][idx]:>8.4f}")

    print(f"\n--- State Occupancy (Pembro) ---")
    print(f"{'Year':>5s} {'RF':>8s} {'LR':>8s} {'DM':>8s} {'Dead':>8s} {'OS':>8s}")
    for yr in [5, 10, 20, 30, 40]:
        idx = int(yr / p.cycle_length)
        if idx < len(t):
            st = states["pembro"]
            print(f"{yr:5d} {st['rf'][idx]:>8.4f} {st['lr'][idx]:>8.4f} {st['dm'][idx]:>8.4f} {st['dead'][idx]:>8.4f} {sur['os_pembro'][idx]:>8.4f}")

    return {
        "name": scenario_name,
        "combo_os": {yr: float(sur[f"os_combo"][min(int(yr/p.cycle_length), len(sur[f"os_combo"])-1)]) for yr in [5, 10, 20, 30, 40]},
        "pembro_os": {yr: float(sur[f"os_pembro"][min(int(yr/p.cycle_length), len(sur[f"os_pembro"])-1)]) for yr in [5, 10, 20, 30, 40]},
        "combo_ly_disc": float(c["ly_disc"]),
        "pembro_ly_disc": float(p_n["ly_disc"]),
        "combo_qaly": float(c["qaly"]),
        "pembro_qaly": float(p_n["qaly"]),
        "combo_cost": float(c["cost"]),
        "pembro_cost": float(p_n["cost"]),
        "dly": float(dly),
        "dqaly": float(dq),
        "dcost": float(dc),
        "icer": float(icer),
        "structural_ok": ok,
        "violations": violations,
    }


# ══════════════════════════════════════════
# RUN ALL 4 SCENARIOS
# ══════════════════════════════════════════

results = []

# Scenario 1: GP Mortality Constrained (PRIMARY)
p = ModelParams(constraint_general_pop=True)
results.append(run_diagnostics("GP Mortality Constrained (PRIMARY)", p))

# Scenario 2: Unconstrained log-normal (methodological comparison)
p = ModelParams(constraint_general_pop=False)
results.append(run_diagnostics("Unconstrained Log-Normal (methodological comparison)", p))

# Scenario 3: Treatment-effect waning
p = ModelParams(constraint_general_pop=True, treatment_waning=True)
results.append(run_diagnostics("Treatment-Effect Waning (conservative)", p))

# Scenario 4: No additional OS hazard benefit beyond 5 years
p = ModelParams(constraint_general_pop=True)
results.append(run_diagnostics("No additional OS hazard benefit beyond 5 years (strong conservative)", p, {"no_os_hazard_benefit": True}))

# Scenario 5: RFS Calibration (anchor pembrolizumab RFS to 49.1% at 60 months)
# Scale factor = 0.913 / 0.753 ≈ 1.212 to achieve RFS(60) = 49.1%
# where 0.913 = 49.1% / 53.8% (target r_LR) and 0.753 = 40.5% / 53.8% (base r_LR)
p = ModelParams(constraint_general_pop=True, rLR_pembro_scale=1.212)
results.append(run_diagnostics("RFS Calibration (pembro anchored to 49.1%)", p))

# Scenario 6: LR cost +50%
p = ModelParams(constraint_general_pop=True, cost_lr_monthly=4_500)
results.append(run_diagnostics("LR cost +50%", p))

# Scenario 7: LR cost -50%
p = ModelParams(constraint_general_pop=True, cost_lr_monthly=1_500)
results.append(run_diagnostics("LR cost -50%", p))

# Scenario 8: LR utility +0.1
p = ModelParams(constraint_general_pop=True, util_lr=0.74)
results.append(run_diagnostics("LR utility +0.1", p))

# Scenario 9: LR utility -0.1
p = ModelParams(constraint_general_pop=True, util_lr=0.54)
results.append(run_diagnostics("LR utility -0.1", p))

# ── Summary table ──
print(f"\n\n{'='*80}")
print("FINAL SUMMARY: ALL 5 SCENARIOS")
print(f"{'='*80}")
print(f"{'Scenario':<40s} {'Combo LY':>8s} {'Pembro LY':>8s} {'ΔLY':>8s} {'Combo Q':>8s} {'Pembro Q':>8s} {'ΔQALY':>8s} {'ΔCost':>10s} {'ICER':>10s}")
print("-" * 100)
for r in results:
    print(f"{r['name']:<40s} {r['combo_ly_disc']:>8.2f} {r['pembro_ly_disc']:>8.2f} {r['dly']:>8.2f} "
          f"{r['combo_qaly']:>8.2f} {r['pembro_qaly']:>8.2f} {r['dqaly']:>8.2f} "
          f"${r['dcost']:>8,.0f} ${r['icer']:>8,.0f}")

# ── OS comparison table ──
print(f"\n{'='*80}")
print("OS AT KEY TIME POINTS BY SCENARIO")
print(f"{'='*80}")
for yr in [5, 10, 20, 30, 40]:
    print(f"\nYear {yr}:")
    print(f"{'Scenario':<40s} {'Combo OS':>10s} {'Pembro OS':>10s} {'Diff':>10s}")
    for r in results:
        cos = r['combo_os'][yr]
        pos = r['pembro_os'][yr]
        print(f"{r['name']:<40s} {cos:>10.4f} {pos:>10.4f} {cos-pos:>10.4f}")

print(f"\n{'='*80}")
print("COMPLETE — Manuscript not modified. Output above for review.")
print(f"{'='*80}")
