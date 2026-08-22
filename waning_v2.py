"""
Waning v2: hazard-convergence implementation.
- 0-60 months: identical to base case (KEYNOTE-942 observed evidence preserved)
- 60+ months: combo's disease-related mortality advantage wanes (5-10y linear),
  both arms converge to same MORTALITY HAZARD (not survival curves).
  r_dm (DMFS/OS) and r_lr (RFS/DMFS) also wane to pembro levels (HR->1 by 10y).
- OS survival gap at 60 months is preserved.
- GP mortality floor continues to apply (auto-satisfied via pembro hazard).
"""
import sys, os, json
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from model import ModelParams, CEAModel, lognorm_surv, weibull_surv
from rerun_primary import CEAModelV2, LT_hazard_monthly


class CEAModelV3(CEAModelV2):
    """V3: Waning v2 — hazard convergence from 60 months (primary OS path preserved)."""

    def __init__(self, p: ModelParams):
        super().__init__(p)
        # Cache base (non-waned) combo/pembro curves for waning reference
        self._base_sur = {}

    def _waning_factor(self, mo):
        """Linear 0->1 from wane_start to wane_end months, 0 before."""
        start = getattr(self.p, 'wane_start_month', 60)
        end = getattr(self.p, 'wane_end_month', 120)
        w = np.clip((mo - start) / (end - start), 0, 1)
        return w

    def _os_fn(self, mo: np.ndarray, arm: str) -> np.ndarray:
        # Non-waned base curves: compute directly (skip V2's waning logic)
        # We need GP-constrained base curves WITHOUT the old waning applied
        prefix = "combo" if arm == "combo" else "pembro"
        if self.p.os_distribution == "weibull":
            scale = getattr(self.p, f"os_weibull_scale_{prefix}")
            shape = getattr(self.p, f"os_weibull_shape_{prefix}")
            os_val = weibull_surv(mo, scale, shape)
        else:
            mu = getattr(self.p, f"os_mu_{prefix}")
            sigma = getattr(self.p, f"os_sigma_{prefix}")
            os_val = lognorm_surv(mo, mu, sigma)
        # Apply GP floor (same as V2)
        if self.p.constraint_general_pop:
            os_val = self._apply_gp_hazard_floor(os_val, mo)

        if not self.p.treatment_waning or arm != "combo":
            return os_val

        # ── Waning v2: hazard convergence for combo ──
        # Get pembro base (also GP-constrained)
        p_mu = getattr(self.p, "os_mu_pembro")
        p_sigma = getattr(self.p, "os_sigma_pembro")
        p_os = lognorm_surv(mo, p_mu, p_sigma)
        if self.p.constraint_general_pop:
            p_os = self._apply_gp_hazard_floor(p_os, mo)

        s_c = os_val  # combo base
        s_p = p_os    # pembro base

        # Monthly hazards
        h_c = np.zeros_like(s_c)
        h_p = np.zeros_like(s_p)
        n = len(s_c)
        for i in range(n - 1):
            h_c[i] = -np.log(max(s_c[i+1] / max(s_c[i], 1e-12), 1e-12))
            h_p[i] = -np.log(max(s_p[i+1] / max(s_p[i], 1e-12), 1e-12))
        h_c[-1] = h_c[-2]
        h_p[-1] = h_p[-2]

        # Waned hazard: combo hazard drifts to pembro hazard from 60mo onward
        w = self._waning_factor(mo)
        h_waned = (1 - w) * h_c + w * h_p

        # Reconstruct: keep base combo OS exactly before 60 months
        idx60 = int(60)
        s_waned = np.zeros_like(s_c)
        s_waned[:idx60+1] = s_c[:idx60+1]
        for i in range(idx60, n - 1):
            s_waned[i+1] = s_waned[i] * np.exp(-h_waned[i])
        for i in range(1, n):
            s_waned[i] = min(s_waned[i], s_waned[i-1])
        return s_waned

    def _survival(self) -> dict:
        """Same as parent but with waning v2 applied to r_dm and r_lr too."""
        mo = self.t * 12
        p = self.p
        sur = {}

        for arm in ["combo", "pembro"]:
            os_val = self._os_fn(mo, arm)
            sur[f"os_{arm}"] = os_val

            prefix = "combo" if arm == "combo" else "pembro"
            rd_mu = getattr(p, f"rd_mu_{prefix}")
            rd_sigma = getattr(p, f"rd_sigma_{prefix}")
            r_dm = lognorm_surv(mo, rd_mu, rd_sigma)
            rl_mu = getattr(p, f"rl_mu_{prefix}")
            rl_sigma = getattr(p, f"rl_sigma_{prefix}")
            r_lr = lognorm_surv(mo, rl_mu, rl_sigma)

            # Waning v2 for combo: r_dm / r_lr drift to pembro levels (HR->1) from 60mo
            if self.p.treatment_waning and arm == "combo":
                # Get pembro ratios
                prd_mu = getattr(p, f"rd_mu_pembro")
                prd_sigma = getattr(p, f"rd_sigma_pembro")
                p_r_dm = lognorm_surv(mo, prd_mu, prd_sigma)
                prl_mu = getattr(p, f"rl_mu_pembro")
                prl_sigma = getattr(p, f"rl_sigma_pembro")
                p_r_lr = lognorm_surv(mo, prl_mu, prl_sigma)
                w = self._waning_factor(mo)
                # Splice: keep base ratios before 60mo, drift after
                w_dm = (1 - w) * r_dm + w * p_r_dm
                w_lr = (1 - w) * r_lr + w * p_r_lr
                # Ensure ratios stay <= 1 and monotone non-increasing
                w_dm = np.minimum(w_dm, 1.0)
                w_lr = np.minimum(w_lr, 1.0)
                for i in range(1, len(w_dm)):
                    w_dm[i] = min(w_dm[i], w_dm[i-1])
                    w_lr[i] = min(w_lr[i], w_lr[i-1])
                r_dm = w_dm
                r_lr = w_lr

            sur[f"dmfs_{arm}"] = os_val * r_dm
            sur[f"rfs_{arm}"] = sur[f"dmfs_{arm}"] * r_lr

        return sur


def run_v3_diagnostics(p: ModelParams, label: str):
    m = CEAModelV3(p)
    t = m.t
    sur = m._survival()
    states = m._state_proportions(sur)
    res = m.run()

    print(f"\n{'='*90}")
    print(f"WANING V2: {label}")
    print(f"{'='*90}")

    # Structure
    viol_rfs_dmfs = viol_dmfs_os = 0
    maxsum = 0.0
    for i in range(len(t)):
        for a in ["combo", "pembro"]:
            if sur[f"rfs_{a}"][i] > sur[f"dmfs_{a}"][i] + 1e-10:
                viol_rfs_dmfs += 1
            if sur[f"dmfs_{a}"][i] > sur[f"os_{a}"][i] + 1e-10:
                viol_dmfs_os += 1
            st = states[a]
            ssum = st["rf"][i] + st["lr"][i] + st["dm"][i] + st["dead"][i]
            maxsum = max(maxsum, abs(ssum - 1.0))
    ok = viol_rfs_dmfs == 0 and viol_dmfs_os == 0 and maxsum < 1e-10
    print(f"RFS<=DMFS: {viol_rfs_dmfs}/480  DMFS<=OS: {viol_dmfs_os}/480  maxsum: {maxsum:.2e} → {'PASS' if ok else 'FAIL'}")

    # OS at key years
    print(f"\n{'Year':>5s} {'Combo OS':>10s} {'Pembro OS':>10s} {'Diff':>10s}")
    for yr in [5, 10, 15, 20, 30, 40]:
        idx = int(yr / p.cycle_length)
        if idx < len(t):
            print(f"{yr:5d} {sur['os_combo'][idx]:>10.4f} {sur['os_pembro'][idx]:>10.4f} {sur['os_combo'][idx]-sur['os_pembro'][idx]:>10.4f}")

    # Cost decomposition
    disc = np.exp(-p.discount_rate * t)
    cl = p.cycle_length
    comps = {}
    for a in ["combo", "pembro"]:
        s = states[a]
        comps[a] = {
            "Intismeran": (p.cost_intismeran * disc[0] * cl * 12) if a == "combo" else 0.0,
            "LR mgmt": np.sum(s["lr"] * p.cost_lr_monthly * 12 * cl * disc),
            "DM mgmt": np.sum(s["dm"] * p.cost_dm_monthly * 12 * cl * disc),
        }
    comps["combo"]["Pembro"] = p.cost_keytruda_annual * disc[:12].sum() * cl
    comps["pembro"]["Pembro"] = p.cost_keytruda_annual * disc[:12].sum() * cl

    print(f"\n{'Component':<15s} {'Combo':>12s} {'Pembro':>12s} {'Incr':>12s}")
    for c in ["Intismeran", "Pembro", "LR mgmt", "DM mgmt"]:
        cv = comps["combo"][c]; pv = comps["pembro"][c]
        print(f"{c:<15s} ${cv:>9,.0f} ${pv:>9,.0f} ${cv-pv:>9,.0f}")

    dc = res["combo"]["cost"] - res["pembro"]["cost"]
    dq = res["combo"]["qaly"] - res["pembro"]["qaly"]
    icer = dc / dq if dq > 0 else float('inf')
    print(f"\nCombo QALY={res['combo']['qaly']:.2f} Pembro QALY={res['pembro']['qaly']:.2f} ΔQALY={dq:.2f}")
    print(f"Combo Cost=${res['combo']['cost']:,.0f} Pembro=${res['pembro']['cost']:,.0f} ΔCost=${dc:,.0f}")
    print(f"ICER = ${icer:,.0f}/QALY")

    return {"icer": icer, "dq": dq, "dc": dc, "ok": ok}


if __name__ == "__main__":
    # Compare: base case vs waning v2 vs OLD waning (v1, for reference)
    print("=" * 90)
    print("BASE CASE (GP constrained) — reference")
    print("=" * 90)
    p = ModelParams(constraint_general_pop=True)
    m = CEAModelV2(p)
    r = m.run()
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    print(f"ΔQALY={dq:.2f} ΔCost=${dc:,.0f} ICER=${m.icer(r):,.0f}/QALY")

    # Waning v2 (5-10y): default window
    p1 = ModelParams(constraint_general_pop=True, treatment_waning=True)
    p1.wane_start_month = 60
    p1.wane_end_month = 120
    r1 = run_v3_diagnostics(p1, "Waning v2 (hazard convergence, 5-10y linear)")

    # Sensitivity: wane 5-15 years
    p2 = ModelParams(constraint_general_pop=True, treatment_waning=True)
    p2.wane_start_month = 60
    p2.wane_end_month = 180
    r2 = run_v3_diagnostics(p2, "Waning v2 extended (5-15y linear)")

    print(f"\n{'='*90}")
    print("COMPARISON")
    print(f"{'='*90}")
    print(f"Base case:              ΔQALY=3.76 ΔCost=$146,961 ICER=$39,105")
    print(f"Waning v2 (5-10y):      ΔQALY={r1['dq']:.2f} ΔCost=${r1['dc']:,.0f} ICER=${r1['icer']:,.0f}")
    print(f"Waning v2 (5-15y):      ΔQALY={r2['dq']:.2f} ΔCost=${r2['dc']:,.0f} ICER=${r2['icer']:,.0f}")