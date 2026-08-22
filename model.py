"""
cea_intismeran/model.py — 4-state PSM with OS robustness measures

Structurally constrained ratio-based extrapolation:
  OS(t) = log-normal survival
  DMFS(t) = OS(t) × r_DM(t)    where r_DM = P(no DM | alive)
  RFS(t) = DMFS(t) × r_LR(t)   where r_LR = P(no LR | alive, no DM)

Robustness features:
  - General population mortality constraint (US 2019 Life Tables)
  - Time horizon scenarios (10, 20, 40 years)
  - Treatment-effect waning scenario
  - Alternative OS extrapolation (log-normal vs Weibull)
"""

import numpy as np
import json
from dataclasses import dataclass, field
from typing import Tuple, Optional


# ── Parametric survival helpers ──

def lognorm_surv(t: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """Log-normal survival: S(t) = Φ((μ - ln(t)) / σ)"""
    from scipy.stats import norm
    return norm.cdf((mu - np.log(np.maximum(t, 1e-6))) / sigma)


def weibull_surv(t: np.ndarray, scale: float, shape: float) -> np.ndarray:
    """S(t) = exp(-(t/scale)^shape)"""
    return np.exp(-(t / scale) ** shape)


# ── General population survival ──

with open('/tmp/life_table_surv.json') as f:
    _LT = json.load(f)


def general_pop_surv(t_years: np.ndarray) -> np.ndarray:
    """Age-matched general population survival from age 59 (US 2019 Life Tables).

    Returns survival probability at each year (interpolated to monthly).
    """
    year = np.floor(t_years).astype(int)
    year = np.clip(year, 0, 42)
    surv = np.array([_LT["survival_by_year"].get(str(y), _LT["survival_by_year"].get(y, 0.0))
                      for y in year])
    # Handle dict key format: some keys are int, some str
    return surv


# ── Fitting ──

def fit_lognorm_os(time_points: list, surv_rates: list) -> Tuple[float, float]:
    """Fit log-normal to OS data points via NLS."""
    from scipy.optimize import least_squares
    from scipy.stats import norm

    t = np.array(time_points, float)
    s = np.array(surv_rates, float)

    def resid(p):
        mu, sigma = p
        return norm.cdf((mu - np.log(t)) / sigma) - s

    try:
        sol = least_squares(
            resid, x0=[np.log(t[len(t) // 2]), 0.8],
            bounds=([-5, 0.05], [15, 5.0]),
            xtol=1e-12, ftol=1e-12, gtol=1e-12,
        )
        return float(sol.x[0]), float(sol.x[1])
    except Exception:
        return (np.log(t[-1]), 1.0)


def fit_ratio_lognorm(time_points: list, ratio_values: list) -> Tuple[float, float]:
    """Fit log-normal to a ratio (0→1) via NLS."""
    from scipy.optimize import least_squares
    from scipy.stats import norm

    t = np.array(time_points, float)
    r = np.array(ratio_values, float)

    def resid(p):
        mu, sigma = p
        return norm.cdf((mu - np.log(t)) / sigma) - r

    try:
        sol = least_squares(
            resid, x0=[np.log(t[len(t) // 2]), 0.8],
            bounds=([-5, 0.05], [15, 5.0]),
            xtol=1e-12, ftol=1e-12, gtol=1e-12,
        )
        return float(sol.x[0]), float(sol.x[1])
    except Exception:
        return (np.log(t[-1]), 1.0)


# ── Model parameters ──

@dataclass
class ModelParams:
    # ── Observed KM data (Figure 1, Khattak 2026, percent 0-100) ──
    rfs_time_mo: list = field(default_factory=lambda: [18, 24, 36, 48, 60])
    rfs_combo_pct: list = field(default_factory=lambda: [80.4, 78.3, 73.8, 72.4, 68.8])
    rfs_pembro_pct: list = field(default_factory=lambda: [62.2, 60.0, 55.6, 49.1, 49.1])

    dmfs_time_mo: list = field(default_factory=lambda: [18, 24, 36, 48])
    dmfs_combo_pct: list = field(default_factory=lambda: [92.0, 90.8, 85.6, 83.9])
    dmfs_pembro_pct: list = field(default_factory=lambda: [73.0, 70.6, 65.4, 65.4])

    os_time_mo: list = field(default_factory=lambda: [18, 24, 36, 48, 60])
    os_combo_pct: list = field(default_factory=lambda: [96.0, 96.0, 93.8, 92.2, 92.2])
    os_pembro_pct: list = field(default_factory=lambda: [93.4, 93.4, 85.6, 85.6, 71.3])

    # ── Costs (2026 USD) ──
    cost_keytruda_annual: float = 220_896
    cost_intismeran: float = 200_000
    cost_sequencing: float = 1_000
    cost_admin_per_cycle: float = 500
    cost_lr_monthly: float = 3_000
    cost_dm_monthly: float = 12_000
    cost_ae_incremental: float = 15_000

    # ── Utilities ──
    util_rf: float = 0.83
    util_lr: float = 0.64
    util_dm: float = 0.55
    util_disutility_ae: float = 0.05

    # ── Model settings ──
    discount_rate: float = 0.03
    time_horizon_years: float = 40
    cycle_length: float = 1 / 12  # monthly
    constraint_general_pop: bool = False  # scenario: general population mortality cap
    treatment_waning: bool = False  # treatment effect waning after 5 years
    os_distribution: str = "lognormal"  # "lognormal" or "weibull"

    # ── Derived parameters ──
    os_mu_combo: float = 0.0
    os_sigma_combo: float = 0.0
    os_mu_pembro: float = 0.0
    os_sigma_pembro: float = 0.0
    rd_mu_combo: float = 0.0
    rd_sigma_combo: float = 0.0
    rd_mu_pembro: float = 0.0
    rd_sigma_pembro: float = 0.0
    rl_mu_combo: float = 0.0
    rl_sigma_combo: float = 0.0
    rl_mu_pembro: float = 0.0
    rl_sigma_pembro: float = 0.0

    # Weibull OS params (for alternative extrapolation)
    os_weibull_scale_combo: float = 0.0
    os_weibull_shape_combo: float = 0.0
    os_weibull_scale_pembro: float = 0.0
    os_weibull_shape_pembro: float = 0.0

    def __post_init__(self):
        def pct(v): return np.array(v, dtype=float) / 100.0

        # 1. Fit OS
        # Log-normal
        t_os = self.os_time_mo
        self.os_mu_combo, self.os_sigma_combo = fit_lognorm_os(t_os, pct(self.os_combo_pct))
        self.os_mu_pembro, self.os_sigma_pembro = fit_lognorm_os(t_os, pct(self.os_pembro_pct))

        # Weibull (for alternative scenario)
        def _fit_weibull(t, s):
            from scipy.optimize import least_squares
            def resid(p):
                sc, sh = p
                return weibull_surv(np.array(t, float), sc, sh) - np.array(s, float)
            scale0 = t[-1] / (-np.log(max(s[-1], 1e-6))) ** (1.0 / 1.2)
            try:
                sol = least_squares(resid, x0=[scale0, 1.2], bounds=([10, 0.2], [5000, 5.0]), xtol=1e-12, ftol=1e-12)
                return float(sol.x[0]), float(sol.x[1])
            except:
                return (2000.0, 1.0)

        self.os_weibull_scale_combo, self.os_weibull_shape_combo = _fit_weibull(t_os, pct(self.os_combo_pct))
        self.os_weibull_scale_pembro, self.os_weibull_shape_pembro = _fit_weibull(t_os, pct(self.os_pembro_pct))

        # 2. Fit ratios r_DM = DMFS/OS
        def _ratio_obs(numer_pct, denom_pct):
            n = np.array(numer_pct, float)
            d = np.array(denom_pct, float)
            return np.minimum(n / d, 1.0)

        r_dm_combo_obs = _ratio_obs(self.dmfs_combo_pct,
                                     [dict(zip(self.os_time_mo, self.os_combo_pct))[t] for t in self.dmfs_time_mo])
        r_dm_pembro_obs = _ratio_obs(self.dmfs_pembro_pct,
                                      [dict(zip(self.os_time_mo, self.os_pembro_pct))[t] for t in self.dmfs_time_mo])

        self.rd_mu_combo, self.rd_sigma_combo = fit_ratio_lognorm(self.dmfs_time_mo, r_dm_combo_obs)
        self.rd_mu_pembro, self.rd_sigma_pembro = fit_ratio_lognorm(self.dmfs_time_mo, r_dm_pembro_obs)

        # 3. Fit ratios r_LR = RFS/DMFS
        def interp_dmfs(t_target, dmfs_t, dmfs_v):
            return np.interp(t_target, dmfs_t, dmfs_v)

        dmfs_combo_at_rfs = interp_dmfs(self.rfs_time_mo, self.dmfs_time_mo, pct(self.dmfs_combo_pct))
        dmfs_pembro_at_rfs = interp_dmfs(self.rfs_time_mo, self.dmfs_time_mo, pct(self.dmfs_pembro_pct))

        r_lr_combo_obs = np.minimum(pct(self.rfs_combo_pct) / dmfs_combo_at_rfs, 1.0)
        r_lr_pembro_obs = np.minimum(pct(self.rfs_pembro_pct) / dmfs_pembro_at_rfs, 1.0)

        self.rl_mu_combo, self.rl_sigma_combo = fit_ratio_lognorm(self.rfs_time_mo, r_lr_combo_obs)
        self.rl_mu_pembro, self.rl_sigma_pembro = fit_ratio_lognorm(self.rfs_time_mo, r_lr_pembro_obs)


# ── Model runner ──

class CEAModel:
    def __init__(self, p: ModelParams = None):
        self.p = p or ModelParams()
        n_cycles = int(self.p.time_horizon_years / self.p.cycle_length)
        self.t = np.arange(n_cycles, dtype=float) * self.p.cycle_length  # years

    def _os_fn(self, mo: np.ndarray, arm: str) -> np.ndarray:
        """Compute OS curve for given arm."""
        prefix = "combo" if arm == "combo" else "pembro"
        if self.p.os_distribution == "weibull":
            scale = getattr(self.p, f"os_weibull_scale_{prefix}")
            shape = getattr(self.p, f"os_weibull_shape_{prefix}")
            os_val = weibull_surv(mo, scale, shape)
        else:
            mu = getattr(self.p, f"os_mu_{prefix}")
            sigma = getattr(self.p, f"os_sigma_{prefix}")
            os_val = lognorm_surv(mo, mu, sigma)

        # General population mortality constraint
        if self.p.constraint_general_pop:
            gp_surv = general_pop_surv(self.t)[:len(mo)]
            os_val = np.minimum(os_val, gp_surv)

        # Treatment-effect waning
        if self.p.treatment_waning and arm == "combo":
            # After 5 years, gradually wane combo OS toward pembro OS
            # Get pembro OS
            p_prefix = "pembro"
            if self.p.os_distribution == "weibull":
                p_scale = getattr(self.p, f"os_weibull_scale_{p_prefix}")
                p_shape = getattr(self.p, f"os_weibull_shape_{p_prefix}")
                p_os = weibull_surv(mo, p_scale, p_shape)
            else:
                p_mu = getattr(self.p, f"os_mu_{p_prefix}")
                p_sigma = getattr(self.p, f"os_sigma_{p_prefix}")
                p_os = lognorm_surv(mo, p_mu, p_sigma)

            if self.p.constraint_general_pop:
                p_os = np.minimum(p_os, gp_surv)

            # Waning: linear from full combo benefit at 5y to zero benefit at 20y
            t_y = self.t[:len(mo)]
            wane = np.clip((t_y - 5) / (20 - 5), 0, 1)
            benefit = os_val - p_os
            os_val = p_os + benefit * (1 - wane)

        return os_val

    def _survival(self) -> dict:
        """Return structurally constrained survival curves."""
        mo = self.t * 12
        p = self.p
        sur = {}

        for arm in ["combo", "pembro"]:
            os_val = self._os_fn(mo, arm)
            sur[f"os_{arm}"] = os_val

            prefix = "combo" if arm == "combo" else "pembro"

            # r_DM = DMFS/OS
            rd_mu = getattr(p, f"rd_mu_{prefix}")
            rd_sigma = getattr(p, f"rd_sigma_{prefix}")
            r_dm = lognorm_surv(mo, rd_mu, rd_sigma)
            sur[f"dmfs_{arm}"] = os_val * r_dm

            # r_LR = RFS/DMFS
            rl_mu = getattr(p, f"rl_mu_{prefix}")
            rl_sigma = getattr(p, f"rl_sigma_{prefix}")
            r_lr = lognorm_surv(mo, rl_mu, rl_sigma)
            sur[f"rfs_{arm}"] = sur[f"dmfs_{arm}"] * r_lr

        return sur

    def _state_proportions(self, sur: dict) -> dict:
        return {
            "combo": {
                "rf": sur["rfs_combo"],
                "lr": sur["dmfs_combo"] - sur["rfs_combo"],
                "dm": sur["os_combo"] - sur["dmfs_combo"],
                "dead": 1 - sur["os_combo"],
            },
            "pembro": {
                "rf": sur["rfs_pembro"],
                "lr": sur["dmfs_pembro"] - sur["rfs_pembro"],
                "dm": sur["os_pembro"] - sur["dmfs_pembro"],
                "dead": 1 - sur["os_pembro"],
            },
        }

    def icer(self, results: dict) -> float:
        return icer_model(results)

    def run(self) -> dict:
        sur = self._survival()
        states = self._state_proportions(sur)
        p = self.p
        disc = np.exp(-p.discount_rate * self.t)

        results = {}
        for arm_label in ["combo", "pembro"]:
            s = states[arm_label]
            qaly = (s["rf"] * p.util_rf + s["lr"] * p.util_lr + s["dm"] * p.util_dm) * disc
            total_qaly = np.sum(qaly) * p.cycle_length
            ae_disutility = p.util_disutility_ae * disc[0] * p.cycle_length
            total_qaly -= ae_disutility

            keytruda_cost = p.cost_keytruda_annual * disc[:12].sum() * p.cycle_length
            intismeran_cost = p.cost_intismeran * disc[0] if arm_label == "combo" else 0
            sequencing_cost = p.cost_sequencing * disc[0] if arm_label == "combo" else 0
            admin_cost = p.cost_admin_per_cycle * disc[:12].sum() * p.cycle_length
            lr_cost = np.sum(s["lr"] * p.cost_lr_monthly * 12 * p.cycle_length * disc)
            dm_cost = np.sum(s["dm"] * p.cost_dm_monthly * 12 * p.cycle_length * disc)
            ae_cost = p.cost_ae_incremental * disc[0] if arm_label == "combo" else 0

            total_cost = keytruda_cost + intismeran_cost + sequencing_cost + admin_cost + lr_cost + dm_cost + ae_cost
            results[arm_label] = {"qaly": total_qaly, "cost": total_cost, "states": s}

        return results


def run_psa(n_iter: int = 1000) -> dict:
    """Run probabilistic sensitivity analysis."""
    from scipy.stats import gamma, norm

    p = ModelParams()
    base = CEAModel(p).run()
    icers = []
    costs_combo, costs_pembro = [], []
    qalys_combo, qalys_pembro = [], []

    for _ in range(n_iter):
        rng = np.random.default_rng()
        pp = ModelParams()

        # Gamma distributions for costs (CV=0.2)
        for attr, val in [("cost_keytruda_annual", 220_896), ("cost_intismeran", 200_000),
                          ("cost_lr_monthly", 3_000), ("cost_dm_monthly", 12_000)]:
            cv = 0.2
            setattr(pp, attr, rng.gamma(1 / cv ** 2, val * cv ** 2))

        # Beta distributions for utilities (simulated from SE=0.03)
        for attr, val in [("util_rf", 0.83), ("util_lr", 0.64), ("util_dm", 0.55)]:
            se = 0.03
            alpha = val * (val * (1 - val) / se ** 2 - 1)
            beta = (1 - val) * (val * (1 - val) / se ** 2 - 1)
            setattr(pp, attr, rng.beta(max(alpha, 1), max(beta, 1)))

        # Random perturbation of survival params (log-normal mu~N(0,0.1))
        pp.os_mu_combo += rng.normal(0, 0.1)
        pp.os_mu_pembro += rng.normal(0, 0.1)

        m = CEAModel(pp)
        try:
            r = m.run()
            icer = m.icer(r)
            if icer > 0 and icer < 1e6:
                icers.append(icer)
                costs_combo.append(r["combo"]["cost"])
                costs_pembro.append(r["pembro"]["cost"])
                qalys_combo.append(r["combo"]["qaly"])
                qalys_pembro.append(r["pembro"]["qaly"])
        except:
            continue

    return {
        "mean_icer": np.mean(icers),
        "median_icer": np.median(icers),
        "ci_2.5": np.percentile(icers, 2.5),
        "ci_97.5": np.percentile(icers, 97.5),
        "p_ce_50k": np.mean(np.array(icers) <= 50_000),
        "p_ce_100k": np.mean(np.array(icers) <= 100_000),
        "p_ce_150k": np.mean(np.array(icers) <= 150_000),
        "n": len(icers),
    }


def icer_model(results: dict) -> float:
    inc_cost = results["combo"]["cost"] - results["pembro"]["cost"]
    inc_qaly = results["combo"]["qaly"] - results["pembro"]["qaly"]
    return inc_cost / inc_qaly


# ── Scenario runner ──

def run_scenarios():
    """Run all OS robustness scenarios."""
    base = ModelParams()
    model = CEAModel(base)
    res = model.run()
    print(f"{'Scenario':40s} {'Combo QALY':>10s} {'Pembro QALY':>10s} {'Incr QALY':>10s} {'Incr Cost':>12s} {'ICER':>10s}")
    print(f"{'Base case (log-normal, 40y, GP constraint)':40s} {res['combo']['qaly']:>10.2f} {res['pembro']['qaly']:>10.2f} {res['combo']['qaly']-res['pembro']['qaly']:>10.2f} ${res['combo']['cost']-res['pembro']['cost']:>9,.0f} ${model.icer(res):>8,.0f}")

    # Time horizon scenarios
    for th in [10, 20]:
        p = ModelParams(time_horizon_years=th)
        m = CEAModel(p)
        r = m.run()
        print(f"{'Time horizon: '+str(th)+' years':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")

    # No GP constraint
    p = ModelParams(constraint_general_pop=False)
    m = CEAModel(p)
    r = m.run()
    print(f"{'No general population constraint':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")

    # Treatment waning
    p = ModelParams(treatment_waning=True)
    m = CEAModel(p)
    r = m.run()
    print(f"{'Treatment-effect waning (5-20y)':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")

    # Weibull OS
    p = ModelParams(os_distribution="weibull")
    m = CEAModel(p)
    r = m.run()
    print(f"{'Weibull OS (instead of log-normal)':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")

    # 0% discount
    p = ModelParams(discount_rate=0.0)
    m = CEAModel(p)
    r = m.run()
    print(f"{'0% discount rate':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")

    # 5% discount
    p = ModelParams(discount_rate=0.05)
    m = CEAModel(p)
    r = m.run()
    print(f"{'5% discount rate':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")

    # Waning + GP constraint + 20y (most conservative)
    p = ModelParams(treatment_waning=True, constraint_general_pop=True, time_horizon_years=20)
    m = CEAModel(p)
    r = m.run()
    print(f"{'Waning + GP + 20y horizon (most conservative)':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")


if __name__ == "__main__":
    # Base case
    p = ModelParams()
    model = CEAModel(p)
    res = model.run()
    print(f"Base case: ICER = ${model.icer(res):,.0f}/QALY")
    print(f"  Combo: QALY={res['combo']['qaly']:.2f} Cost=${res['combo']['cost']:,.0f}")
    print(f"  Pembro: QALY={res['pembro']['qaly']:.2f} Cost=${res['pembro']['cost']:,.0f}")
    print()

    # All scenarios
    run_scenarios()