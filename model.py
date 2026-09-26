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
import os
from dataclasses import dataclass, field
from typing import Tuple, Optional


# ── Regimen definitions (Methods: pembrolizumab 200 mg every 3 weeks, up to 18
#    cycles; intismeran autogene 9 doses on the same 3-week schedule; $500 per
#    infusion) ──
N_PEMBRO_DOSES = 18
N_INTISMERAN_DOSES = 9
PEMBRO_CYCLE_WEEKS = 3



def discount_factors(rate: float, t) -> np.ndarray:
    """Discrete annual discount factor (1+r)^-t, t in years (Second Panel reference case)."""
    return (1.0 + rate) ** (-np.asarray(t, dtype=float))


# ── Parametric survival helpers ──

def lognorm_surv(t: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """Log-normal survival: S(t) = Φ((μ - ln(t)) / σ)"""
    from scipy.stats import norm
    return norm.cdf((mu - np.log(np.maximum(t, 1e-6))) / sigma)


def weibull_surv(t: np.ndarray, scale: float, shape: float) -> np.ndarray:
    """S(t) = exp(-(t/scale)^shape)"""
    return np.exp(-(t / scale) ** shape)


def loglogistic_surv(t: np.ndarray, shape: float, scale: float) -> np.ndarray:
    """Log-logistic survival: S(t) = 1 / (1 + (t/scale)^shape)"""
    return 1.0 / (1.0 + (np.maximum(t, 1e-6) / scale) ** shape)


# ── General population survival ──
# Repo copy first (data/life_table_surv.json); /tmp kept as a legacy fallback.
_LT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "data", "life_table_surv.json")
if not os.path.exists(_LT_PATH):
    _LT_PATH = "/tmp/life_table_surv.json"

with open(_LT_PATH) as f:
    _LT = json.load(f)


def general_pop_surv(t_years: np.ndarray) -> np.ndarray:
    """Age-matched general population survival from baseline age 61 (US 2019 Life Tables).

    Returns log-linear interpolated survival at requested years (floats allowed).
    """
    # Build annual survival array (year 0..42)
    surv_year = np.array([_LT["survival_by_year"].get(str(y), 0.0) for y in range(0, 43)])
    # Add year 43 as ~0 to bound interpolation beyond 42
    x = np.arange(0, 43, dtype=float)
    y = np.log(np.maximum(surv_year, 1e-12))
    log_surv = np.interp(np.asarray(t_years, dtype=float), x, y)
    return np.exp(log_surv)


def survival_at_month(surv: np.ndarray, month: float) -> float:
    """S(month) from a monthly survival array whose index i is month i.

    Single entry point for every landmark read (Supplementary Tables S5-S9, all
    figures): the array index equals the month number, so a landmark labelled
    "60 months" reads index 60. Values beyond the array are clamped to the last
    point (the model's projection horizon).
    """
    m = int(round(float(month)))
    if m < 0:
        raise ValueError(f"month must be >= 0, got {month!r}")
    if len(surv) == 0:
        raise ValueError("empty survival array")
    return float(surv[min(m, len(surv) - 1)])


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
    # The two CMS fee-schedule amounts below are 2025 rates inflated to 2026 with the
    # same medical-care CPI factor used for the other price years, so every cost in the
    # model is expressed in 2026 US dollars (576.029 -> 592.116, factor 1.027926).
    # Values published in earlier years were inflated to 2026 US dollars with the
    # medical-care component of the Consumer Price Index (Methods, Costs section).
    cost_keytruda_annual: float = 220_896
    cost_intismeran: float = 200_000
    cost_sequencing: float = 1_000
    # Administration, by route (CMS Physician Fee Schedule, CPT codes below).
    # KEYNOTE-942: intismeran (mRNA-4157) 1 mg intramuscularly, maximum nine doses;
    # pembrolizumab 200 mg intravenously, maximum 18 doses; both every 3 weeks.
    cost_admin_iv_infusion: float = 122.62     # per IV infusion (pembrolizumab, CPT 96413)
    cost_admin_im_injection: float = 14.3     # per IM injection (intismeran; CPT 96372 proxy)
    cost_lr_monthly: float = 3_000
    cost_dm_monthly: float = 12_000
    cost_ae_incremental: float = 15_000

    # ── Utilities ──
    util_rf: float = 0.83
    util_lr: float = 0.64
    util_dm: float = 0.55

    # ── Model settings ──
    discount_rate: float = 0.03
    time_horizon_years: float = 40
    cycle_length: float = 1 / 12  # monthly
    constraint_general_pop: bool = False  # scenario: general population mortality cap
    gp_floor_start_months: int = 60  # months after which GP floor applies
    treatment_waning: bool = False  # treatment effect waning after 5 years
    os_hr_scaling: float = 1.0  # Scale OS hazard ratio (1.0 = base HR 0.471; 0.35 = lower bound HR 0.165; 2.86 = upper bound HR 1.345)
    os_distribution: str = "lognormal"  # "lognormal" or "weibull" (applies to both arms unless per-arm overrides set)
    os_distribution_combo: str = ""  # if non-empty, overrides os_distribution for combo arm
    os_distribution_pembro: str = ""  # if non-empty, overrides os_distribution for pembro arm
    rDM_distribution: str = "lognormal"  # "lognormal", "weibull", "loglogistic", "gengamma" (default for both arms)
    rDM_distribution_combo: str = ""  # if non-empty, overrides rDM_distribution for combo arm
    rDM_distribution_pembro: str = ""  # if non-empty, overrides rDM_distribution for pembro arm
    rLR_distribution: str = "lognormal"  # "lognormal", "weibull", "loglogistic", "gengamma" (default for both arms)
    rLR_distribution_combo: str = ""  # if non-empty, overrides rLR_distribution for combo arm
    rLR_distribution_pembro: str = ""  # if non-empty, overrides rLR_distribution for pembro arm

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

    # Weibull rDM params
    rd_weibull_scale_combo: float = 0.0
    rd_weibull_shape_combo: float = 0.0
    rd_weibull_scale_pembro: float = 0.0
    rd_weibull_shape_pembro: float = 0.0

    # Weibull rLR params
    rl_weibull_scale_combo: float = 0.0
    rl_weibull_shape_combo: float = 0.0
    rl_weibull_scale_pembro: float = 0.0
    rl_weibull_shape_pembro: float = 0.0

    # Log-logistic rDM params
    rd_loglogistic_shape_combo: float = 0.0
    rd_loglogistic_scale_combo: float = 0.0
    rd_loglogistic_shape_pembro: float = 0.0
    rd_loglogistic_scale_pembro: float = 0.0

    # Log-logistic rLR params
    rl_loglogistic_shape_combo: float = 0.0
    rl_loglogistic_scale_combo: float = 0.0
    rl_loglogistic_shape_pembro: float = 0.0
    rl_loglogistic_scale_pembro: float = 0.0

    # GenGamma rDM params
    rd_gengamma_a_combo: float = 0.0
    rd_gengamma_c_combo: float = 0.0
    rd_gengamma_b_combo: float = 0.0
    rd_gengamma_a_pembro: float = 0.0
    rd_gengamma_c_pembro: float = 0.0
    rd_gengamma_b_pembro: float = 0.0

    # GenGamma rLR params
    rl_gengamma_a_combo: float = 0.0
    rl_gengamma_c_combo: float = 0.0
    rl_gengamma_b_combo: float = 0.0
    rl_gengamma_a_pembro: float = 0.0
    rl_gengamma_c_pembro: float = 0.0
    rl_gengamma_b_pembro: float = 0.0

    def __post_init__(self):
        def pct(v): return np.array(v, dtype=float) / 100.0

        # Adjust OS combo data for HR sensitivity (before fitting)
        # Base HR = 0.471; target HR = base_HR * os_hr_scaling
        os_combo_pct_adj = list(self.os_combo_pct)  # copy
        if abs(self.os_hr_scaling - 1.0) > 0.01:
            # Convert pembro survival to hazard, apply scaled HR to get combo hazard, back to survival
            base_hr = 0.471
            target_hr = base_hr * self.os_hr_scaling
            # Use log-normal approximation: if S_pembro(t) = Φ((μ_p - ln(t))/σ_p),
            # then hazard ratio ≈ exp( (μ_c - μ_p)/σ ) for same σ
            # We need to find μ_c such that HR = target_hr
            # S_c(t) = S_p(t) ^ (hazard_ratio) approx
            for i, t_mo in enumerate(self.os_time_mo):
                t_yr = t_mo / 12.0
                s_pembro = pct(self.os_pembro_pct)[i]
                # Approximate: S_combo = S_pembro ^ (target_HR / base_HR_effect)
                # But base case already has S_combo, so we scale relative to base
                s_combo_base = pct(self.os_combo_pct)[i]
                # New combo survival: scale hazard
                # h_new = h_base * (target_hr / 0.471)
                # S_new = exp(-H_new) = exp(-H_base * target_hr / 0.471) = S_base ^ (target_hr / 0.471)
                hr_ratio = target_hr / 0.471
                s_combo_new = s_combo_base ** hr_ratio
                os_combo_pct_adj[i] = s_combo_new * 100.0

        # 1. Fit OS
        # Log-normal
        t_os = self.os_time_mo
        self.os_mu_combo, self.os_sigma_combo = fit_lognorm_os(t_os, pct(os_combo_pct_adj))
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

        # Weibull rDM
        self.rd_weibull_scale_combo, self.rd_weibull_shape_combo = _fit_weibull(self.dmfs_time_mo, r_dm_combo_obs)
        self.rd_weibull_scale_pembro, self.rd_weibull_shape_pembro = _fit_weibull(self.dmfs_time_mo, r_dm_pembro_obs)

        # Log-logistic rDM
        def _fit_loglogistic(t, s):
            from scipy.optimize import least_squares
            def resid(p):
                sh, sc = p
                return loglogistic_surv(np.array(t, float), sh, sc) - np.array(s, float)
            try:
                sol = least_squares(resid, x0=[1.5, 50.0], bounds=([0.05, 1], [15, 5000]), xtol=1e-12, ftol=1e-12)
                return float(sol.x[0]), float(sol.x[1])
            except:
                return (1.5, 50.0)

        self.rd_loglogistic_shape_combo, self.rd_loglogistic_scale_combo = _fit_loglogistic(self.dmfs_time_mo, r_dm_combo_obs)
        self.rd_loglogistic_shape_pembro, self.rd_loglogistic_scale_pembro = _fit_loglogistic(self.dmfs_time_mo, r_dm_pembro_obs)

        # GenGamma rDM
        def _fit_gengamma(t, s):
            from scipy.optimize import least_squares
            from scipy.stats import gengamma
            def resid(p):
                a, c, b = p
                return gengamma.sf(np.array(t, float), a=np.exp(a), c=c, scale=np.exp(b)) - np.array(s, float)
            try:
                sol = least_squares(resid, x0=[0.5, 1.0, np.log(50)], bounds=([-2, 0.1, 1], [5, 5, 300]), xtol=1e-10, ftol=1e-10, max_nfev=5000)
                return float(np.exp(sol.x[0])), float(sol.x[1]), float(np.exp(sol.x[2]))
            except:
                return (1.0, 1.0, 50.0)

        self.rd_gengamma_a_combo, self.rd_gengamma_c_combo, self.rd_gengamma_b_combo = _fit_gengamma(self.dmfs_time_mo, r_dm_combo_obs)
        self.rd_gengamma_a_pembro, self.rd_gengamma_c_pembro, self.rd_gengamma_b_pembro = _fit_gengamma(self.dmfs_time_mo, r_dm_pembro_obs)

        # 3. Fit ratios r_LR = RFS/DMFS
        def interp_dmfs(t_target, dmfs_t, dmfs_v):
            return np.interp(t_target, dmfs_t, dmfs_v)

        dmfs_combo_at_rfs = interp_dmfs(self.rfs_time_mo, self.dmfs_time_mo, pct(self.dmfs_combo_pct))
        dmfs_pembro_at_rfs = interp_dmfs(self.rfs_time_mo, self.dmfs_time_mo, pct(self.dmfs_pembro_pct))

        r_lr_combo_obs = np.minimum(pct(self.rfs_combo_pct) / dmfs_combo_at_rfs, 1.0)
        r_lr_pembro_obs = np.minimum(pct(self.rfs_pembro_pct) / dmfs_pembro_at_rfs, 1.0)

        self.rl_mu_combo, self.rl_sigma_combo = fit_ratio_lognorm(self.rfs_time_mo, r_lr_combo_obs)
        self.rl_mu_pembro, self.rl_sigma_pembro = fit_ratio_lognorm(self.rfs_time_mo, r_lr_pembro_obs)

        # Weibull rLR
        self.rl_weibull_scale_combo, self.rl_weibull_shape_combo = _fit_weibull(self.rfs_time_mo, r_lr_combo_obs)
        self.rl_weibull_scale_pembro, self.rl_weibull_shape_pembro = _fit_weibull(self.rfs_time_mo, r_lr_pembro_obs)

        # Log-logistic rLR
        self.rl_loglogistic_shape_combo, self.rl_loglogistic_scale_combo = _fit_loglogistic(self.rfs_time_mo, r_lr_combo_obs)
        self.rl_loglogistic_shape_pembro, self.rl_loglogistic_scale_pembro = _fit_loglogistic(self.rfs_time_mo, r_lr_pembro_obs)

        # GenGamma rLR
        self.rl_gengamma_a_combo, self.rl_gengamma_c_combo, self.rl_gengamma_b_combo = _fit_gengamma(self.rfs_time_mo, r_lr_combo_obs)
        self.rl_gengamma_a_pembro, self.rl_gengamma_c_pembro, self.rl_gengamma_b_pembro = _fit_gengamma(self.rfs_time_mo, r_lr_pembro_obs)


# ── Model runner ──

class CEAModel:
    def __init__(self, p: ModelParams = None):
        self.p = p or ModelParams()
        n_cycles = int(self.p.time_horizon_years / self.p.cycle_length)
        self.t = np.arange(n_cycles, dtype=float) * self.p.cycle_length  # years

    def _os_fn(self, mo: np.ndarray, arm: str) -> np.ndarray:
        """Compute OS curve for given arm."""
        prefix = "combo" if arm == "combo" else "pembro"
        # Determine distribution for this arm
        arm_dist = getattr(self.p, f"os_distribution_{prefix}", "")
        if not arm_dist:
            arm_dist = self.p.os_distribution
        if arm_dist == "weibull":
            scale = getattr(self.p, f"os_weibull_scale_{prefix}")
            shape = getattr(self.p, f"os_weibull_shape_{prefix}")
            os_val = weibull_surv(mo, scale, shape)
        else:
            mu = getattr(self.p, f"os_mu_{prefix}")
            sigma = getattr(self.p, f"os_sigma_{prefix}")
            os_val = lognorm_surv(mo, mu, sigma)

        # General population mortality hazard floor (from 60 months)
        if self.p.constraint_general_pop:
            gp_surv = general_pop_surv(self.t)[:len(mo)]
            # Convert OS to hazard, floor at GP hazard, reconstruct
            n = len(os_val)
            h_model = np.zeros(n)
            for i in range(n - 1):
                h_model[i] = -np.log(max(os_val[i + 1] / max(os_val[i], 1e-12), 1e-12))
            mask = np.arange(n) >= self.p.gp_floor_start_months  # months
            gp_h = np.zeros(n)
            for i in range(n):
                gp_h[i] = -np.log(max(gp_surv[min(i+1, len(gp_surv)-1)] / max(gp_surv[i], 1e-12), 1e-12))
            h_final = np.where(mask, np.maximum(h_model, gp_h), h_model)
            h_final[-1] = h_final[-2]
            # Reconstruct S(i) = exp(-sum_{j<i} h_j): index i is month i.
            cum_h = np.concatenate([[0.0], np.cumsum(h_final[:-1])])
            os_val = np.exp(-cum_h)
            # Ensure monotonic
            for i in range(1, n):
                os_val[i] = min(os_val[i], os_val[i-1])

        # Treatment-effect waning
        if self.p.treatment_waning and arm == "combo":
            # After 5 years, gradually wane combo OS toward pembro OS
            # Get pembro OS
            p_prefix = "pembro"
            p_arm_dist = getattr(self.p, f"os_distribution_{p_prefix}", "")
            if not p_arm_dist:
                p_arm_dist = self.p.os_distribution
            if p_arm_dist == "weibull":
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

        # OS HR sensitivity (clinical efficacy uncertainty)
        # Scale hazard ratio: os_hr_scaling = test_HR / base_HR
        # base HR = 0.471; lower bound = 0.165 (scale = 0.35); upper bound = 1.345 (scale = 2.86)
        if arm == "combo" and abs(self.p.os_hr_scaling - 1.0) > 0.01:
            # Convert survival to hazard, scale, then back to survival
            n = len(os_val)
            h_model = np.zeros(n)
            for i in range(n - 1):
                h_model[i] = -np.log(max(os_val[i + 1] / max(os_val[i], 1e-12), 1e-12))
            # Apply HR scaling to combo arm hazard (higher scale = worse survival)
            h_model = h_model * self.p.os_hr_scaling
            # Reconstruct survival
            cum_h = np.concatenate([[0.0], np.cumsum(h_model[:-1])])
            os_val = np.exp(-cum_h)
            # Ensure monotonic
            for i in range(1, n):
                os_val[i] = min(os_val[i], os_val[i-1])

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
            rdm_dist = getattr(p, f"rDM_distribution_{prefix}", "")
            if not rdm_dist:
                rdm_dist = p.rDM_distribution
            if rdm_dist == "weibull":
                rd_sc = getattr(p, f"rd_weibull_scale_{prefix}")
                rd_sh = getattr(p, f"rd_weibull_shape_{prefix}")
                r_dm = weibull_surv(mo, rd_sc, rd_sh)
            elif rdm_dist == "loglogistic":
                rd_sh = getattr(p, f"rd_loglogistic_shape_{prefix}")
                rd_sc = getattr(p, f"rd_loglogistic_scale_{prefix}")
                r_dm = loglogistic_surv(mo, rd_sh, rd_sc)
            elif rdm_dist == "gengamma":
                from scipy.stats import gengamma as gengamma_dist
                rd_a = getattr(p, f"rd_gengamma_a_{prefix}")
                rd_c = getattr(p, f"rd_gengamma_c_{prefix}")
                rd_b = getattr(p, f"rd_gengamma_b_{prefix}")
                r_dm = gengamma_dist.sf(mo, a=rd_a, c=rd_c, scale=rd_b)
            else:
                rd_mu = getattr(p, f"rd_mu_{prefix}")
                rd_sigma = getattr(p, f"rd_sigma_{prefix}")
                r_dm = lognorm_surv(mo, rd_mu, rd_sigma)
            sur[f"dmfs_{arm}"] = os_val * r_dm

            # r_LR = RFS/DMFS
            rlr_dist = getattr(p, f"rLR_distribution_{prefix}", "")
            if not rlr_dist:
                rlr_dist = p.rLR_distribution
            if rlr_dist == "weibull":
                rl_sc = getattr(p, f"rl_weibull_scale_{prefix}")
                rl_sh = getattr(p, f"rl_weibull_shape_{prefix}")
                r_lr = weibull_surv(mo, rl_sc, rl_sh)
            elif rlr_dist == "loglogistic":
                rl_sh = getattr(p, f"rl_loglogistic_shape_{prefix}")
                rl_sc = getattr(p, f"rl_loglogistic_scale_{prefix}")
                r_lr = loglogistic_surv(mo, rl_sh, rl_sc)
            elif rlr_dist == "gengamma":
                from scipy.stats import gengamma as gengamma_dist
                rl_a = getattr(p, f"rl_gengamma_a_{prefix}")
                rl_c = getattr(p, f"rl_gengamma_c_{prefix}")
                rl_b = getattr(p, f"rl_gengamma_b_{prefix}")
                r_lr = gengamma_dist.sf(mo, a=rl_a, c=rl_c, scale=rl_b)
            else:
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
        disc = discount_factors(p.discount_rate, self.t)

        results = {}
        for arm_label in ["combo", "pembro"]:
            s = states[arm_label]
            qaly = (s["rf"] * p.util_rf + s["lr"] * p.util_lr + s["dm"] * p.util_dm) * disc
            total_qaly = np.sum(qaly) * p.cycle_length

            keytruda_cost = p.cost_keytruda_annual * disc[:12].sum() * p.cycle_length
            intismeran_cost = p.cost_intismeran * disc[0] if arm_label == "combo" else 0
            sequencing_cost = p.cost_sequencing * disc[0] if arm_label == "combo" else 0
            # Administration, charged by route: the 18 pembrolizumab doses in both
            # arms are intravenous infusions (CPT 96365); the combination arm's 9
            # intismeran doses are intramuscular injections (CPT 96372, 1 mg) and are
            # charged at the injection rate, not the infusion rate.
            pembro_t = np.arange(N_PEMBRO_DOSES) * (PEMBRO_CYCLE_WEEKS / 52)
            admin_cost = p.cost_admin_iv_infusion * float(
                discount_factors(p.discount_rate, pembro_t).sum())
            if arm_label == "combo":
                ismeran_t = np.arange(N_INTISMERAN_DOSES) * (PEMBRO_CYCLE_WEEKS / 52)
                admin_cost += p.cost_admin_im_injection * float(
                    discount_factors(p.discount_rate, ismeran_t).sum())
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

        # Gamma distributions for costs (CV=0.2) - intismeran held fixed per Methods
        for attr, val in [("cost_keytruda_annual", 220_896),
                          ("cost_lr_monthly", 3_000), ("cost_dm_monthly", 12_000)]:
            cv = 0.2
            setattr(pp, attr, rng.gamma(1 / cv ** 2, val * cv ** 2))
        # Intismeran price held fixed at $200,000 (price uncertainty explored in threshold analysis)

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

    # OS HR sensitivity (clinical efficacy uncertainty bounds)
    # HR 95% CI: 0.165 (better) to 1.345 (worse/no benefit)
    print()
    print("OS HR sensitivity (clinical efficacy uncertainty):")
    # Lower bound (HR=0.165, better survival)
    p = ModelParams(os_hr_scaling=0.35, constraint_general_pop=True)
    m = CEAModel(p)
    r = m.run()
    print(f"{'  OS HR lower bound (0.165, better)':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")
    # Upper bound (HR=1.345, no benefit)
    p = ModelParams(os_hr_scaling=2.86, constraint_general_pop=True)
    m = CEAModel(p)
    r = m.run()
    print(f"{'  OS HR upper bound (1.345, no benefit)':40s} {r['combo']['qaly']:>10.2f} {r['pembro']['qaly']:>10.2f} {r['combo']['qaly']-r['pembro']['qaly']:>10.2f} ${r['combo']['cost']-r['pembro']['cost']:>9,.0f} ${m.icer(r):>8,.0f}")


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