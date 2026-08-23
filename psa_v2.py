"""
PSA v2 — exclusion-free probabilistic sensitivity analysis.

KEY CHANGE vs v1: draws are NEVER excluded based on the ICER value.
Negative ΔQALY, dominated cases, ICER > $1M are legitimate parts of the
uncertainty distribution and are all retained.

Validity is judged ONLY on trajectory-level structural criteria:
  - survival functions numerically computable
  - all survival probabilities within [0, 1]
  - monotonicity of RFS, DMFS, OS
  - ordering constraint RFS ≤ DMFS ≤ OS at every monthly cycle
  - state occupancies non-negative and sum to unity
  - no numerical solver failure

Iterations are resampled until N_VALID structurally valid iterations are
collected (fixed seed → fully reproducible).

Primary outputs (as recommended):
  ΔCost, ΔQALY, incremental NMB, CEAC  —  NOT an ICER 95% CI.
"""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
import json
from scipy.optimize import least_squares
from scipy.stats import norm, multivariate_normal
from model import ModelParams
from rerun_primary import CEAModelV2

N_VALID = 2000
SEED = 42

# ── Fit all 6 survival curves with Jacobian covariance ──
def fit_lognorm_cov(t, s):
    def resid(p):
        return norm.cdf((p[0] - np.log(t)) / p[1]) - s
    sol = least_squares(resid, x0=[np.log(t[len(t)//2]), 0.8],
                        bounds=([-5, 0.05], [15, 5.0]), xtol=1e-12, ftol=1e-12)
    J = sol.jac
    n, p_ = len(s), 2
    s2 = np.sum(sol.fun**2) / (n - p_)
    try:
        cov = np.linalg.inv(J.T @ J) * s2
    except Exception:
        cov = np.eye(2) * 0.01
    return float(sol.x[0]), float(sol.x[1]), cov

def fit_ratio_lognorm_cov(t, r):
    def resid(p):
        return norm.cdf((p[0] - np.log(t)) / p[1]) - r
    sol = least_squares(resid, x0=[np.log(t[len(t)//2]), 0.8],
                        bounds=([-5, 0.05], [15, 5.0]), xtol=1e-12, ftol=1e-12)
    J = sol.jac
    n, p_ = len(r), 2
    s2 = np.sum(sol.fun**2) / (n - p_)
    try:
        cov = np.linalg.inv(J.T @ J) * s2
    except Exception:
        cov = np.eye(2) * 0.01
    return float(sol.x[0]), float(sol.x[1]), cov

p0 = ModelParams()
def pct(v): return np.array(v, dtype=float) / 100.0

t_os = np.array(p0.os_time_mo)
t_dmfs = np.array(p0.dmfs_time_mo)
t_rfs = np.array(p0.rfs_time_mo)

os_mu_c, os_sig_c, os_cov_c = fit_lognorm_cov(t_os, pct(p0.os_combo_pct))
os_mu_p, os_sig_p, os_cov_p = fit_lognorm_cov(t_os, pct(p0.os_pembro_pct))

rd_c = np.clip(pct(p0.dmfs_combo_pct) / np.interp(t_dmfs, t_os, pct(p0.os_combo_pct)), 0.01, 0.99)
rd_p = np.clip(pct(p0.dmfs_pembro_pct) / np.interp(t_dmfs, t_os, pct(p0.os_pembro_pct)), 0.01, 0.99)
rd_mu_c, rd_sig_c, rd_cov_c = fit_ratio_lognorm_cov(t_dmfs, rd_c)
rd_mu_p, rd_sig_p, rd_cov_p = fit_ratio_lognorm_cov(t_dmfs, rd_p)

rl_c = np.clip(pct(p0.rfs_combo_pct) / np.interp(t_rfs, t_dmfs, pct(p0.dmfs_combo_pct)), 0.01, 0.99)
rl_p = np.clip(pct(p0.rfs_pembro_pct) / np.interp(t_rfs, t_dmfs, pct(p0.dmfs_pembro_pct)), 0.01, 0.99)
rl_mu_c, rl_sig_c, rl_cov_c = fit_ratio_lognorm_cov(t_rfs, rl_c)
rl_mu_p, rl_sig_p, rl_cov_p = fit_ratio_lognorm_cov(t_rfs, rl_p)

SURV = [
    ("os_mu_combo", "os_sigma_combo", os_mu_c, os_sig_c, os_cov_c),
    ("os_mu_pembro", "os_sigma_pembro", os_mu_p, os_sig_p, os_cov_p),
    ("rd_mu_combo", "rd_sigma_combo", rd_mu_c, rd_sig_c, rd_cov_c),
    ("rd_mu_pembro", "rd_sigma_pembro", rd_mu_p, rd_sig_p, rd_cov_p),
    ("rl_mu_combo", "rl_sigma_combo", rl_mu_c, rl_sig_c, rl_cov_c),
    ("rl_mu_pembro", "rl_sigma_pembro", rl_mu_p, rl_sig_p, rl_cov_p),
]

def sample_params(rng):
    pp = ModelParams(constraint_general_pop=True)
    for am, as_, mh, sh, cov in SURV:
        try:
            s = multivariate_normal.rvs([mh, sh], cov, random_state=rng)
        except Exception:
            s = [mh + rng.normal(0, 0.1), max(sh + rng.normal(0, 0.05), 0.05)]
        setattr(pp, am, float(s[0]))
        setattr(pp, as_, float(max(s[1], 0.05)))
    for attr, val in [("cost_keytruda_annual", 220_896), ("cost_intismeran", 200_000),
                      ("cost_lr_monthly", 3_000), ("cost_dm_monthly", 12_000)]:
        cv = 0.2
        setattr(pp, attr, float(rng.gamma(1/cv**2, val*cv**2)))
    for attr, val in [("util_rf", 0.83), ("util_lr", 0.64), ("util_dm", 0.55)]:
        se = 0.03
        alpha = val * (val*(1-val)/se**2 - 1)
        beta_p = (1-val) * (val*(1-val)/se**2 - 1)
        setattr(pp, attr, float(rng.beta(max(alpha, 1), max(beta_p, 1))))
    return pp

def structurally_valid(m):
    """Trajectory-level structural check. Returns (valid, reason)."""
    try:
        sur = m._survival()
    except Exception:
        return False, "survival computation failed"
    for arm in ["combo", "pembro"]:
        rfs = np.asarray(sur.get(f"rfs_{arm}"))
        dmfs = np.asarray(sur.get(f"dmfs_{arm}"))
        osv = np.asarray(sur.get(f"os_{arm}"))
        if len(rfs) == 0:
            return False, f"{arm}: missing survival"
        # in [0, 1]
        for name, arr in [("RFS", rfs), ("DMFS", dmfs), ("OS", osv)]:
            if np.any(arr < -1e-9) or np.any(arr > 1 + 1e-9):
                return False, f"{arm}: {name} outside [0,1]"
            if np.any(np.isnan(arr)) or np.any(np.isinf(arr)):
                return False, f"{arm}: {name} NaN/Inf"
        # monotonic non-increasing
        for name, arr in [("RFS", rfs), ("DMFS", dmfs), ("OS", osv)]:
            if np.any(np.diff(arr) > 1e-9):
                return False, f"{arm}: {name} not monotonic"
        # ordering RFS ≤ DMFS ≤ OS
        if np.any(rfs > dmfs + 1e-9):
            return False, f"{arm}: RFS > DMFS"
        if np.any(dmfs > osv + 1e-9):
            return False, f"{arm}: DMFS > OS"
        # state occupancy: RF, LR=RFS-DMFS, DM=DMFS-OS non-negative by ordering; sum = 1
        lr = dmfs - rfs
        dm = osv - dmfs
        states = osv + (1 - osv)  # rfs + lr + dm + dead = osv + (1-osv) = 1 by construction
        if np.any(np.abs(states - 1.0) > 1e-6):
            return False, f"{arm}: state sum != 1"
        if np.any(lr < -1e-9) or np.any(dm < -1e-9):
            return False, f"{arm}: negative state occupancy"
    return True, "valid"

rng = np.random.default_rng(SEED)
results = []          # (dqaly, dcost, dly)
rejects = {}          # reason -> count
n_draws = 0

while len(results) < N_VALID:
    n_draws += 1
    pp = sample_params(rng)
    try:
        m = CEAModelV2(pp)
        r = m.run()
    except Exception as e:
        rejects.setdefault("numerical failure", 0)
        rejects["numerical failure"] += 1
        continue
    ok, reason = structurally_valid(m)
    if not ok:
        rejects[reason] = rejects.get(reason, 0) + 1
        continue
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    sur = m._survival()
    t = m.t
    disc = np.exp(-pp.discount_rate * t)
    cl = pp.cycle_length
    dly = np.sum(sur["os_combo"] * disc * cl) - np.sum(sur["os_pembro"] * disc * cl)
    results.append((dq, dc, dly))

results = np.array(results)
dq = results[:, 0]
dc = results[:, 1]
dly = results[:, 2]

print(f"PSA v2 — exclusion-free, {len(results)} structurally valid iterations "
      f"(from {n_draws} draws, seed={SEED})")
print("Rejection reasons (structural only):")
for k, v in sorted(rejects.items(), key=lambda x: -x[1]):
    print(f"  {k}: {v}")
print(f"\nΔQALY: mean={dq.mean():.2f} median={np.median(dq):.2f} "
      f"95% CI [{np.percentile(dq,2.5):.2f}, {np.percentile(dq,97.5):.2f}]")
print(f"ΔLY:   mean={dly.mean():.2f} median={np.median(dly):.2f}")
print(f"ΔCost: mean=${dc.mean():,.0f} median=${np.median(dc):,.0f} "
      f"95% CI [${np.percentile(dc,2.5):,.0f}, ${np.percentile(dc,97.5):,.0f}]")

nmbs = {wtp: -dc + wtp * dq for wtp in [50_000, 100_000, 150_000]}
for wtp, nmb in nmbs.items():
    print(f"NMB @ ${wtp:,}: mean=${nmb.mean():,.0f}  P(NMB>0)={np.mean(nmb > 0):.4f}")

pce = {wtp: float(np.mean(-dc + wtp*dq > 0)) for wtp in [50_000, 100_000, 150_000]}
print(f"\nP(CE): $50K={pce[50_000]:.3f}  $100K={pce[100_000]:.3f}  $150K={pce[150_000]:.3f}")

# Fraction of cost-saving / dominated draws (legitimate uncertainty mass)
print(f"\nFraction ΔQALY < 0:     {np.mean(dq < 0):.3f}")
print(f"Fraction ΔCost < 0:      {np.mean(dc < 0):.3f}")
print(f"Fraction net cost-saving: {np.mean((dc < 0) & (dq > 0)):.3f}")

# Save full draws for CE plane / convergence
np.savez("/tmp/figures/psa_v2_draws.npz",
         dq=dq, dc=dc, dly=dly, seed=SEED, n_valid=len(results), n_draws=n_draws,
         pce50=pce[50_000], pce100=pce[100_000], pce150=pce[150_000])
print("Saved /tmp/figures/psa_v2_draws.npz")

summary = {
    "n_valid": len(results), "n_draws": n_draws,
    "dq_mean": dq.mean(), "dq_median": float(np.median(dq)),
    "dc_mean": dc.mean(), "dc_median": float(np.median(dc)),
    "pce50": pce[50_000], "pce100": pce[100_000], "pce150": pce[150_000],
    "nmb100_mean": nmbs[100_000].mean(),
    "frac_saving": float(np.mean((dc < 0) & (dq > 0))),
}
print("SUMMARY_JSON=" + json.dumps(summary))