"""
Proper PSA with full survival uncertainty propagation.
- Fits log-normal curves with Jacobian-based covariance for (mu, sigma)
- Jointly samples (mu, sigma) for all 6 survival curves (OS×2, r_DM×2, r_LR×2)
- GP-constrained model (CEAModelV2)
- Reports ΔQALY, ΔCost, ICER, P(CE)
"""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm, multivariate_normal
from model import ModelParams
from rerun_primary import CEAModelV2

def fit_lognorm_cov(time_points, surv_rates):
    """Fit log-normal, return (mu, sigma, cov_matrix)."""
    t = np.array(time_points, float)
    s = np.array(surv_rates, float)
    def resid(p):
        mu, sigma = p
        return norm.cdf((mu - np.log(t)) / sigma) - s
    sol = least_squares(resid, x0=[np.log(t[len(t)//2]), 0.8],
                        bounds=([-5, 0.05], [15, 5.0]),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    mu, sigma = sol.x
    # Covariance from Jacobian: cov = s² * (JᵀJ)⁻¹
    J = sol.jac
    n = len(s)
    p = 2
    s2 = np.sum(sol.fun**2) / (n - p)  # residual variance
    try:
        cov = np.linalg.inv(J.T @ J) * s2
    except:
        cov = np.eye(2) * 0.01
    return mu, sigma, cov

def fit_ratio_lognorm_cov(time_points, ratio_values):
    """Fit log-normal ratio, return (mu, sigma, cov_matrix)."""
    t = np.array(time_points, float)
    r = np.array(ratio_values, float)
    def resid(p):
        mu, sigma = p
        return norm.cdf((mu - np.log(t)) / sigma) - r
    sol = least_squares(resid, x0=[np.log(t[len(t)//2]), 0.8],
                        bounds=([-5, 0.05], [15, 5.0]),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    mu, sigma = sol.x
    J = sol.jac
    n = len(r)
    p = 2
    s2 = np.sum(sol.fun**2) / (n - p)
    try:
        cov = np.linalg.inv(J.T @ J) * s2
    except:
        cov = np.eye(2) * 0.01
    return mu, sigma, cov

# ── Fit all 6 curves with covariance ──
p0 = ModelParams()
def pct(v): return np.array(v, dtype=float) / 100.0

# OS
os_mu_c, os_sigma_c, os_cov_c = fit_lognorm_cov(p0.os_time_mo, pct(p0.os_combo_pct))
os_mu_p, os_sigma_p, os_cov_p = fit_lognorm_cov(p0.os_time_mo, pct(p0.os_pembro_pct))

# r_DM = DMFS/OS
# Compute r_DM at each time point
t_os = np.array(p0.os_time_mo)
t_dmfs = np.array(p0.dmfs_time_mo)
# Interpolate OS to DMFS time points for ratio
os_c_dmfs = np.interp(t_dmfs, t_os, pct(p0.os_combo_pct))
os_p_dmfs = np.interp(t_dmfs, t_os, pct(p0.os_pembro_pct))
rd_c = np.array(pct(p0.dmfs_combo_pct)) / os_c_dmfs
rd_p = np.array(pct(p0.dmfs_pembro_pct)) / os_p_dmfs
rd_c = np.clip(rd_c, 0.01, 0.99)
rd_p = np.clip(rd_p, 0.01, 0.99)
rd_mu_c, rd_sigma_c, rd_cov_c = fit_ratio_lognorm_cov(p0.dmfs_time_mo, rd_c)
rd_mu_p, rd_sigma_p, rd_cov_p = fit_ratio_lognorm_cov(p0.dmfs_time_mo, rd_p)

# r_LR = RFS/DMFS
rfs_c = np.array(pct(p0.rfs_combo_pct))
rfs_p = np.array(pct(p0.rfs_pembro_pct))
t_rfs = np.array(p0.rfs_time_mo)
dmfs_c_rfs = np.interp(t_rfs, t_dmfs, pct(p0.dmfs_combo_pct))
dmfs_p_rfs = np.interp(t_rfs, t_dmfs, pct(p0.dmfs_pembro_pct))
rl_c = rfs_c / dmfs_c_rfs
rl_p = rfs_p / dmfs_p_rfs
rl_c = np.clip(rl_c, 0.01, 0.99)
rl_p = np.clip(rl_p, 0.01, 0.99)
rl_mu_c, rl_sigma_c, rl_cov_c = fit_ratio_lognorm_cov(p0.rfs_time_mo, rl_c)
rl_mu_p, rl_sigma_p, rl_cov_p = fit_ratio_lognorm_cov(p0.rfs_time_mo, rl_p)

print("=== Fitted parameters with covariance ===")
for label, mu, sigma, cov in [
    ("OS combo", os_mu_c, os_sigma_c, os_cov_c),
    ("OS pembro", os_mu_p, os_sigma_p, os_cov_p),
    ("r_DM combo", rd_mu_c, rd_sigma_c, rd_cov_c),
    ("r_DM pembro", rd_mu_p, rd_sigma_p, rd_cov_p),
    ("r_LR combo", rl_mu_c, rl_sigma_c, rl_cov_c),
    ("r_LR pembro", rl_mu_p, rl_sigma_p, rl_cov_p),
]:
    se_mu = np.sqrt(cov[0,0])
    se_sigma = np.sqrt(cov[1,1])
    corr = cov[0,1] / np.sqrt(cov[0,0]*cov[1,1])
    print(f"  {label:<20s} mu={mu:>6.3f}(SE={se_mu:.3f}) sigma={sigma:>6.3f}(SE={se_sigma:.3f}) r={corr:.3f}")

# ── PSA: 1,000 iterations with joint (mu, sigma) sampling ──
print("\n=== PSA (1,000 iterations, joint survival sampling) ===")
rng = np.random.default_rng(42)
icers = []
costs_diff = []
qalys_diff = []
n_ok = 0
os_10y = []
os_20y = []

for it in range(1000):
    pp = ModelParams(constraint_general_pop=True)

    # Sample survival params jointly (mu, sigma)
    for attr_mu, attr_sigma, mu_hat, sigma_hat, cov in [
        ("os_mu_combo", "os_sigma_combo", os_mu_c, os_sigma_c, os_cov_c),
        ("os_mu_pembro", "os_sigma_pembro", os_mu_p, os_sigma_p, os_cov_p),
        ("rd_mu_combo", "rd_sigma_combo", rd_mu_c, rd_sigma_c, rd_cov_c),
        ("rd_mu_pembro", "rd_sigma_pembro", rd_mu_p, rd_sigma_p, rd_cov_p),
        ("rl_mu_combo", "rl_sigma_combo", rl_mu_c, rl_sigma_c, rl_cov_c),
        ("rl_mu_pembro", "rl_sigma_pembro", rl_mu_p, rl_sigma_p, rl_cov_p),
    ]:
        try:
            sampled = multivariate_normal.rvs([mu_hat, sigma_hat], cov, random_state=rng)
        except:
            sampled = [mu_hat + rng.normal(0, 0.1), max(sigma_hat + rng.normal(0, 0.05), 0.05)]
        setattr(pp, attr_mu, sampled[0])
        setattr(pp, attr_sigma, max(sampled[1], 0.05))

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

    m = CEAModelV2(pp)
    try:
        r = m.run()
        ic = m.icer(r)
        if 0 < ic < 1e6:
            icers.append(ic)
            costs_diff.append(r["combo"]["cost"] - r["pembro"]["cost"])
            qalys_diff.append(r["combo"]["qaly"] - r["pembro"]["qaly"])
            n_ok += 1
            # Track OS at 10y and 20y
            idx10 = int(10 / m.p.cycle_length)
            idx20 = int(20 / m.p.cycle_length)
            sur = m._survival()
            os_10y.append((sur["os_combo"][idx10], sur["os_pembro"][idx10]))
            os_20y.append((sur["os_combo"][idx20], sur["os_pembro"][idx20]))
    except:
        continue

ic_arr = np.array(icers)
cd_arr = np.array(costs_diff)
qd_arr = np.array(qalys_diff)
os10_arr = np.array(os_10y)
os20_arr = np.array(os_20y)

print(f"Valid iterations: {n_ok}/1000")
print(f"\nΔQALY: mean={np.mean(qd_arr):.2f} median={np.median(qd_arr):.2f} "
      f"95% CI [{np.percentile(qd_arr,2.5):.2f}–{np.percentile(qd_arr,97.5):.2f}]")
print(f"ΔCost: mean=${np.mean(cd_arr):,.0f} median=${np.median(cd_arr):,.0f} "
      f"95% CI [${np.percentile(cd_arr,2.5):,.0f}–${np.percentile(cd_arr,97.5):,.0f}]")
print(f"ICER:  mean=${np.mean(ic_arr):,.0f} median=${np.median(ic_arr):,.0f} "
      f"95% CI [${np.percentile(ic_arr,2.5):,.0f}–${np.percentile(ic_arr,97.5):,.0f}]")
for thresh, label in [(50_000, "$50K"), (100_000, "$100K"), (150_000, "$150K")]:
    print(f"  P(CE at {label}): {np.mean(ic_arr <= thresh):.4f}")

# Survival uncertainty propagation check
print(f"\nCombo OS at 10y: mean={np.mean(os10_arr[:,0]):.3f} "
      f"95% CI [{np.percentile(os10_arr[:,0],2.5):.3f}–{np.percentile(os10_arr[:,0],97.5):.3f}]")
print(f"Pembro OS at 10y: mean={np.mean(os10_arr[:,1]):.3f} "
      f"95% CI [{np.percentile(os10_arr[:,1],2.5):.3f}–{np.percentile(os10_arr[:,1],97.5):.3f}]")
print(f"Combo OS at 20y: mean={np.mean(os20_arr[:,0]):.3f} "
      f"95% CI [{np.percentile(os20_arr[:,0],2.5):.3f}–{np.percentile(os20_arr[:,0],97.5):.3f}]")
print(f"Pembro OS at 20y: mean={np.mean(os20_arr[:,1]):.3f} "
      f"95% CI [{np.percentile(os20_arr[:,1],2.5):.3f}–{np.percentile(os20_arr[:,1],97.5):.3f}]")