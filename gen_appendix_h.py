"""Generate all Appendix H data: DSA table, PSA convergence/scenarios, CE plane, structural scenarios."""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.stats import norm, multivariate_normal, beta as beta_dist
from model import ModelParams
from model import discount_factors  # noqa: E402
from rerun_primary import CEAModelV2, general_pop_surv, LT_hazard_monthly

plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "figure.dpi": 300, "savefig.dpi": 300})

# ──────────────────────────────────
# H.1 DSA Table
# ──────────────────────────────────
print("=" * 80)
print("H.1 DSA — ONE-WAY SENSITIVITY ANALYSIS")
print("=" * 80)

def run_model(**kwargs):
    p = ModelParams(constraint_general_pop=True)
    for k, v in kwargs.items():
        setattr(p, k, v)
    m = CEAModelV2(p)
    r = m.run()
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    icer = dc / dq if dq > 0 else float('inf')
    nmb_100k = dc - 100_000 * dq
    nmb_150k = dc - 150_000 * dq
    return dc, dq, icer, nmb_100k, nmb_150k

base_dc, base_dq, base_icer, base_nmb100k, base_nmb150k = run_model()

dsa_params = {
    "Utility RF": ("util_rf", 0.83, 0.80, 0.86),
    "Utility LR": ("util_lr", 0.64, 0.61, 0.67),
    "Utility DM": ("util_dm", 0.55, 0.52, 0.58),
    "Keytruda cost": ("cost_keytruda_annual", 220_896, 176_717, 265_075),
    "Intismeran cost": ("cost_intismeran", 200_000, 160_000, 240_000),
    "LR monthly cost": ("cost_lr_monthly", 3_000, 2_400, 3_600),
    "DM monthly cost": ("cost_dm_monthly", 12_000, 9_600, 14_400),
    "Discount rate": ("discount_rate", 0.03, 0.0, 0.05),
    "OS mu combo": ("os_mu_combo", 4.428, 4.228, 4.628),
    "OS mu pembro": ("os_mu_pembro", 3.725, 3.525, 3.925),
}

print(f"{'Parameter':<22s} {'Base':>10s} {'Low':>10s} {'High':>10s} {'Low ICER':>10s} {'High ICER':>10s} {'ΔNMB@100K':>10s}")
print("-" * 82)
for label, (attr, base_val, lo, hi) in sorted(dsa_params.items()):
    _, _, icer_lo, nmb_lo, _ = run_model(**{attr: lo})
    _, _, icer_hi, nmb_hi, _ = run_model(**{attr: hi})
    dnmb_100k = abs(nmb_hi - nmb_lo)
    print(f"{label:<22s} {base_val:>10} {lo:>10} {hi:>10} {icer_lo:>10,.0f} {icer_hi:>10,.0f} {dnmb_100k:>10,.0f}")

# ──────────────────────────────────
# H.2 PSA — convergence test
# ──────────────────────────────────
print("\n" + "=" * 80)
print("H.2 PSA CONVERGENCE")
print("=" * 80)

# Fit survival params with cov
p0 = ModelParams()
def pct(v): return np.array(v, dtype=float) / 100.0

def fit_cov(t, s):
    def resid(p):
        return norm.cdf((p[0] - np.log(t)) / p[1]) - s
    sol = least_squares(resid, x0=[np.log(t[len(t)//2]), 0.8], bounds=([-5, 0.05], [15, 5.0]), xtol=1e-12, ftol=1e-12)
    J = sol.jac
    s2 = np.sum(sol.fun**2) / (len(s) - 2)
    try:
        cov = np.linalg.inv(J.T @ J) * s2
    except:
        cov = np.eye(2) * 0.01
    return sol.x[0], sol.x[1], cov

t_os = np.array(p0.os_time_mo)
t_dmfs = np.array(p0.dmfs_time_mo)
t_rfs = np.array(p0.rfs_time_mo)

os_mu_c, os_sig_c, os_cov_c = fit_cov(t_os, pct(p0.os_combo_pct))
os_mu_p, os_sig_p, os_cov_p = fit_cov(t_os, pct(p0.os_pembro_pct))

os_c_dmfs = np.interp(t_dmfs, t_os, pct(p0.os_combo_pct))
os_p_dmfs = np.interp(t_dmfs, t_os, pct(p0.os_pembro_pct))
rd_c = np.clip(pct(p0.dmfs_combo_pct) / os_c_dmfs, 0.01, 0.99)
rd_p = np.clip(pct(p0.dmfs_pembro_pct) / os_p_dmfs, 0.01, 0.99)
rd_mu_c, rd_sig_c, rd_cov_c = fit_cov(t_dmfs, rd_c)
rd_mu_p, rd_sig_p, rd_cov_p = fit_cov(t_dmfs, rd_p)

dmfs_c_rfs = np.interp(t_rfs, t_dmfs, pct(p0.dmfs_combo_pct))
dmfs_p_rfs = np.interp(t_rfs, t_dmfs, pct(p0.dmfs_pembro_pct))
rl_c = np.clip(pct(p0.rfs_combo_pct) / dmfs_c_rfs, 0.01, 0.99)
rl_p = np.clip(pct(p0.rfs_pembro_pct) / dmfs_p_rfs, 0.01, 0.99)
rl_mu_c, rl_sig_c, rl_cov_c = fit_cov(t_rfs, rl_c)
rl_mu_p, rl_sig_p, rl_cov_p = fit_cov(t_rfs, rl_p)

surv_params = [
    ("os_mu_combo", "os_sigma_combo", os_mu_c, os_sig_c, os_cov_c),
    ("os_mu_pembro", "os_sigma_pembro", os_mu_p, os_sig_p, os_cov_p),
    ("rd_mu_combo", "rd_sigma_combo", rd_mu_c, rd_sig_c, rd_cov_c),
    ("rd_mu_pembro", "rd_sigma_pembro", rd_mu_p, rd_sig_p, rd_cov_p),
    ("rl_mu_combo", "rl_sigma_combo", rl_mu_c, rl_sig_c, rl_cov_c),
    ("rl_mu_pembro", "rl_sigma_pembro", rl_mu_p, rl_sig_p, rl_cov_p),
]

n_iters_list = [100, 500, 1000, 2500]
rng = np.random.default_rng(42)
all_q = []
all_c = []
all_ce = []

for n_iter in n_iters_list:
    dqs, dcs = [], []
    for it in range(n_iter):
        pp = ModelParams(constraint_general_pop=True)
        for attr_mu, attr_sig, mu_h, sig_h, cov in surv_params:
            try:
                s = np.random.multivariate_normal([mu_h, sig_h], cov)
            except:
                s = [mu_h + rng.normal(0, 0.1), max(sig_h + rng.normal(0, 0.05), 0.05)]
            setattr(pp, attr_mu, s[0])
            setattr(pp, attr_sig, max(s[1], 0.05))
        for attr, val in [("cost_keytruda_annual", 220_896), ("cost_intismeran", 200_000),
                          ("cost_lr_monthly", 3_000), ("cost_dm_monthly", 12_000)]:
            cv = 0.2
            setattr(pp, attr, rng.gamma(1/cv**2, val*cv**2))
        for attr, val in [("util_rf", 0.83), ("util_lr", 0.64), ("util_dm", 0.55)]:
            se = 0.03
            alpha = val * (val*(1-val)/se**2 - 1)
            beta_p = (1-val) * (val*(1-val)/se**2 - 1)
            setattr(pp, attr, beta_dist.rvs(max(alpha,1), max(beta_p,1), random_state=rng))
        try:
            m = CEAModelV2(pp)
            r = m.run()
            dc = r["combo"]["cost"] - r["pembro"]["cost"]
            dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
            ic = dc / dq if dq > 0 else float('inf')
            if 0 < ic < 1e6:
                dqs.append(dq)
                dcs.append(dc)
        except:
            pass

    all_q.extend(dqs)
    all_c.extend(dcs)
    nmb = np.array(dcs) - 100_000 * np.array(dqs)
    pce = np.mean(nmb < 0)
    print(f"N={n_iter:>5d}: ΔQALY={np.mean(dqs):.2f} ΔCost=${np.mean(dcs):>,.0f} NMB=${np.mean(nmb):>,.0f} P(CE@100K)={pce:.3f}")

# ──────────────────────────────────
# Figure S6: PSA convergence
# ──────────────────────────────────
# Use cumulative means
cum_q = np.cumsum(all_q) / np.arange(1, len(all_q)+1)
cum_c = np.cumsum(all_c) / np.arange(1, len(all_c)+1)
cum_nmb = np.array(all_c) - 100_000 * np.array(all_q)
cum_nmb_mean = np.cumsum(cum_nmb) / np.arange(1, len(cum_nmb)+1)

fig, axes = plt.subplots(2, 2, figsize=(10, 7))
for ax, data, label, unit in [
    (axes[0,0], cum_q, "Mean ΔQALY", "QALYs"),
    (axes[0,1], cum_c, "Mean ΔCost", "USD"),
    (axes[1,0], cum_nmb_mean, "Mean NMB ($100K WTP)", "USD"),
    (axes[1,1], np.array([np.mean(np.array(all_c[:i+1]) - 100_000*np.array(all_q[:i+1]) < 0) for i in range(len(all_q))]), "P(CE) @ $100K", "Probability"),
]:
    ax.plot(range(1, len(data)+1), data, "-", color="#2166ac", lw=1)
    # 1,000 iteration mark
    ax.axvline(1000, color="gray", ls=":", alpha=0.5)
    ax.set_xlabel("PSA Iterations")
    ax.set_ylabel(label)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("/tmp/figures/fig_s6_psa_convergence.png", dpi=300)
plt.savefig("/tmp/figures/fig_s6_psa_convergence.pdf")
print("\nSaved Figure S6: PSA convergence")

# ──────────────────────────────────
# Figure S7: CE plane
# ──────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(all_q, all_c, s=5, alpha=0.3, color="#2166ac", label="PSA iterations")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
ax.axvline(0, color="gray", ls=":", alpha=0.5)
# WTP lines
for wtp, color, ls in [(50_000, "#e41a1c", "-"), (100_000, "#4daf4a", "-"), (150_000, "#984ea3", "-")]:
    x = np.linspace(0, 10, 100)
    ax.plot(x, wtp * x, color=color, ls=ls, lw=1, label=f"WTP=${wtp:,}/QALY")
ax.set_xlabel("Incremental QALYs")
ax.set_ylabel("Incremental Cost (USD)")
ax.legend(fontsize=8, loc="upper left")
ax.grid(True, alpha=0.3)
# Base case point
ax.scatter([base_dq], [base_dc], s=80, color="black", marker="D", zorder=10, label="Base case")
ax.set_xlim(0, 8)
ax.set_ylim(0, 600_000)
plt.tight_layout()
plt.savefig("/tmp/figures/fig_s7_ce_plane.png", dpi=300)
plt.savefig("/tmp/figures/fig_s7_ce_plane.pdf")
print("Saved Figure S7: CE plane")

# ──────────────────────────────────
# H.3 Structural Scenarios
# ──────────────────────────────────
print("\n" + "=" * 80)
print("H.3 STRUCTURAL SCENARIOS")
print("=" * 80)

def run_scenario(**kwargs):
    p = ModelParams(constraint_general_pop=True)
    for k, v in kwargs.items():
        setattr(p, k, v)
    m = CEAModelV2(p)
    r = m.run()
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    # LY
    sur = m._survival()
    t = m.t
    disc = discount_factors(p.discount_rate, t)
    cl = p.cycle_length
    ly_c = np.sum(sur["os_combo"] * disc * cl)
    ly_p = np.sum(sur["os_pembro"] * disc * cl)
    icer = dc / dq if dq > 0 else float('inf')
    nmb_100k = -dc + 100_000 * dq
    nmb_150k = -dc + 150_000 * dq
    return ly_c - ly_p, dq, dc, icer, nmb_100k, nmb_150k

# GP floor variants: need to hack _apply_gp_hazard_floor
class CEAModelFlex(CEAModelV2):
    def _apply_gp_hazard_floor(self, os_val, mo):
        n = len(os_val)
        h_model = np.zeros(n)
        for i in range(n - 1):
            h_model[i] = -np.log(max(os_val[i+1] / max(os_val[i], 1e-12), 1e-12))
        gp_start = getattr(self.p, 'gp_floor_start_month', 60)
        mask = mo >= gp_start
        gp_h = np.interp(np.arange(n), np.arange(480), LT_hazard_monthly)
        h_final = np.where(mask, np.maximum(h_model, gp_h), h_model)
        h_final[-1] = h_final[-2]
        cum_h = np.cumsum(h_final)
        new_os = np.exp(-cum_h)
        for i in range(1, n):
            new_os[i] = min(new_os[i], new_os[i-1])
        return new_os

# Waning V2 (hazard convergence)
class CEAModelWaningV2(CEAModelV2):
    def _os_fn(self, mo, arm):
        prefix = "combo" if arm == "combo" else "pembro"
        mu = getattr(self.p, f"os_mu_{prefix}")
        sigma = getattr(self.p, f"os_sigma_{prefix}")
        os_val = norm.cdf((mu - np.log(np.maximum(mo, 1e-6))) / sigma)
        if self.p.constraint_general_pop:
            os_val = self._apply_gp_hazard_floor(os_val, mo)
        if arm == "combo":
            p_prefix = "pembro"
            p_mu = getattr(self.p, f"os_mu_{p_prefix}")
            p_sigma = getattr(self.p, f"os_sigma_{p_prefix}")
            p_os = norm.cdf((p_mu - np.log(np.maximum(mo, 1e-6))) / p_sigma)
            if self.p.constraint_general_pop:
                p_os = self._apply_gp_hazard_floor(p_os, mo)
            # Hazard convergence: wane combo hazard to pembro hazard (5-10y)
            n = len(mo)
            h_c = np.zeros(n)
            h_p = np.zeros(n)
            for i in range(n - 1):
                h_c[i] = -np.log(max(os_val[i+1] / max(os_val[i], 1e-12), 1e-12))
                h_p[i] = -np.log(max(p_os[i+1] / max(p_os[i], 1e-12), 1e-12))
            wane = np.clip((mo - 60) / (120 - 60), 0, 1)
            h_final = (1 - wane) * h_c + wane * h_p
            h_final[-1] = h_final[-2]
            cum_h = np.cumsum(h_final)
            os_val = np.exp(-cum_h)
            for i in range(1, n):
                os_val[i] = min(os_val[i], os_val[i-1])
        return os_val

scenarios = [
    ("Base case (GP-constrained, log-normal)", lambda: run_scenario(constraint_general_pop=True)),
    ("Weibull OS", lambda: run_scenario(constraint_general_pop=True, os_distribution="weibull")),
    ("Unconstrained (no GP)", lambda: run_scenario(constraint_general_pop=False)),
    ("GP floor at 48 mo (4 yr)", lambda: (
        lambda: (lambda p: (setattr(p, 'gp_floor_start_month', 48),
         CEAModelFlex(p).run()))(ModelParams(constraint_general_pop=True)))[0]()),
    ("GP floor at 60 mo (5 yr)", lambda: run_scenario(constraint_general_pop=True)),
    ("GP floor at 72 mo (6 yr)", lambda: (
        lambda: (lambda p: (setattr(p, 'gp_floor_start_month', 72),
         CEAModelFlex(p).run()))(ModelParams(constraint_general_pop=True)))[0]()),
    ("10-year horizon", lambda: run_scenario(constraint_general_pop=True, time_horizon_years=10)),
    ("20-year horizon", lambda: run_scenario(constraint_general_pop=True, time_horizon_years=20)),
    ("40-year horizon (base)", lambda: run_scenario(constraint_general_pop=True, time_horizon_years=40)),
    ("Price $100K", lambda: run_scenario(constraint_general_pop=True, cost_intismeran=100_000)),
    ("Price $200K (base)", lambda: run_scenario(constraint_general_pop=True, cost_intismeran=200_000)),
    ("Price $300K", lambda: run_scenario(constraint_general_pop=True, cost_intismeran=300_000)),
    ("Price $500K", lambda: run_scenario(constraint_general_pop=True, cost_intismeran=500_000)),
    ("Recurrence-effect waning (V1)", lambda: run_scenario(constraint_general_pop=True, treatment_waning=True)),
    ("Hazard-convergence waning (V2)", lambda: (
        lambda: (lambda p, m: (m.run()))(ModelParams(constraint_general_pop=True),
         CEAModelWaningV2(ModelParams(constraint_general_pop=True))))[0]()),
    ("No direct OS benefit", lambda: run_scenario(constraint_general_pop=True, no_direct_os_benefit=True)),
]

print(f"{'Scenario':<42s} {'ΔLY':>6s} {'ΔQALY':>7s} {'ΔCost':>10s} {'ICER':>10s} {'NMB@100K':>10s} {'NMB@150K':>10s}")
print("-" * 95)
for label, fn in scenarios:
    try:
        dy, dq, dc, icer, nmb100, nmb150 = fn()
        icer_str = f"${icer:>,.0f}" if icer < 1e6 else "Dom"
        print(f"{label:<42s} {dy:>6.2f} {dq:>7.2f} ${dc:>8,.0f} {icer_str:>10s} ${nmb100:>9,.0f} ${nmb150:>9,.0f}")
    except Exception as e:
        print(f"{label:<42s} ERROR: {e}")

# Save CE plane data
np.savez("/tmp/figures/ce_plane_data.npz", dq=all_q, dc=all_c, base_dq=base_dq, base_dc=base_dc)
print("\nSaved CE plane data")
