"""
Survival model selection: exponential, Weibull, log-normal, log-logistic, generalized gamma.
Fits OS, r_DM (=DMFS/OS), r_LR (=RFS/DMFS) curves; checks structure + plausibility.
"""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm
from model import ModelParams

p0 = ModelParams()
def pct(v): return np.array(v, dtype=float) / 100.0

t_os = np.array(p0.os_time_mo, float)
t_dmfs = np.array(p0.dmfs_time_mo, float)
t_rfs = np.array(p0.rfs_time_mo, float)

def interp_os_at(t_target, os_pct, os_times):
    """Interpolate OS percent to target times."""
    return np.interp(t_target, os_times, pct(os_pct))

# ── Data at observation times ──
# OS observed
os_c_obs = pct(p0.os_combo_pct)
os_p_obs = pct(p0.os_pembro_pct)
# r_DM = DMFS/OS at DMFS times
os_c_at_dmfs = interp_os_at(t_dmfs, p0.os_combo_pct, t_os)
os_p_at_dmfs = interp_os_at(t_dmfs, p0.os_pembro_pct, t_os)
rd_c_obs = np.clip(np.array(pct(p0.dmfs_combo_pct)) / os_c_at_dmfs, 0.01, 0.99)
rd_p_obs = np.clip(np.array(pct(p0.dmfs_pembro_pct)) / os_p_at_dmfs, 0.01, 0.99)
# r_LR = RFS/DMFS at RFS times
dmfs_c_at_rfs = np.interp(t_rfs, t_dmfs, pct(p0.dmfs_combo_pct))
dmfs_p_at_rfs = np.interp(t_rfs, t_dmfs, pct(p0.dmfs_pembro_pct))
rl_c_obs = np.clip(np.array(pct(p0.rfs_combo_pct)) / dmfs_c_at_rfs, 0.01, 0.99)
rl_p_obs = np.clip(np.array(pct(p0.rfs_pembro_pct)) / dmfs_p_at_rfs, 0.01, 0.99)

# ── Distributions ──
def s_expon(t, rate):
    return np.exp(-rate * t)

def s_weibull(t, scale, shape):
    return np.exp(-(t / scale) ** shape)

def s_lognorm(t, mu, sigma):
    return norm.cdf((mu - np.log(np.maximum(t, 1e-6))) / sigma)

def s_loglogistic(t, shape, scale):
    # S(t) = 1/(1 + (t/scale)^shape)
    return 1.0 / (1.0 + (t / scale) ** shape)

def s_gengamma(t, mu, sigma, q):
    # Generalized gamma via survival = 1 - F; use scipy still usable form:
    # Use the standard parameterization: T ~ GenGamma(shape=a, scale=b, c)
    from scipy.stats import gengamma
    # gengamma(a, c, scale=b).sf(t)
    # Fit using a = k, c, scale = b
    # We'll fit (a, c, scale) in log space
    return gengamma.sf(t, a=np.exp(mu), c=sigma, scale=np.exp(mu + q))  # placeholder

DISTS = ["exponential", "weibull", "lognormal", "loglogistic", "gengamma"]

def fit_dist(dist, t, s):
    """Fit distribution to (t, s), return (pred_fn, n_params, resid_sse)."""
    if dist == "exponential":
        def resid(p):
            return s_expon(t, p[0]) - s
        sol = least_squares(resid, x0=[0.01], bounds=([1e-6], [10]), xtol=1e-10, ftol=1e-10)
        pred = lambda tt: s_expon(np.array(tt, float), sol.x[0])
        return pred, 1, np.sum(sol.fun**2)
    if dist == "weibull":
        def resid(p):
            return s_weibull(t, p[0], p[1]) - s
        sol = least_squares(resid, x0=[t[-1]/(-np.log(max(s[-1],1e-6)))**(1/1.2), 1.2],
                            bounds=([10, 0.2], [5000, 5.0]), xtol=1e-10, ftol=1e-10)
        pred = lambda tt: s_weibull(np.array(tt, float), sol.x[0], sol.x[1])
        return pred, 2, np.sum(sol.fun**2)
    if dist == "lognormal":
        def resid(p):
            return s_lognorm(t, p[0], p[1]) - s
        sol = least_squares(resid, x0=[np.log(t[len(t)//2]), 0.8], bounds=([-5, 0.05], [15, 5.0]),
                            xtol=1e-12, ftol=1e-12, gtol=1e-12)
        pred = lambda tt: s_lognorm(np.array(tt, float), sol.x[0], sol.x[1])
        return pred, 2, np.sum(sol.fun**2)
    if dist == "loglogistic":
        def resid(p):
            return s_loglogistic(t, p[0], p[1]) - s
        sol = least_squares(resid, x0=[1.5, 50.0],
                            bounds=([0.05, 1], [15, 5000]), xtol=1e-10, ftol=1e-10)
        pred = lambda tt: s_loglogistic(np.array(tt, float), sol.x[0], sol.x[1])
        return pred, 2, np.sum(sol.fun**2)
    if dist == "gengamma":
        from scipy.stats import gengamma
        def resid(p):
            a, c, b = p
            return gengamma.sf(t, a=np.exp(a), c=c, scale=np.exp(b)) - s
        try:
            sol = least_squares(resid, x0=[0.5, 1.0, np.log(50)], bounds=([-2, 0.1, 1], [5, 5, 300]),
                                xtol=1e-10, ftol=1e-10, max_nfev=5000)
            a_, c_, b_ = sol.x
            pred = lambda tt: gengamma.sf(np.array(tt, float), a=np.exp(a_), c=c_, scale=np.exp(b_))
            return pred, 3, np.sum(sol.fun**2)
        except Exception:
            return None, 3, np.inf

# ── Fit all curves × dists ──
n_obs_os = len(t_os)
n_obs_dmfs = len(t_dmfs)
n_obs_rfs = len(t_rfs)

fits = {}  # key: (curve, dist) -> pred
print("=" * 110)
print("SURVIVAL MODEL SELECTION — INDIVIDUAL FIT (RMSE, AIC)")
print("=" * 110)

curves = [
    ("OS combo", t_os, os_c_obs, n_obs_os),
    ("OS pembro", t_os, os_p_obs, n_obs_os),
    ("r_DM combo", t_dmfs, rd_c_obs, n_obs_dmfs),
    ("r_DM pembro", t_dmfs, rd_p_obs, n_obs_dmfs),
    ("r_LR combo", t_rfs, rl_c_obs, n_obs_rfs),
    ("r_LR pembro", t_rfs, rl_p_obs, n_obs_rfs),
]

individual = {}
for cname, t_pts, s_obs, n in curves:
    print(f"\n--- {cname} (n={n} points) ---")
    print(f"{'Distribution':<16s} {'RMSE':>10s} {'AIC':>10s} {'params':>6s}")
    best = None
    for dist in DISTS:
        res = fit_dist(dist, t_pts, s_obs)
        if res is None:
            print(f"{dist:<16s} {'FAIL':>10s}")
            continue
        pred, k, sse = res
        # RMSE
        rmse = np.sqrt(sse / n)
        # AIC = n*ln(SSE/n) + 2k (exploratory; landmarks are serially correlated)
        aic = n * np.log(sse / n) + 2 * k
        individual[(cname, dist)] = pred
        print(f"{dist:<16s} {rmse:>10.4f} {aic:>10.2f} {k:>6d}")
        if best is None or aic < best[1]:
            best = (dist, aic)
    # Store best dist per curve
    if best:
        fits[cname] = best[0]
        print(f"  → Best AIC: {best[0]}")

# ── Composite model evaluation (choose dist per curve, check 40y structure) ──
print("\n" + "=" * 110)
print("COMPOSITE EVALUATION: pick dist for each curve, build 40y structure")
print("=" * 110)

mo = np.arange(1, 481, 1.0)  # 1..480 months
t_years = mo / 12.0

def build_model(dist_os_c, dist_os_p, dist_rd_c, dist_rd_p, dist_rl_c, dist_rl_p):
    """Return dict with os/dmfs/rfs arrays at mo, or None if structure fails."""
    try:
        pred_os_c = individual[("OS combo", dist_os_c)]
        pred_os_p = individual[("OS pembro", dist_os_p)]
        pred_rd_c = individual[("r_DM combo", dist_rd_c)]
        pred_rd_p = individual[("r_DM pembro", dist_rd_p)]
        pred_rl_c = individual[("r_LR combo", dist_rl_c)]
        pred_rl_p = individual[("r_LR pembro", dist_rl_p)]

        os_c = np.clip(pred_os_c(mo), 0, 1)
        os_p = np.clip(pred_os_p(mo), 0, 1)
        rd_c = np.clip(pred_rd_c(mo), 0, 1)
        rd_p = np.clip(pred_rd_p(mo), 0, 1)
        rl_c = np.clip(pred_rl_c(mo), 0, 1)
        rl_p = np.clip(pred_rl_p(mo), 0, 1)

        dmfs_c = os_c * rd_c
        dmfs_p = os_p * rd_p
        rfs_c = dmfs_c * rl_c
        rfs_p = dmfs_p * rl_p

        # Structural check at every month
        ok = True
        for a, rfs, dmfs, osv in [("combo", rfs_c, dmfs_c, os_c), ("pembro", rfs_p, dmfs_p, os_p)]:
            if np.any(rfs > dmfs + 1e-9) or np.any(dmfs > osv + 1e-9):
                ok = False
        # Monotonic
        for arr in [rfs_c, dmfs_c, os_c, rfs_p, dmfs_p, os_p]:
            if np.any(np.diff(arr) > 1e-9):
                ok = False

        return {"os_c": os_c, "os_p": os_p, "rfs_c": rfs_c, "rfs_p": rfs_p,
                "dmfs_c": dmfs_c, "dmfs_p": dmfs_p, "ordering": ok}
    except KeyError:
        return None

print(f"{'OS combo':<10s} {'OS pembro':<10s} {'rDM c':<7s} {'rDM p':<7s} {'rLR c':<7s} {'rLR p':<7s} {'Order 40y':>10s} {'OS40 combo':>10s} {'OS40 pembro':>10s}")

# All combos: each curve ∈ {weibull, lognormal} first (survey), then report
dist_options = ["exponential", "weibull", "lognormal", "loglogistic"]
results = []
for oc in dist_options:
    for op in dist_options:
        for rc in dist_options:
            for rp in dist_options:
                for lc in dist_options:
                    for lp in dist_options:
                        m = build_model(oc, op, rc, rp, lc, lp)
                        if m is None:
                            continue
                        os40_c = m["os_c"][-1]
                        os40_p = m["os_p"][-1]
                        results.append({
                            "oc": oc, "op": op, "rc": rc, "rp": rp, "lc": lc, "lp": lp,
                            "ordering": m["ordering"], "os40c": os40_c, "os40p": os40_p,
                            "rfs_c": m["rfs_c"], "os_c": m["os_c"], "dmfs_c": m["dmfs_c"],
                        })

# Print all combos ordering-valid
valid = [r for r in results if r["ordering"]]
print(f"\nOrdering-valid 40y combinations: {len(valid)} / {len(results)}")
print(f"\n{'#':<4s} {'OScombo':<10s} {'OSpembro':<10s} {'rDMc':<7s} {'rDMp':<7s} {'rLRc':<7s} {'rLRp':<7s} {'OS40c':>8s} {'OS40p':>8s} {'GP-ok':>6s}")

def gp_ok(os40c, os40p):
    # GP survival at year 40 (age 101) ~ 0.0; require OS40 < 0.10
    return os40c < 0.10 and os40p < 0.10

shown = 0
for i, r in enumerate(valid):
    gpok = gp_ok(r["os40c"], r["os40p"])
    mark = "✓" if (r["oc"] == "lognormal" and r["op"] == "lognormal" and r["rc"] == "lognormal" and
                   r["rp"] == "lognormal" and r["lc"] == "lognormal" and r["lp"] == "lognormal") else ""
    print(f"{i:<4d} {r['oc']:<10s} {r['op']:<10s} {r['rc']:<7s} {r['rp']:<7s} {r['lc']:<7s} {r['lp']:<7s} {r['os40c']:>8.3f} {r['os40p']:>8.3f} {'✓' if gpok else '✗':>6s} {mark}")
    shown += 1
    if shown >= 40:
        print("  ... (truncated)")
        break

# ── The log-normal baseline combo ──
print("\n" + "=" * 110)
print("LOG-NORMAL BASELINE VERIFICATION")
print("=" * 110)
ln_model = build_model("lognormal", "lognormal", "lognormal", "lognormal", "lognormal", "lognormal")
if ln_model:
    print(f"Ordering valid 40y: {'YES' if ln_model['ordering'] else 'NO'}")
    print(f"OS40 combo: {ln_model['os_c'][-1]:.4f}  OS40 pembro: {ln_model['os_p'][-1]:.4f}")
    # GP check: compare to life table at year 40
    print(f"GP survival year 40 (age 101): ~0.000")
    print(f"Clinical plausibility: {'OK - OS decays to GP floor' if ln_model['os_c'][-1] < 0.1 else 'SUSPECT'}")
else:
    print("Log-normal combo FAILED to build!")