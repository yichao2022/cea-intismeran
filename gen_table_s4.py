"""Generate Table S4 — Joint model selection: 15,625 combinations, top 20 by total AIC + structural validity."""
import sys
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm, gengamma
from model import ModelParams

p = ModelParams()
def pct(v): return np.array(v, dtype=float) / 100.0

t_os = np.array(p.os_time_mo, float)
t_dmfs = np.array(p.dmfs_time_mo, float)
t_rfs = np.array(p.rfs_time_mo, float)
os_c, os_p = pct(p.os_combo_pct), pct(p.os_pembro_pct)
dmfs_c, dmfs_p = pct(p.dmfs_combo_pct), pct(p.dmfs_pembro_pct)
rfs_c, rfs_p = pct(p.rfs_combo_pct), pct(p.rfs_pembro_pct)

# Interpolate OS at DMFS times for r_DM
def interp(t_target, t_obs, s_obs):
    return np.interp(t_target, t_obs, s_obs)
os_c_at_dmfs = interp(t_dmfs, t_os, os_c)
os_p_at_dmfs = interp(t_dmfs, t_os, os_p)
rd_c = np.clip(dmfs_c / os_c_at_dmfs, 0.01, 0.99)
rd_p = np.clip(dmfs_p / os_p_at_dmfs, 0.01, 0.99)
dmfs_c_at_rfs = np.interp(t_rfs, t_dmfs, dmfs_c)
dmfs_p_at_rfs = np.interp(t_rfs, t_dmfs, dmfs_p)
rl_c = np.clip(rfs_c / dmfs_c_at_rfs, 0.01, 0.99)
rl_p = np.clip(rfs_p / dmfs_p_at_rfs, 0.01, 0.99)

CURVES = {
    "OS_c": (t_os, os_c, 5), "OS_p": (t_os, os_p, 5),
    "rDM_c": (t_dmfs, rd_c, 4), "rDM_p": (t_dmfs, rd_p, 4),
    "rLR_c": (t_rfs, rl_c, 5), "rLR_p": (t_rfs, rl_p, 5),
}
DISTS = ["exponential", "weibull", "lognormal", "loglogistic", "gengamma"]

def s_expon(t, r): return np.exp(-r * t)
def s_weibull(t, a, b): return np.exp(-(t / a) ** b)
def s_lognorm(t, mu, sig): return norm.cdf((mu - np.log(np.maximum(t, 1e-6))) / sig)
def s_loglog(t, a, b): return 1.0 / (1.0 + (t / b) ** a)

def fit_dist(dist, t, s):
    n = len(t)
    if dist == "exponential":
        sol = least_squares(lambda pp: s_expon(t, pp[0]) - s, x0=[0.01], bounds=([1e-6], [10]))
        pred = lambda tt: s_expon(np.array(tt, float), sol.x[0])
        k = 1
    elif dist == "weibull":
        scale0 = t[-1] / (-np.log(max(s[-1], 1e-6))) ** (1 / 1.2)
        sol = least_squares(lambda pp: s_weibull(t, pp[0], pp[1]) - s, x0=[scale0, 1.2], bounds=([10, 0.2], [5000, 5]))
        pred = lambda tt: s_weibull(np.array(tt, float), sol.x[0], sol.x[1])
        k = 2
    elif dist == "lognormal":
        sol = least_squares(lambda pp: s_lognorm(t, pp[0], pp[1]) - s, x0=[np.log(t[len(t) // 2]), 0.8], bounds=([-5, 0.05], [15, 5]))
        pred = lambda tt: s_lognorm(np.array(tt, float), sol.x[0], sol.x[1])
        k = 2
    elif dist == "loglogistic":
        sol = least_squares(lambda pp: s_loglog(t, pp[0], pp[1]) - s, x0=[1.5, 50.0], bounds=([0.05, 1], [15, 5000]))
        pred = lambda tt: s_loglog(np.array(tt, float), sol.x[0], sol.x[1])
        k = 2
    elif dist == "gengamma":
        sol = least_squares(lambda pp: gengamma.sf(t, a=np.exp(pp[0]), c=pp[1], scale=np.exp(pp[2])) - s,
                            x0=[0.5, 1.0, np.log(50)], bounds=([-2, 0.1, 1], [5, 5, 300]),
                            max_nfev=5000)
        pred = lambda tt: gengamma.sf(np.array(tt, float), a=np.exp(sol.x[0]), c=sol.x[1], scale=np.exp(sol.x[2]))
        k = 3
    sse = np.sum((pred(t) - s) ** 2)
    aic = n * np.log(sse / n) + 2 * k
    rmse = np.sqrt(sse / n)
    return pred, aic, rmse

# Fit all curves x dists
fits = {}
for cname, (t, s, n) in CURVES.items():
    for dist in DISTS:
        pred, aic, rmse = fit_dist(dist, t, s)
        fits[(cname, dist)] = (pred, aic, rmse)

# Joint model: choose dist per curve, compute total AIC + 40y structural validity
mo = np.arange(1, 481, 1.0)
t_yrs = mo / 12.0

def build_and_score(combo):
    oc, op, dc, dp, lc, lp = combo
    preds = {}
    for key, dist in zip(["OS_c", "OS_p", "rDM_c", "rDM_p", "rLR_c", "rLR_p"], combo):
        preds[key] = fits[(key, dist)][0]
    os_c40 = np.clip(preds["OS_c"](mo), 0, 1)
    os_p40 = np.clip(preds["OS_p"](mo), 0, 1)
    rd_c40 = np.clip(preds["rDM_c"](mo), 0, 1)
    rd_p40 = np.clip(preds["rDM_p"](mo), 0, 1)
    rl_c40 = np.clip(preds["rLR_c"](mo), 0, 1)
    rl_p40 = np.clip(preds["rLR_p"](mo), 0, 1)
    dmfs_c40 = os_c40 * rd_c40
    dmfs_p40 = os_p40 * rd_p40
    rfs_c40 = dmfs_c40 * rl_c40
    rfs_p40 = dmfs_p40 * rl_p40

    # Structural validity: RFS <= DMFS <= OS, monotonic
    valid = True
    for rfs, dmfs, osv in [(rfs_c40, dmfs_c40, os_c40), (rfs_p40, dmfs_p40, os_p40)]:
        if np.any(rfs > dmfs + 1e-9) or np.any(dmfs > osv + 1e-9):
            valid = False
    for arr in [rfs_c40, dmfs_c40, os_c40, rfs_p40, dmfs_p40, os_p40]:
        if np.any(np.diff(arr) > 1e-9):
            valid = False

    # Total AIC = sum of component AICs
    total_aic = sum(fits[(key, dist)][1] for key, dist in zip(["OS_c", "OS_p", "rDM_c", "rDM_p", "rLR_c", "rLR_p"], combo))
    return {
        "combo": combo, "valid": valid,
        "os40c": float(os_c40[-1]), "os40p": float(os_p40[-1]),
        "total_aic": total_aic,
        "rfs_c": rfs_c40, "dmfs_c": dmfs_c40, "os_c": os_c40,
    }

results = []
for oc in DISTS:
    for op in DISTS:
        for dc in DISTS:
            for dp in DISTS:
                for lc in DISTS:
                    for lp in DISTS:
                        results.append(build_and_score((oc, op, dc, dp, lc, lp)))

valid = [r for r in results if r["valid"]]
invalid = [r for r in results if not r["valid"]]
print(f"Total combinations: {len(results)}")
print(f"Structurally valid (40y, RFS<=DMFS<=OS, monotonic): {len(valid)}")
print(f"Invalid (crossing/monotonic violations): {len(invalid)}")
print(f"Clinical plausibility (OS40<10%): {sum(1 for r in valid if r['os40c'] < 0.10 and r['os40p'] < 0.10)}")

# Sort by total AIC (smaller is better), rank
valid_sorted = sorted(valid, key=lambda r: r["total_aic"])
print("\nTOP 20 by total AIC:")
print(f"{'Rank':<5s} {'OSc':<12s} {'OSp':<12s} {'rDMc':<10s} {'rDMp':<10s} {'rLRc':<10s} {'rLRp':<10s} {'TotalAIC':>9s} {'40yOK':>5s} {'OS40c':>7s} {'OS40p':>7s}")
for rank, r in enumerate(valid_sorted[:20], 1):
    gpok = r["os40c"] < 0.10 and r["os40p"] < 0.10
    print(f"{rank:<5d} {r['combo'][0]:<12s} {r['combo'][1]:<12s} {r['combo'][2]:<10s} {r['combo'][3]:<10s} {r['combo'][4]:<10s} {r['combo'][5]:<10s} {r['total_aic']:>9.2f} {'OK' if gpok else '--':>5s} {r['os40c']*100:>6.1f}% {r['os40p']*100:>6.1f}%")

# Where does the log-normal baseline rank?
ln_combo = ("lognormal",)*6
ln_rank = None
for i, r in enumerate(valid_sorted, 1):
    if r["combo"] == ln_combo:
        ln_rank = i
        break
print(f"\nLog-normal baseline (all 6 curves lognormal) rank by total AIC: {ln_rank}/{len(valid_sorted)}")
print(f"  Total AIC: {next(r['total_aic'] for r in valid_sorted if r['combo'] == ln_combo):.2f}")

# Save full results to CSV for GitHub
import csv
with open("/Users/cary/Documents/cea-intismeran/output/survival_combos_15625.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "OS_combo", "OS_pembro", "rDM_combo", "rDM_pembro", "rLR_combo", "rLR_pembro",
                "total_AIC", "structurally_valid_40y", "OS40_combo", "OS40_pembro", "GP_plausible"])
    for i, r in enumerate(valid_sorted, 1):
        w.writerow([i, *r["combo"], round(r["total_aic"], 3), int(r["valid"]),
                    round(r["os40c"], 5), round(r["os40p"], 5), int(r["os40c"] < 0.10 and r["os40p"] < 0.10)])
print(f"\nSaved full 15625-row results to output/survival_combos_15625.csv")