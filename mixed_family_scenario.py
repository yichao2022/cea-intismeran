"""
Best-joint-AIC mixed-family scenario (Rank 1 combination from Table S4):
  OS:      combo=lognormal, pembro=weibull
  rDM:     combo=lognormal, pembro=loglogistic
  rLR:     combo=loglogistic, pembro=weibull
Runs the model with these per-component distributions and reports ICER.
"""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm
from model import ModelParams, CEAModel, lognorm_surv, weibull_surv, loglogistic_surv

p0 = ModelParams()
def pct(v): return np.array(v, dtype=float) / 100.0

t_os = np.array(p0.os_time_mo)
t_dmfs = np.array(p0.dmfs_time_mo)
t_rfs = np.array(p0.rfs_time_mo)

os_combo_pct = pct(p0.os_combo_pct)
os_pembro_pct = pct(p0.os_pembro_pct)
dmfs_combo_pct = pct(p0.dmfs_combo_pct)
dmfs_pembro_pct = pct(p0.dmfs_pembro_pct)
rfs_combo_pct = pct(p0.rfs_combo_pct)
rfs_pembro_pct = pct(p0.rfs_pembro_pct)

os_combo_at_dmfs = np.interp(t_dmfs, t_os, os_combo_pct)
os_pembro_at_dmfs = np.interp(t_dmfs, t_os, os_pembro_pct)
dmfs_combo_at_rfs = np.interp(t_rfs, t_dmfs, dmfs_combo_pct)
dmfs_pembro_at_rfs = np.interp(t_rfs, t_dmfs, dmfs_pembro_pct)

# ratios
rd_c = np.clip(dmfs_combo_pct / np.maximum(os_combo_at_dmfs, 1e-12), 0.01, 0.99)
rd_p = np.clip(dmfs_pembro_pct / np.maximum(os_pembro_at_dmfs, 1e-12), 0.01, 0.99)
rl_c = np.clip(rfs_combo_pct / np.maximum(dmfs_combo_at_rfs, 1e-12), 0.01, 0.99)
rl_p = np.clip(rfs_pembro_pct / np.maximum(dmfs_pembro_at_rfs, 1e-12), 0.01, 0.99)

def fit(dstr, t, s, pname):
    if dstr == "lognormal":
        f = lambda p: norm.cdf((p[0] - np.log(t)) / p[1]) - s
        x0, bnd = [np.log(t[len(t)//2]), 0.8], ([0.1, 0.01], [15, 5])
    elif dstr == "weibull":
        f = lambda p: np.exp(-(t / p[0]) ** p[1]) - s
        x0, bnd = [t[-1]/(-np.log(max(s[-1],1e-6)))**(1/1.2), 1.2], ([10, 0.2], [5000, 5.0])
    elif dstr == "loglogistic":
        f = lambda p: 1.0 / (1.0 + (t / p[1]) ** p[0]) - s
        x0, bnd = [1.5, 50.0], ([0.05, 1], [15, 5000])
    sol = least_squares(f, x0=x0, bounds=bnd, xtol=1e-12, ftol=1e-12, max_nfev=20000)
    return sol.x

# Rank-1 combination
fits = {
    "os_combo": fit("lognormal", t_os, os_combo_pct, "os_combo"),
    "os_pembro": fit("weibull", t_os, os_pembro_pct, "os_pembro"),
    "rd_combo": fit("lognormal", t_dmfs, rd_c, "rd_combo"),
    "rd_pembro": fit("loglogistic", t_dmfs, rd_p, "rd_pembro"),
    "rl_combo": fit("loglogistic", t_rfs, rl_c, "rl_combo"),
    "rl_pembro": fit("weibull", t_rfs, rl_p, "rl_pembro"),
}
print("Fitted parameters:")
for k, v in fits.items():
    print(f"  {k}: {np.round(v, 4)}")

class CEAModelMixed(CEAModel):
    """Mixed-family: per-component distributions (Rank-1 AIC combination)."""
    def _survival(self):
        mo = self.t * 12
        p = self.p
        sur = {}
        for arm in ["combo", "pembro"]:
            prefix = "combo" if arm == "combo" else "pembro"
            # OS
            if prefix == "combo":
                os_val = lognorm_surv(mo, fits["os_combo"][0], fits["os_combo"][1])
            else:
                os_val = weibull_surv(mo, fits["os_pembro"][0], fits["os_pembro"][1])
            if p.constraint_general_pop:
                gp_surv = general_pop_surv(self.t)[:len(mo)]
                os_val = np.minimum(os_val, gp_surv)   # ponytail: V2 hazard floor is in CEAModelV2; here min-cap to test AIC-sensitivity only
            sur[f"os_{arm}"] = os_val
            # rDM
            if prefix == "combo":
                r_dm = lognorm_surv(mo, fits["rd_combo"][0], fits["rd_combo"][1])
            else:
                r_dm = loglogistic_surv(mo, fits["rd_pembro"][0], fits["rd_pembro"][1])
            sur[f"dmfs_{arm}"] = os_val * r_dm
            # rLR
            if prefix == "combo":
                r_lr = loglogistic_surv(mo, fits["rl_combo"][0], fits["rl_combo"][1])
            else:
                r_lr = weibull_surv(mo, fits["rl_pembro"][0], fits["rl_pembro"][1])
            sur[f"rfs_{arm}"] = sur[f"dmfs_{arm}"] * r_lr
        # monotonic + ordering enforcement
        for arm in ["combo", "pembro"]:
            for key in [f"os_{arm}", f"dmfs_{arm}", f"rfs_{arm}"]:
                arr = sur[key]
                for i in range(1, len(arr)):
                    arr[i] = min(arr[i], arr[i-1])
            sur[f"rfs_{arm}"] = np.minimum(sur[f"rfs_{arm}"], sur[f"dmfs_{arm}"])
            sur[f"dmfs_{arm}"] = np.minimum(sur[f"dmfs_{arm}"], sur[f"os_{arm}"])
        return sur

from model import general_pop_surv

for gp in [False, True]:
    p = ModelParams(constraint_general_pop=gp)
    m = CEAModelMixed(p)
    r = m.run()
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    icer = dc / dq if dq > 0 else float("inf")
    print(f"\nCOMBINATION (Rank-1 AIC mixed): GP={gp}")
    print(f"  ΔQALY={dq:.2f}  ΔCost=${dc:,.0f}  ICER=${icer:,.0f}")