"""Full DSA — canonical source for Figure 3, Table S12, and main Results sentence.
Runs V2 model over all Table S12 parameter ranges, writes JSON with ICER + NMB@150K.
"""
import json
import numpy as np
from model import ModelParams
from rerun_primary import CEAModelV2

WTP = 150_000

params = ModelParams(constraint_general_pop=True)
base = CEAModelV2(params).run()
base_icer = (base["combo"]["cost"] - base["pembro"]["cost"]) / (base["combo"]["qaly"] - base["pembro"]["qaly"])
print(f"base_icer = {base_icer:.2f}\n")

# Generic runner: returns (dC, dQ) under a perturbed param set
def run_perturbed(**kwargs):
    pp = ModelParams(constraint_general_pop=True)
    for k, v in kwargs.items():
        setattr(pp, k, v)
    r = CEAModelV2(pp).run()
    dC = r["combo"]["cost"] - r["pembro"]["cost"]
    dQ = r["combo"]["qaly"] - r["pembro"]["qaly"]
    return dC, dQ

def run_os_mu(which, factor):
    pp = ModelParams(constraint_general_pop=True)
    if which == "combo":
        pp.os_mu_combo *= factor
    else:
        pp.os_mu_pembro *= factor
    r = CEAModelV2(pp).run()
    dC = r["combo"]["cost"] - r["pembro"]["cost"]
    dQ = r["combo"]["qaly"] - r["pembro"]["qaly"]
    return dC, dQ

def to_row(param, base_disp, low_lbl, high_lbl, lo_dC, lo_dQ, hi_dC, hi_dQ):
    lo_icer = lo_dC / lo_dQ if lo_dQ > 0 else None
    hi_icer = hi_dC / hi_dQ if hi_dQ > 0 else None
    lo_nmb = WTP * lo_dQ - lo_dC
    hi_nmb = WTP * hi_dQ - hi_dC
    row = {
        "param": param, "base": base_disp,
        "low_lbl": low_lbl, "high_lbl": high_lbl,
        "lo_icer": "Dominant" if lo_icer is None or lo_icer < 0 else f"{round(lo_icer):,}",
        "hi_icer": "Dominant" if hi_icer is None or hi_icer < 0 else f"{round(hi_icer):,}",
        "lo_icer_raw": None if lo_icer is None else round(lo_icer),
        "hi_icer_raw": None if hi_icer is None else round(hi_icer),
        "lo_nmb": f"{round(lo_nmb):,}",
        "hi_nmb": f"{round(hi_nmb):,}",
        "change_nmb": f"{abs(round(hi_nmb - lo_nmb)):,}",
    }
    print(f"  {param:28s} lo={row['lo_icer']:>9s} hi={row['hi_icer']:>9s} |dNMB|={row['change_nmb']:>9s}")
    return row

t12_rows = []
print("=== Table S12 rows (V2) ===")

# Utilities — beta ±0.03
for key, base_disp, lo_v, hi_v, lo_lbl, hi_lbl in [
    ("util_rf", "0.83", 0.80, 0.86, "0.80", "0.86"),
    ("util_lr", "0.64", 0.61, 0.67, "0.61", "0.67"),
    ("util_dm", "0.55", 0.52, 0.58, "0.52", "0.58"),
]:
    lo_dC, lo_dQ = run_perturbed(**{key: lo_v})
    hi_dC, hi_dQ = run_perturbed(**{key: hi_v})
    t12_rows.append(to_row(key.upper(), base_disp, lo_lbl, hi_lbl, lo_dC, lo_dQ, hi_dC, hi_dQ))

# Costs
for key, base_disp, lo_v, hi_v, lo_lbl, hi_lbl in [
    ("cost_keytruda_annual", "$220,896", 176_717, 265_075, "$176,717", "$265,075"),
    ("cost_intismeran",      "$200,000", 160_000, 240_000, "$160,000", "$240,000"),
    ("cost_lr_monthly",      "$3,000",   2_400,   3_600,   "$2,400",   "$3,600"),
    ("cost_dm_monthly",      "$12,000",  9_600,   14_400,  "$9,600",   "$14,400"),
]:
    lo_dC, lo_dQ = run_perturbed(**{key: lo_v})
    hi_dC, hi_dQ = run_perturbed(**{key: hi_v})
    t12_rows.append(to_row(key, base_disp, lo_lbl, hi_lbl, lo_dC, lo_dQ, hi_dC, hi_dQ))

# Discount rate 0% / 5%
lo_dC, lo_dQ = run_perturbed(discount_rate=0.00)
hi_dC, hi_dQ = run_perturbed(discount_rate=0.05)
t12_rows.append(to_row("discount_rate", "3%", "0%", "5%", lo_dC, lo_dQ, hi_dC, hi_dQ))

# OS mu combo ±40%
lo_dC, lo_dQ = run_os_mu("combo", 0.6)
hi_dC, hi_dQ = run_os_mu("combo", 1.4)
t12_rows.append(to_row("OS mu combo", "8.396", "-40%", "+40%", lo_dC, lo_dQ, hi_dC, hi_dQ))

# OS mu pembro ±20%
lo_dC, lo_dQ = run_os_mu("pembro", 0.8)
hi_dC, hi_dQ = run_os_mu("pembro", 1.2)
t12_rows.append(to_row("OS mu pembro", "4.841", "-20%", "+20%", lo_dC, lo_dQ, hi_dC, hi_dQ))

out = {"base_icer": round(base_icer), "wtp": WTP, "t12_rows": t12_rows}
with open("output/dsa_v2_full.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved output/dsa_v2_full.json")
