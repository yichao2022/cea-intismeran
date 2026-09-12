"""Generate canonical results JSON — single source of truth for all tables & text."""
import sys, os, json
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from model import ModelParams
from rerun_primary import CEAModelV2

p = ModelParams(constraint_general_pop=True)
m = CEAModelV2(p)
r = m.run()
c, pn = r["combo"], r["pembro"]
t = m.t
disc = discount_factors(p.discount_rate, t)
cl = p.cycle_length
sur = m._survival()

ly_c = float(np.sum(sur["os_combo"] * disc) * cl)
ly_p = float(np.sum(sur["os_pembro"] * disc) * cl)
dq = float(c["qaly"] - pn["qaly"])
dc = float(c["cost"] - pn["cost"])
dly = float(ly_c - ly_p)
icer = round(dc / dq)

results = {
    "base_case": {
        "combo": {"qaly": round(c["qaly"], 2), "cost": round(c["cost"], 0), "ly": round(ly_c, 2)},
        "pembro": {"qaly": round(pn["qaly"], 2), "cost": round(pn["cost"], 0), "ly": round(ly_p, 2)},
        "incremental": {
            "dqaly": round(dq, 2),
            "dly": round(dly, 2),
            "dcost": round(dc, 0),
            "icer": icer,
            "nmb_100k": round(-dc + 100_000 * dq, 0),
            "nmb_150k": round(-dc + 150_000 * dq, 0),
        },
        "precision": {
            "dqaly_full": round(dq, 6),
            "dcost_full": round(dc, 4),
            "icer_full": round(dc/dq, 4),
        },
    }
}

print(json.dumps(results, indent=2))
with open("/Users/cary/cea-intismeran/output/base_case_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("Saved to output/base_case_results.json")