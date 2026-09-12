"""Canonical S10 (cost decomposition) + S11 (LY/QALY decomposition).
From V2 model state proportions. Matches model.py discount formula exactly.
Single source of truth → output/decomposition_v2.json.
"""
import json, sys, io
import numpy as np
from model import ModelParams
from rerun_primary import CEAModelV2

params = ModelParams(constraint_general_pop=True)
m = CEAModelV2(params)
old = sys.stdout; sys.stdout = io.StringIO()
r = m.run()
sur = m._survival()
states = m._state_proportions(sur)
sys.stdout = old

# Match model.py line 336: disc = np.exp(-p.discount_rate * self.t)
disc = np.exp(-params.discount_rate * m.t)

out = {"ly_qaly": {}, "cost_decomp": {}}
for arm in ["combo", "pembro"]:
    s = states[arm]
    rf, lr, dm = s["rf"], s["lr"], s["dm"]
    alive = rf + lr + dm

    ly_undisc = float(alive.sum() * params.cycle_length)
    ly_disc = float((alive * disc).sum() * params.cycle_length)

    qaly_rf = float((rf * disc).sum() * params.cycle_length * params.util_rf)
    qaly_lr = float((lr * disc).sum() * params.cycle_length * params.util_lr)
    qaly_dm = float((dm * disc).sum() * params.cycle_length * params.util_dm)
    qaly_ae = params.util_disutility_ae * disc[0] * params.cycle_length
    total_qaly = qaly_rf + qaly_lr + qaly_dm - qaly_ae

    keytruda = params.cost_keytruda_annual * disc[:12].sum() * params.cycle_length
    intismeran = params.cost_intismeran * disc[0] if arm == "combo" else 0
    sequencing = params.cost_sequencing * disc[0] if arm == "combo" else 0
    admin = params.cost_admin_per_cycle * float(np.exp(-params.discount_rate * np.arange(18) * (3 / 52)).sum())  # 18 infusions q3w
    lr_cost = float((lr * disc).sum() * params.cycle_length * params.cost_lr_monthly * 12)
    dm_cost = float((dm * disc).sum() * params.cycle_length * params.cost_dm_monthly * 12)
    ae_cost = params.cost_ae_incremental * disc[0] if arm == "combo" else 0
    total_cost = keytruda + intismeran + sequencing + admin + lr_cost + dm_cost + ae_cost

    out["ly_qaly"][arm] = {
        "ly_undisc": round(ly_undisc, 2), "ly_disc": round(ly_disc, 2),
        "qaly_rf": round(qaly_rf, 2), "qaly_lr": round(qaly_lr, 2),
        "qaly_dm": round(qaly_dm, 2), "qaly_total": round(total_qaly, 2),
    }
    out["cost_decomp"][arm] = {
        "pembrolizumab": round(keytruda), "intismeran": round(intismeran),
        "sequencing": round(sequencing), "administration": round(admin),
        "adverse_events": round(ae_cost),
        "lr_management": round(lr_cost), "dm_management": round(dm_cost),
        "total": round(total_cost),
    }

for cat in ["ly_qaly", "cost_decomp"]:
    out[cat]["incremental"] = {}
    for k in out[cat]["combo"]:
        out[cat]["incremental"][k] = round(out[cat]["combo"][k] - out[cat]["pembro"][k], 2)

with open("output/decomposition_v2.json", "w") as f:
    json.dump(out, f, indent=2)

print("=== S11 LY/QALY Decomposition ===")
for arm in ["combo", "pembro"]:
    d = out["ly_qaly"][arm]
    print(f"  {arm:8s} LY(ud)={d['ly_undisc']:>6.2f}  LY(d)={d['ly_disc']:>6.2f}  RF={d['qaly_rf']:.2f}  LR={d['qaly_lr']:.2f}  DM={d['qaly_dm']:.2f}  Tot={d['qaly_total']:.2f}")
inc = out["ly_qaly"]["incremental"]
print(f"  {'inc':8s} LY(ud)={inc['ly_undisc']:>6.2f}  LY(d)={inc['ly_disc']:>6.2f}  RF={inc['qaly_rf']:.2f}  LR={inc['qaly_lr']:.2f}  DM={inc['qaly_dm']:.2f}  Tot={inc['qaly_total']:.2f}")

print("\n=== S10 Cost Decomposition ===")
for k in out["cost_decomp"]["combo"]:
    c, p, d = out["cost_decomp"]["combo"][k], out["cost_decomp"]["pembro"][k], out["cost_decomp"]["incremental"][k]
    print(f"  {k:20s} combo=${c:>10,}  pembro=${p:>10,}  inc=${d:>10,}")

base_ic = out["cost_decomp"]["incremental"]["total"]
base_dq = out["ly_qaly"]["incremental"]["qaly_total"]
print(f"\nCross-check: dCost=${base_ic:,}  dQALY={base_dq}  ICER=${base_ic/base_dq:,.0f}")
print("Saved output/decomposition_v2.json")
