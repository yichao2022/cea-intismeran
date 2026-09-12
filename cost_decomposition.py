"""
Cost decomposition: base case (GP constrained) vs waning scenario.
Per arm, per cost component.
"""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
from model import ModelParams
from model import discount_factors  # noqa: E402
from rerun_primary import CEAModelV2

def cost_decomposition(p: ModelParams, label: str):
    """Return per-component costs for both arms."""
    m = CEAModelV2(p)
    t = m.t
    sur = m._survival()
    states = m._state_proportions(sur)
    disc = discount_factors(p.discount_rate, t)
    cl = p.cycle_length

    print(f"\n{'='*95}")
    print(f"COST DECOMPOSITION: {label}")
    print(f"{'='*95}")

    components = ["Intismeran", "Pembrolizumab", "Sequencing", "Administration",
                  "Adverse events", "LR management", "DM management"]
    rows = {arm: {} for arm in ["combo", "pembro"]}

    for arm in ["combo", "pembro"]:
        s = states[arm]
        rows[arm]["Intismeran"] = (p.cost_intismeran * disc[0] * cl * 12) if arm == "combo" else 0.0
        rows[arm]["Pembrolizumab"] = p.cost_keytruda_annual * disc[:12].sum() * cl
        rows[arm]["Sequencing"] = (p.cost_sequencing * disc[0]) if arm == "combo" else 0.0
        rows[arm]["Administration"] = p.cost_admin_iv_infusion * float(discount_factors(p.discount_rate, np.arange(18) * (3 / 52)).sum())  # 18 IV infusions q3w
        if arm == "combo":   # 9 intismeran doses are IM injections, not infusions
            rows[arm]["Administration"] += p.cost_admin_im_injection * float(discount_factors(p.discount_rate, np.arange(9) * (3 / 52)).sum())
        rows[arm]["Adverse events"] = (p.cost_ae_incremental * disc[0]) if arm == "combo" else 0.0
        rows[arm]["LR management"] = np.sum(s["lr"] * p.cost_lr_monthly * 12 * cl * disc)
        rows[arm]["DM management"] = np.sum(s["dm"] * p.cost_dm_monthly * 12 * cl * disc)

    # Print table
    hdr = f"{'Cost component':<22s} {'Combo':>12s} {'Pembro':>12s} {'Incremental':>12s}"
    print(hdr)
    print("-" * 95)
    total_c = total_p = 0.0
    for comp in components:
        c = rows["combo"][comp]
        pv = rows["pembro"][comp]
        total_c += c
        total_p += pv
        print(f"{comp:<22s} ${c:>9,.0f} ${pv:>9,.0f} ${c-pv:>9,.0f}")
    total_inc = total_c - total_p
    print("-" * 95)
    print(f"{'TOTAL':<22s} ${total_c:>9,.0f} ${total_p:>9,.0f} ${total_inc:>9,.0f}")

    # QALYs
    qc = np.sum((states["combo"]["rf"]*p.util_rf + states["combo"]["lr"]*p.util_lr + states["combo"]["dm"]*p.util_dm) * disc) * cl
    qp = np.sum((states["pembro"]["rf"]*p.util_rf + states["pembro"]["lr"]*p.util_lr + states["pembro"]["dm"]*p.util_dm) * disc) * cl

    icer = total_inc / (qc - qp) if (qc - qp) > 0 else float('inf')
    print(f"\nCombo QALY: {qc:.2f} | Pembro QALY: {qp:.2f} | ΔQALY: {qc-qp:.2f}")
    print(f"ICER = ${total_inc:,.0f} / {qc-qp:.2f} = ${icer:,.0f}/QALY")

    return {"total_c": total_c, "total_p": total_p, "inc": total_inc,
            "qc": qc, "qp": qp, "icer": icer}

# Base case (GP constrained)
p1 = ModelParams(constraint_general_pop=True)
r1 = cost_decomposition(p1, "BASE CASE (GP constrained)")

# Waning
p2 = ModelParams(constraint_general_pop=True, treatment_waning=True)
r2 = cost_decomposition(p2, "WANING SCENARIO (GP constrained)")

# Compare
print(f"\n{'='*95}")
print("BASE vs WANING EXPLANATION")
print(f"{'='*95}")
print(f"{'':>22s} {'Base':>12s} {'Waning':>12s} {'Change':>12s}")
for field, label, fmt in [("inc", "ΔCost", "${:,.0f}"), ("qc", "Combo QALY", "{:.2f}"),
                           ("qp", "Pembro QALY", "{:.2f}"), ("icer", "ICER", "${:,.0f}")]:
    v1, v2 = r1[field], r2[field]
    print(f"{label:<22s} {fmt.format(v1):>12s} {fmt.format(v2):>12s} {fmt.format(v2-v1):>12s}")

d_cost = r2["inc"] - r1["inc"]
d_qaly = (r2["qc"]-r2["qp"]) - (r1["qc"]-r1["qp"])
print(f"\nΔCost change: ${d_cost:,.0f} ({d_cost/max(r1['inc'],1e-9)*100:.0f}% of base ΔCost)")
print(f"ΔQALY change: {d_qaly:.2f} ({d_qaly/(r1['qc']-r1['qp'])*100:.0f}% of base ΔQALY)")
print(f"\n→ Waning cuts ΔCost by {abs(d_cost)/max(r1['inc'],1e-9)*100:.0f}% but ΔQALY by only {abs(d_qaly)/(r1['qc']-r1['qp'])*100:.0f}%")
print(f"  (Cost drops MORE than QALY → ICER falls)")
