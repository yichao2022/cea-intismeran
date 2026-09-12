"""Tornado diagram (Figure 3) — generated from canonical output/dsa_v2_full.json.
Replaces figures.py plot_tornado (hardcoded 5-param subset with stale labels).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

with open("output/dsa_canonical.json") as f:
    dsa = json.load(f)

base_icer = dsa["base_icer"]
rows = dsa["rows"]

# Human-readable labels (keys are ModelParams attribute names)
LABELS = {
    "util_rf": "Utility: recurrence-free",
    "util_lr": "Utility: locoregional recurrence",
    "util_dm": "Utility: distant metastasis",
    "cost_keytruda_annual": "Pembrolizumab annual cost",
    "cost_intismeran": "Intismeran acquisition cost",
    "cost_lr_monthly": "LR management cost",
    "cost_dm_monthly": "DM management cost",
    "discount_rate": "Discount rate",
    "os_mu_combo": "OS log-mean, combination arm",
    "os_mu_pembro": "OS log-mean, pembrolizumab arm",
}

# Build tornado data: (label, low_delta, high_delta)
items = []
for r in rows:
    label = LABELS.get(r["param"], r["param"])
    lo = r.get("lo_icer_raw")
    hi = r.get("hi_icer_raw")
    # Dominant side (None) = no finite ICER; cap at plot max, annotate below
    items.append((label, lo, hi, r["lo_icer"] if lo is not None else None, r["hi_icer"] if hi is not None else None))

# Sort by max swing (largest at top); Dominant treated as full-width
def swing(it):
    lo, hi = it[1], it[2]
    v = []
    if lo is not None: v.append(abs(lo - base_icer))
    if hi is not None: v.append(abs(hi - base_icer))
    if lo is None or hi is None: v.append(1e12)  # dominant = biggest
    return max(v)
items.sort(key=swing)
items.reverse()

# Plot cap: max finite ICER value, plus headroom
finite_icers = [v for it in items for v in (it[1], it[2]) if v is not None]
xmax = max(max(finite_icers), base_icer) * 1.12

labels = [it[0] for it in items]
lows, highs = [], []
dominant_flags = []  # (side, is_dominant)
for it in items:
    label, lo, hi, lo_txt, hi_txt = it
    lows.append(None if lo is None else lo - base_icer)
    highs.append(None if hi is None else hi - base_icer)
    dominant_flags.append(("lo" if lo is None else None, "hi" if hi is None else None))

fig, ax = plt.subplots(figsize=(8, 6))
y = np.arange(len(items))
for i, (l, h) in enumerate(zip(lows, highs)):
    if l is not None:
        ax.barh(i, l, left=base_icer, color="#2166ac", alpha=0.85)
    else:
        ax.barh(i, xmax - base_icer, left=base_icer, color="#2166ac", alpha=0.45, hatch="//")
        ax.text(xmax + 2000, i, "Dominant", fontsize=8, va="center", color="#2166ac")
    if h is not None:
        ax.barh(i, h, left=base_icer, color="#b2182b", alpha=0.85)
    else:
        ax.barh(i, xmax - base_icer, left=base_icer, color="#b2182b", alpha=0.45, hatch="//")
        ax.text(xmax + 2000, i, "Dominant", fontsize=8, va="center", color="#b2182b")
ax.axvline(base_icer, color="black", lw=1.5, label=f"Base case: ${base_icer:,.0f}")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.invert_yaxis()  # largest swing at top (journal convention)
ax.set_xlabel("ICER (US$/QALY)")
ax.set_xlim(0, xmax + 45000)
ax.xaxis.set_major_formatter(lambda x, pos: f"${x/1000:.0f}K" if x else "0")
ax.legend(fontsize=8, loc="lower right")
ax.grid(axis="x", alpha=0.3)
plt.tight_layout()
plt.savefig("fig3_tornado.png", dpi=300)
plt.savefig("fig3_tornado.pdf")
plt.savefig("output/fig3_tornado.png", dpi=300)
plt.savefig("output/fig3_tornado.pdf")
print(f"Saved fig3_tornado.png/pdf ({len(items)} parameters, base ICER ${base_icer:,.0f})")