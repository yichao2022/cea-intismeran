"""Appendix I: Value-based pricing — price threshold table + Figure S8."""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import ModelParams
from rerun_primary import CEAModelV2

plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "figure.dpi": 300, "savefig.dpi": 300,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.1})

def icer_at_price(price):
    p = ModelParams(constraint_general_pop=True, cost_intismeran=price)
    m = CEAModelV2(p)
    r = m.run()
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    icer = dc / dq if dq > 0 else float('inf')
    nmb100 = -dc + 100_000 * dq
    nmb150 = -dc + 150_000 * dq
    return icer, nmb100, nmb150

prices = [0, 50_000, 100_000, 150_000, 200_000, 250_000, 300_000, 350_000, 400_000,
          429_077, 450_000, 500_000, 550_000, 600_000, 616_455, 650_000, 700_000, 750_000, 800_000]

print(f"{'Price':>10s} {'ICER':>10s} {'NMB@100K':>12s} {'NMB@150K':>12s}")
print("-" * 48)
for pr in prices:
    icer, nmb100, nmb150 = icer_at_price(pr)
    print(f"${pr:>8,.0f} ${icer:>8,.0f} ${nmb100:>10,.0f} ${nmb150:>10,.0f}")

# Solve for exact max prices
def find_price_for_icer(target, low=0, high=5_000_000, tol=100):
    for _ in range(60):
        mid = (low + high) / 2
        icer, _, _ = icer_at_price(mid)
        if abs(icer - target) < tol:
            break
        if icer < target:
            low = mid
        else:
            high = mid
    return mid, icer

p100, ic100 = find_price_for_icer(100_000)
p150, ic150 = find_price_for_icer(150_000)
print(f"\nMax price @ $100K/QALY: ${p100:,.0f} (ICER ${ic100:,.0f})")
print(f"Max price @ $150K/QALY: ${p150:,.0f} (ICER ${ic150:,.0f})")

# Figure S8
pgrid = np.arange(0, 800_001, 25_000)
icers = [icer_at_price(p)[0] for p in pgrid]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(pgrid / 1000, icers, "-", color="#2166ac", lw=2, label="ICER vs. acquisition price")
ax.axhline(100_000, color="#e41a1c", ls="--", lw=1, label="$100,000/QALY threshold")
ax.axhline(150_000, color="#4daf4a", ls="--", lw=1, label="$150,000/QALY threshold")
ax.axvline(p100 / 1000, color="#e41a1c", ls=":", lw=1, alpha=0.6)
ax.axvline(p150 / 1000, color="#4daf4a", ls=":", lw=1, alpha=0.6)
ax.scatter([p100 / 1000], [100_000], color="#e41a1c", s=50, zorder=10)
ax.scatter([p150 / 1000], [150_000], color="#4daf4a", s=50, zorder=10)
ax.annotate(f"${p100:,.0f}", (p100 / 1000, 100_000), textcoords="offset points", xytext=(-10, 12), fontsize=9, color="#e41a1c")
ax.annotate(f"${p150:,.0f}", (p150 / 1000, 150_000), textcoords="offset points", xytext=(-10, 12), fontsize=9, color="#4daf4a")
ax.set_xlabel("Intismeran acquisition price ($000s)")
ax.set_ylabel("ICER ($/QALY)")
ax.set_xlim(0, 800)
ax.set_ylim(0, 200_000)
ax.legend(fontsize=9, loc="upper left")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/tmp/figures/fig_s8_price_threshold.png", dpi=300)
plt.savefig("/tmp/figures/fig_s8_price_threshold.pdf")
print("Saved Figure S8")