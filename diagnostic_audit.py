"""
Diagnostic audit of the 4-state PSM: verify $51,753 is a real result.
Outputs: 40-year diagnostic table + new Figure 1 with RFS/DMFS/OS.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np

from model import ModelParams, CEAModel

# ── Run model ──
p = ModelParams()
m = CEAModel(p)
sur = m._survival()
states = m._state_proportions(sur)
res = m.run()

# ── 1. Annual diagnostic table ──
print("=" * 130)
print(f"DIAGNOSTIC TABLE: 4-State PSM — Yearly Snapshot")
print(f"Model: log-normal OS, {p.time_horizon_years}y horizon, {p.discount_rate*100:.0f}% discount")
print(f"GP constraint: {p.constraint_general_pop}, Waning: {p.treatment_waning}")
print("=" * 130)
hdr = f"{'Yr':>3s} | {'Arm':>6s} | {'RFS':>6s} {'DMFS':>6s} {'OS':>6s} | {'RF':>6s} {'LR':>6s} {'DM':>6s} {'Dead':>6s} | {'Sum':>6s} | {'RFS≤DMFS':>8s} {'DMFS≤OS':>8s} {'Clamp?':>6s}"
print(hdr)
print("-" * 130)

t = m.t
n_years = int(p.time_horizon_years) + 1

artifacts = []
for yr in range(n_years):
    idx = int(yr / p.cycle_length)  # monthly index for this year
    if idx >= len(t):
        break
    for arm in ["combo", "pembro"]:
        prefix = "combo" if arm == "combo" else "pembro"
        rfs = sur[f"rfs_{arm}"][idx]
        dmfs = sur[f"dmfs_{arm}"][idx]
        os_val = sur[f"os_{arm}"][idx]
        st = states[arm]
        rf = st["rf"][idx]
        lr = st["lr"][idx]
        dm = st["dm"][idx]
        dead = st["dead"][idx]
        ssum = rf + lr + dm + dead

        # Checks
        rfs_ok = rfs <= dmfs + 1e-10
        dmfs_ok = dmfs <= os_val + 1e-10
        clamp = ""

        # Detect clamping: check if DMFS == OS or RFS == DMFS
        if abs(dmfs - os_val) < 1e-10:
            clamp = "DM=OS"
        if abs(rfs - dmfs) < 1e-10:
            clamp = "RF=DM" if clamp else "RF=DMFS"
        if abs(rf) < 1e-10 and yr > 0:
            clamp = "RF→0"

        if clamp:
            artifacts.append((yr, arm, clamp))

        flag = "+" if clamp else ""
        s = f"{yr:3d} | {arm:>6s} | {rfs:>6.4f} {dmfs:>6.4f} {os_val:>6.4f} | {rf:>6.4f} {lr:>6.4f} {dm:>6.4f} {dead:>6.4f} | {ssum:>6.4f} | {str(rfs_ok):>8s} {str(dmfs_ok):>8s} | {flag:>6s}"
        if yr in [0, 1, 2, 3, 4, 5, 10, 15, 20, 30, 40] or clamp:
            print(s)

print("-" * 130)

# ── 2. Summary statistics ──
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
ic = m.icer(res)
print(f"ICER: ${ic:,.0f}/QALY")
print(f"  Combo:  QALY={res['combo']['qaly']:.2f}  Cost=${res['combo']['cost']:,.0f}")
print(f"  Pembro: QALY={res['pembro']['qaly']:.2f}  Cost=${res['pembro']['cost']:,.0f}")
print(f"  ΔQALY={res['combo']['qaly']-res['pembro']['qaly']:.2f}  ΔCost=${res['combo']['cost']-res['pembro']['cost']:,.0f}")

# ── 3. Key year OS ──
print(f"\n{'='*60}")
print("OS AT KEY TIME POINTS")
print(f"{'='*60}")
print(f"{'Year':>5s} {'Combo OS':>10s} {'Pembro OS':>10s} {'Diff':>10s}")
for yr in [5, 10, 15, 20, 30, 40]:
    idx = int(yr / p.cycle_length)
    if idx < len(t):
        cos = sur[f"os_combo"][idx]
        pos = sur[f"os_pembro"][idx]
        print(f"{yr:5d} {cos:>10.4f} {pos:>10.4f} {cos-pos:>10.4f}")

# ── 4. Estimated life-years ──
print(f"\n{'='*60}")
print("LIFE-YEARS (undiscounted OS integral)")
print(f"{'='*60}")
for arm in ["combo", "pembro"]:
    ly = np.sum(sur[f"os_{arm}"]) * p.cycle_length
    print(f"  {arm}: {ly:.2f} life-years")

# ── 5. Check for structural issues ──
print(f"\n{'='*60}")
print("STRUCTURAL INTEGRITY CHECK")
print(f"{'='*60}")
if artifacts:
    print(f"Artifacts found: {len(artifacts)} instances")
    for yr, arm, issue in artifacts:
        print(f"  Year {yr}, {arm}: {issue}")
else:
    print("No clamping/forced equality detected")

# Check monotonicity
print(f"\nMonotonicity checks:")
for arm in ["combo", "pembro"]:
    for curve in ["rfs", "dmfs", "os"]:
        vals = sur[f"{curve}_{arm}"]
        is_mono = np.all(np.diff(vals) <= 1e-10)
        if not is_mono:
            violations = np.where(np.diff(vals) > 1e-10)[0]
            print(f"  ⚠ {arm} {curve}: NOT monotonic — {len(violations)} violations at idx {violations[:5]}")
        else:
            print(f"  ✓ {arm} {curve}: monotonic")

# Check state sum = 1
max_dev = 0
for arm in ["combo", "pembro"]:
    st = states[arm]
    ssum = st["rf"] + st["lr"] + st["dm"] + st["dead"]
    dev = np.max(np.abs(ssum - 1.0))
    max_dev = max(max_dev, dev)
    if dev > 1e-10:
        print(f"  ⚠ {arm}: max state sum deviation = {dev:.2e}")
print(f"  Max state sum deviation: {max_dev:.2e} {'✓' if max_dev < 1e-10 else '⚠'}")

# ── 6. New Figure 1 with RFS, DMFS, OS ──
print(f"\n{'='*60}")
print("GENERATING NEW FIGURE 1 (RFS + DMFS + OS)")
print(f"{'='*60}")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, arm, label, color in [(axes[0], "combo", "Combo", "#2196F3"),
                                (axes[1], "pembro", "Pembrolizumab", "#FF5722")]:
    ax.plot(t, sur[f"rfs_{arm}"], label="RFS", color=color, ls="--", lw=1.5)
    ax.plot(t, sur[f"dmfs_{arm}"], label="DMFS", color=color, ls="-.", lw=1.5)
    ax.plot(t, sur[f"os_{arm}"], label="OS", color=color, ls="-", lw=2)
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Years")
    ax.set_ylabel("Survival Probability")
    ax.set_title(label)
    ax.legend(loc="lower left")
    ax.grid(True, alpha=0.3)

    # Mark data points
    p = m.p
    data_yrs = [18/12, 24/12, 36/12, 48/12, 60/12]
    if arm == "combo":
        rfs_data = [v/100 for v in p.rfs_combo_pct]
        os_data = [v/100 for v in p.os_combo_pct]
        dmfs_data = [v/100 for v in p.dmfs_combo_pct] + [None]
        dmfs_yrs = [18/12, 24/12, 36/12, 48/12, None]
    else:
        rfs_data = [v/100 for v in p.rfs_pembro_pct]
        os_data = [v/100 for v in p.os_pembro_pct]
        dmfs_data = [v/100 for v in p.dmfs_pembro_pct] + [None]
        dmfs_yrs = [18/12, 24/12, 36/12, 48/12, None]

    ax.scatter(data_yrs, rfs_data, color=color, marker="o", s=30, zorder=5)
    ax.scatter(data_yrs, os_data, color=color, marker="s", s=30, zorder=5)
    dmfs_data_clean = [v for v in dmfs_data if v is not None]
    dmfs_yrs_clean = [v for v in dmfs_yrs if v is not None]
    if dmfs_data_clean:
        ax.scatter(dmfs_yrs_clean, dmfs_data_clean, color=color, marker="^", s=30, zorder=5)

fig.suptitle("Figure 1: Partitioned Survival Model — RFS, DMFS, and OS Curves\n(Log-normal extrapolation, 40-year horizon)", fontsize=12)
plt.tight_layout()
plt.savefig("/Users/cary/cea-intismeran/output/fig1_diagnostic_survival.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: output/fig1_diagnostic_survival.png")

# ── 7. Cost breakdown ──
print(f"\n{'='*60}")
print("COST BREAKDOWN")
print(f"{'='*60}")
for arm in ["combo", "pembro"]:
    r = res[arm]
    print(f"  {arm}: Total=${r['cost']:,.0f}")

print(f"\n{'='*60}")
print("DIAGNOSTIC COMPLETE")
print(f"{'='*60}")