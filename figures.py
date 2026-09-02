"""
Generate figures for CEA manuscript.
Uses current CEAModelV2 API.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import contextlib
import io
from model import ModelParams
from rerun_primary import CEAModelV2

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.titlesize": 13, "axes.labelsize": 12,
    "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.1,
})
C = {"combo": "#2166ac", "pembro": "#b2182b"}

# ── 1. Survival Curves ──
def plot_survival():
    params = ModelParams()
    model = CEAModelV2(params)
    sur = model._survival()
    t = model.t

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    for arm_name in ["combo", "pembro"]:
        label = "Intismeran + Pembro" if arm_name == "combo" else "Pembrolizumab"
        # RFS = RFS curve
        rfs = sur[f"rfs_{arm_name}"]
        os = sur[f"os_{arm_name}"]
        ax1.plot(t, rfs, color=C[arm_name], label=label)
        ax2.plot(t, os, color=C[arm_name], label=label)

    # 5-year data points from KEYNOTE-942
    ax1.plot(5, 0.688, "o", color=C["combo"], ms=6)
    ax1.plot(5, 0.491, "o", color=C["pembro"], ms=6)
    ax2.plot(5, 0.922, "o", color=C["combo"], ms=6)
    ax2.plot(5, 0.713, "o", color=C["pembro"], ms=6)

    for ax, title in [(ax1, "Recurrence-Free Survival"), (ax2, "Overall Survival")]:
        ax.set_xlabel("Years"); ax.set_ylabel("Survival Probability")
        ax.set_title(title); ax.legend(fontsize=9)
        ax.set_xlim(0, 40); ax.set_ylim(0, 1); ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig("output/fig1_survival.pdf"); plt.savefig("output/fig1_survival.png")
    plt.close(); print("OK fig1")

# ── 2. CEAC ──
def plot_ceac():
    # Read canonical PSA draws (6000 iterations) to ensure consistency with manuscript text
    import json
    try:
        draws = np.load("psa_v2_draws.npz")
        dq = draws["dq"]
        dc = draws["dc"]
        icers = dc[dq > 0] / dq[dq > 0]
    except FileNotFoundError:
        from scipy.stats import gamma, norm
        rng = np.random.default_rng(42)
        icers_list = []
        for _ in range(6000):
            pp = ModelParams()
            for attr, val in [("cost_keytruda_annual", 220_896), ("cost_intismeran", 200_000),
                              ("cost_lr_monthly", 3_000), ("cost_dm_monthly", 12_000)]:
                setattr(pp, attr, rng.gamma(25, val/25))
            for attr, val in [("util_rf", 0.83), ("util_lr", 0.64), ("util_dm", 0.55)]:
                setattr(pp, attr, np.clip(rng.normal(val, 0.03), 0, 1))
            pp.os_mu_combo += rng.normal(0, 0.1)
            pp.os_mu_pembro += rng.normal(0, 0.1)
            with contextlib.redirect_stdout(io.StringIO()):
                r = CEAModelV2(pp).run()
            dc = r["combo"]["cost"] - r["pembro"]["cost"]
            dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
            if dq > 0: icers_list.append(dc/dq)
        icers = np.array(icers_list)

    thresholds = np.arange(0, 400_001, 10_000)
    prob_ce = [np.mean(icers <= t) for t in thresholds]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(thresholds / 1000, prob_ce, color="#2166ac", lw=2)
    ax.axvline(100, color="gray", ls="--", lw=1, label="$100K/QALY")
    ax.axvline(150, color="gray", ls=":", lw=1, label="$150K/QALY")
    ax.set_xlabel("Willingness-to-Pay Threshold ($1000s/QALY)")
    ax.set_ylabel("Probability Cost-Effective")
    ax.set_title("Cost-Effectiveness Acceptability Curve")
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/fig2_ceac.pdf"); plt.savefig("output/fig2_ceac.png")
    plt.close(); print("OK fig2")

# ── 3. Tornado ──
def plot_tornado():
    params = ModelParams(constraint_general_pop=True)
    base = CEAModelV2(params).run()
    base_icer = (base["combo"]["cost"] - base["pembro"]["cost"]) / (base["combo"]["qaly"] - base["pembro"]["qaly"])

    tornado = {}
    # Each: (param_getter, low_factor, high_factor, label)
    for name, getter, low, high, label in [
        ("rfs_combo", lambda p: None, 0.9, 1.1, "5-year RFS, combo"),
        ("cost_intismeran", lambda p: p.cost_intismeran, 0.5, 2.0, "Intismeran cost"),
        ("cost_dm_monthly", lambda p: p.cost_dm_monthly, 0.5, 1.5, "DM monthly cost"),
        ("util_rf", lambda p: p.util_rf, 0.9, 1.1, "Utility RF"),
        ("discount_rate", lambda p: p.discount_rate, 0.0, 0.05, "Discount rate"),
    ]:
        icers = []
        for factor in [low, high]:
            pp = ModelParams(constraint_general_pop=True)
            if name == "rfs_combo":
                # Perturb RFS by adjusting mu
                pp.os_mu_combo += np.log(1/factor) if factor < 1 else np.log(factor)
            elif name == "cost_intismeran":
                pp.cost_intismeran = 200_000 * factor
            elif name == "cost_dm_monthly":
                pp.cost_dm_monthly = 12_000 * factor
            elif name == "util_rf":
                pp.util_rf = np.clip(0.83 * factor, 0, 1)
            elif name == "discount_rate":
                pp.discount_rate = factor
            r = CEAModelV2(pp).run()
            icer = (r["combo"]["cost"] - r["pembro"]["cost"]) / (r["combo"]["qaly"] - r["pembro"]["qaly"])
            icers.append(icer)
        tornado[label] = (base_icer - icers[0], icers[1] - base_icer)

    # Sort by magnitude
    sorted_items = sorted(tornado.items(), key=lambda x: max(abs(x[1][0]), abs(x[1][1])))
    fig, ax = plt.subplots(figsize=(8, 5))
    y_pos = range(len(sorted_items))
    low_vals = [v[0] for _, v in sorted_items]
    high_vals = [v[1] for _, v in sorted_items]
    ax.barh(y_pos, low_vals, left=base_icer, color="#2166ac", alpha=0.8)
    ax.barh(y_pos, high_vals, left=base_icer, color="#b2182b", alpha=0.8)
    ax.axvline(base_icer, color="black", lw=1.5)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels([l for l, _ in sorted_items])
    ax.set_xlabel("ICER ($/QALY)")
    ax.set_title("One-Way Sensitivity Analysis (Tornado)")
    plt.tight_layout()
    plt.savefig("output/fig3_tornado.pdf"); plt.savefig("output/fig3_tornado.png")
    plt.close(); print("OK fig3")

# ── 4. CE Plane ──
def plot_ce_plane():
    from scipy.stats import gamma, norm
    rng = np.random.default_rng(42)
    delta_costs, delta_qalys = [], []
    for _ in range(1000):
        pp = ModelParams()
        for attr, val in [("cost_keytruda_annual", 220_896), ("cost_intismeran", 200_000),
                          ("cost_lr_monthly", 3_000), ("cost_dm_monthly", 12_000)]:
            setattr(pp, attr, rng.gamma(25, val/25))
        for attr, val in [("util_rf", 0.83), ("util_lr", 0.64), ("util_dm", 0.55)]:
            setattr(pp, attr, np.clip(rng.normal(val, 0.03), 0, 1))
        pp.os_mu_combo += rng.normal(0, 0.1)
        pp.os_mu_pembro += rng.normal(0, 0.1)
        r = CEAModelV2(pp).run()
        dc = r["combo"]["cost"] - r["pembro"]["cost"]
        dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
        if dq > 0:
            delta_costs.append(dc)
            delta_qalys.append(dq)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(delta_qalys, delta_costs, alpha=0.4, s=10, color="#2166ac")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    wtp = np.linspace(0, max(delta_qalys), 100)
    for w, label in [(100_000, "$100K/QALY"), (150_000, "$150K/QALY")]:
        ax.plot(wtp, w * wtp, "--", lw=1, label=label)
    ax.set_xlabel("Incremental QALYs"); ax.set_ylabel("Incremental Cost ($)")
    ax.set_title("Cost-Effectiveness Plane")
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/fig4_ce_plane.pdf"); plt.savefig("output/fig4_ce_plane.png")
    plt.close(); print("OK fig4")

if __name__ == "__main__":
    plot_survival()
    plot_ceac()
    # Tornado (Figure 3) is generated by gen_tornado.py from canonical
    # output/dsa_v2_full.json — plot_tornado() below is a stale 5-param
    # subset and must NOT be used (kept only for reference).
    # plot_tornado()
    plot_ce_plane()
    print("All figures complete.")