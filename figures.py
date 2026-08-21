"""
cea_intismeran/figures.py
Generate publication-quality figures for the CEA manuscript.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
from model import ModelParams, PartitionedSurvivalModel, run_psa

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
    model = PartitionedSurvivalModel(params)
    t = model.t
    r = model.run()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    for arm_name in ["combo", "pembro"]:
        d = r[arm_name]
        label = "Intismeran + Pembro" if arm_name == "combo" else "Pembrolizumab"
        ax1.plot(t, d["pfs"], color=C[arm_name], label=label)
        ax2.plot(t, d["os"], color=C[arm_name], label=label)

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
    psa = run_psa(1000)
    thresholds = np.arange(0, 400_001, 10_000)
    prob_ce = [np.mean(psa["icers"] <= t) for t in thresholds]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(thresholds / 1000, prob_ce, color="#2166ac", lw=2)
    ax.axvline(100, color="gray", ls="--", alpha=0.5, label="$100K/QALY")
    ax.axvline(150, color="gray", ls=":", alpha=0.5, label="$150K/QALY")
    ax.set_xlabel("WTP Threshold ($1000/QALY)"); ax.set_ylabel("P(Cost-Effective)")
    ax.set_title("Cost-Effectiveness Acceptability Curve"); ax.legend()
    ax.set_ylim(0, 1); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/fig2_ceac.pdf"); plt.savefig("output/fig2_ceac.png")
    plt.close(); print("OK fig2")


# ── 3. Tornado ──
def plot_tornado():
    params = ModelParams()
    base = PartitionedSurvivalModel(params).run()
    base_icer = base["icer"]

    params_to_test = [
        ("Intismeran Price", "cost_intismeran_course", 100_000, 300_000),
        ("Keytruda Annual", "cost_keytruda_year", 150_000, 300_000),
        ("Utility DFS", "util_dfs", 0.75, 0.90),
        ("Utility DM", "util_dm", 0.55, 0.75),
        ("Metastatic Cost/Mo", "cost_dm_month", 8_000, 18_000),
        ("5yr RFS Combo", "rfs_5yr_combo", 0.60, 0.78),
        ("5yr RFS Pembro", "rfs_5yr_pembro", 0.40, 0.58),
        ("Discount Rate", "discount_rate", 0.0, 0.05),
    ]

    results = []
    for label, attr, low, high in params_to_test:
        for val in (low, high):
            # Pass via constructor so __post_init__ recomputes derived scales
            p = ModelParams(**{attr: val})
            icer = PartitionedSurvivalModel(p).run()["icer"]
            results.append((label, val, icer))

    grouped = defaultdict(list)
    for label, val, icer in results:
        grouped[label].append((val, icer))

    sorted_params = sorted(grouped.items(),
                           key=lambda x: abs(x[1][0][1] - x[1][1][1]),
                           reverse=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(len(sorted_params))

    for i, (label, values) in enumerate(sorted_params):
        low_val = min(v[1] for v in values)
        high_val = max(v[1] for v in values)
        low_param = [v[0] for v in values if v[1] == low_val][0]
        high_param = [v[0] for v in values if v[1] == high_val][0]

        left = base_icer - (base_icer - low_val)
        width = high_val - low_val
        color = "#2166ac" if low_val < base_icer else "#b2182b"
        ax.barh(i, width, left=left, height=0.6, color=color, alpha=0.8)
        ax.text(low_val - 4000, i, f"${low_val:,.0f}", va="center", ha="right", fontsize=8)
        ax.text(high_val + 4000, i, f"${high_val:,.0f}", va="center", ha="left", fontsize=8)

    ax.axvline(base_icer, color="black", lw=1.5,
               label=f"Base: ${base_icer:,.0f}")
    ax.set_yticks(y_pos)
    ax.set_yticklabels([p[0] for p in sorted_params])
    ax.set_xlabel("ICER ($/QALY)")
    ax.set_title("Tornado Diagram: One-Way Sensitivity Analysis")
    ax.legend(fontsize=9); ax.grid(alpha=0.3, axis="x")
    plt.tight_layout()
    plt.savefig("output/fig3_tornado.pdf"); plt.savefig("output/fig3_tornado.png")
    plt.close(); print("OK fig3")


# ── 4. CE Plane ──
def plot_ce_plane():
    psa = run_psa(1000)
    dc = np.array(psa["costs"]["combo"]) - np.array(psa["costs"]["pembro"])
    dq = np.array(psa["qalys"]["combo"]) - np.array(psa["qalys"]["pembro"])

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(dq, dc / 1e6, alpha=0.3, s=10, color="#2166ac")

    wtp_line = np.linspace(0, max(dq), 100)
    for thr, label, style in [(100_000, "$100K/QALY", "--"),
                               (150_000, "$150K/QALY", ":")]:
        ax.plot(wtp_line, thr * wtp_line / 1e6, color="gray",
                ls=style, alpha=0.5, label=label)

    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Incremental QALYs"); ax.set_ylabel("Incremental Cost ($M)")
    ax.set_title("Cost-Effectiveness Plane (1,000 PSA Iterations)")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/fig4_ce_plane.pdf"); plt.savefig("output/fig4_ce_plane.png")
    plt.close(); print("OK fig4")


if __name__ == "__main__":
    print("Generating figures...")
    plot_survival(); plot_ceac(); plot_tornado(); plot_ce_plane()
    print("All figures saved to output/")