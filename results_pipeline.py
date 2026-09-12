#!/usr/bin/env python3
"""
results_pipeline.py - THE single canonical source for every reported quantity in
the intismeran CEA manuscript.

Design rules
------------
1. One model implementation only: ModelParams + CEAModelV2 (rerun_primary).
2. Nothing downstream is hand-typed: every reported number is computed here and
   exported as JSON + LaTeX fragments.
3. Algebraic invariants are asserted BEFORE export; any failure aborts the run:
     sum(cost components)            == model total cost        (per arm)
     sum(state QALYs) - AE disutility == model total QALY      (per arm)
     ICER                            == dCost / dQALY
     NMB(W)                          == W*dQALY - dCost        (every WTP)
     dCost(P)                        == dCost(P0) + (P - P0)   (affine in price)
     ICER(analytic threshold)        == WTP exactly
     P(CE at W)                      == mean(NMB(W) > 0)       (PSA)
4. Analytic threshold prices, not tolerance-based bisection outputs.

Usage
-----
  python results_pipeline.py                 # deterministic quantities
  python results_pipeline.py --psa           # + PSA (delegates to psa_v2.py; slow)
  python results_pipeline.py --check-only    # run checks, no export
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import subprocess
import sys
from dataclasses import dataclass, fields

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from model import (ModelParams, N_PEMBRO_DOSES, N_INTISMERAN_DOSES,   # noqa: E402
                   PEMBRO_CYCLE_WEEKS, discount_factors)
from rerun_primary import CEAModelV2               # noqa: E402

DEFAULT_PRICE = 200_000.0
WTP_GRID = (50_000, 100_000, 150_000, 200_000)
PRICE_GRID = list(range(0, 800_001, 50_000))
QDP = 3                     # QALY display precision; 2 dp is not row-additive

_FIELDS = {f.name for f in fields(ModelParams)}
# Attributes consumed by the model via getattr() but not declared on ModelParams
DYNAMIC_KEYS = {"no_post_trial_os_benefit", "wane_start_month", "wane_end_month"}
# Scenarios that are identical to the base case by construction
EXPECTED_EQUAL_TO_BASE = {"Base case", "GP floor 60 mo (base)", "40-year horizon (base)"}
FAILURES: list[str] = []
CASES: dict[str, "Case"] = {}


# ---------------------------------------------------------------- helpers
def build_params(**kw) -> ModelParams:
    kw.setdefault("constraint_general_pop", True)
    unknown = set(kw) - _FIELDS - DYNAMIC_KEYS
    if unknown:                     # a dropped kwarg silently yields the base case
        raise ValueError(f"unknown model parameter(s): {sorted(unknown)}")
    p = ModelParams(**{k: v for k, v in kw.items() if k in _FIELDS})
    # ModelParams.__post_init__ refits every survival parameter from the landmark data
    # and overwrites whatever the constructor was given, so re-apply every override
    # afterwards (same trick psa_v2.py uses). Without this, os_mu_* / os_weibull_*
    # perturbations are silently discarded and the DSA row is a no-op.
    for k, v in kw.items():
        setattr(p, k, v)
    return p


def quiet(fn, *a, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def check(label: str, got: float, want: float, tol: float = 1e-6) -> None:
    if not np.isfinite(got) or abs(got - want) > tol:
        FAILURES.append(f"{label}: got {got!r} want {want!r} (tol {tol:g})")


def money(v: float) -> str:
    return f"\\${v:,.0f}" if v >= 0 else "$-$" + f"{abs(v):,.0f}"


ZERO = "\\$0"


def bf(s: str) -> str:
    return f"\\textbf{{{s}}}"


def num(v: float, nd: int = 2) -> str:
    return f"{v:.{nd}f}"


# ---------------------------------------------------------------- one model run
@dataclass
class Case:
    label: str
    params: dict
    dcost: float
    dqaly: float
    dly: float
    dly_undisc: float
    arm: dict            # arm -> {cost, qaly, ly_disc, ly_undisc, cost_parts, qaly_parts}
    W: tuple = WTP_GRID

    @property
    def icer(self) -> float:
        return self.dcost / self.dqaly

    def nmb(self, w: float) -> float:
        return w * self.dqaly - self.dcost


def run_case(label: str, model_cls=CEAModelV2, **kw) -> Case:
    p = build_params(**kw)
    m = model_cls(p)
    res = quiet(m.run)
    sur = m._survival()
    states = m._state_proportions(sur)
    disc = discount_factors(p.discount_rate, m.t)
    cl = p.cycle_length

    arm = {}
    for a in ("combo", "pembro"):
        s = states[a]
        cost_parts = {
            "Intismeran": p.cost_intismeran * disc[0] if a == "combo" else 0.0,
            "Pembrolizumab": p.cost_keytruda_annual * disc[:12].sum() * cl,
            "Sequencing": p.cost_sequencing * disc[0] if a == "combo" else 0.0,
            "Administration": admin_cost_for(p, a),
            "Adverse events": p.cost_ae_incremental * disc[0] if a == "combo" else 0.0,
            "LR management": float(np.sum(s["lr"] * p.cost_lr_monthly * 12 * cl * disc)),
            "DM management": float(np.sum(s["dm"] * p.cost_dm_monthly * 12 * cl * disc)),
        }
        qaly_parts = {
            "RF": float(np.sum(s["rf"] * p.util_rf * disc) * cl),
            "LR": float(np.sum(s["lr"] * p.util_lr * disc) * cl),
            "DM": float(np.sum(s["dm"] * p.util_dm * disc) * cl),
            "AE disutility": -p.util_disutility_ae * disc[0] * cl,
        }
        # additivity: the exported table must reproduce the model's own totals
        check(f"{label}/{a} cost additivity", sum(cost_parts.values()), res[a]["cost"], 1e-6)
        check(f"{label}/{a} QALY additivity", sum(qaly_parts.values()), res[a]["qaly"], 1e-9)
        arm[a] = {
            "cost": res[a]["cost"],
            "qaly": res[a]["qaly"],
            "ly_disc": float(np.sum(sur[f"os_{a}"] * disc) * cl),
            "ly_undisc": float(np.sum(sur[f"os_{a}"]) * cl),
            "cost_parts": cost_parts,
            "qaly_parts": qaly_parts,
        }

    c, b = arm["combo"], arm["pembro"]
    return Case(
        label=label, params={k: v for k, v in kw.items()},
        dcost=c["cost"] - b["cost"], dqaly=c["qaly"] - b["qaly"],
        dly=c["ly_disc"] - b["ly_disc"], dly_undisc=c["ly_undisc"] - b["ly_undisc"],
        arm=arm,
    )


# ---------------------------------------------------------------- scenario registry
def admin_iv_cost(p, arm: str) -> float:
    """Pembrolizumab administration: 18 intravenous infusions (CPT 96365), both arms."""
    t = np.arange(N_PEMBRO_DOSES) * (PEMBRO_CYCLE_WEEKS / 52)
    return p.cost_admin_iv_infusion * float(discount_factors(p.discount_rate, t).sum())


def admin_im_cost(p, arm: str) -> float:
    """Intismeran administration: 9 intramuscular injections (CPT 96372), combination
    arm only -- intismeran is given intramuscularly, not as an intravenous infusion."""
    if arm != "combo":
        return 0.0
    t = np.arange(N_INTISMERAN_DOSES) * (PEMBRO_CYCLE_WEEKS / 52)
    return p.cost_admin_im_injection * float(discount_factors(p.discount_rate, t).sum())


def admin_cost_for(p, arm: str) -> float:
    """Total administration: IV infusions in every arm plus IM injections in combo."""
    return admin_iv_cost(p, arm) + admin_im_cost(p, arm)


def os_dists(kw: dict) -> tuple:
    """Effective (combo, pembro) OS distributions implied by a kwargs dict."""
    shared = kw.get("os_distribution", "lognormal")
    return (kw.get("os_distribution_combo", "") or shared,
            kw.get("os_distribution_pembro", "") or shared)


def rank1_params() -> dict:
    """Lowest total-AIC distribution combination from the 15,625-combination search."""
    cols = {"OS_combo": "os_distribution_combo", "OS_pembro": "os_distribution_pembro",
            "rDM_combo": "rDM_distribution_combo", "rDM_pembro": "rDM_distribution_pembro",
            "rLR_combo": "rLR_distribution_combo", "rLR_pembro": "rLR_distribution_pembro"}
    csv = os.path.join(HERE, "output", "survival_combos_15625.csv")
    with open(csv) as fh:
        header = fh.readline().strip().split(",")
        best = min(fh, key=lambda ln: float(ln.split(",")[header.index("total_AIC")]))
    row = dict(zip(header, best.strip().split(",")))
    return {cols[k]: row[k] for k in cols}


def waning_v2_cls():
    try:
        from waning_v2 import CEAModelV3
        return CEAModelV3
    except Exception as exc:                     # pragma: no cover
        print(f"  [warn] waning V2 unavailable ({exc}); skipped")
        return None


SCENARIOS = [
    ("Base case", {}, "GP-constrained, 40-year horizon, 3% discounting"),
    ("Weibull OS", dict(os_distribution="weibull"), "OS extrapolated with Weibull (both arms)"),
    ("Rank #1 joint model", "RANK1", "lowest total AIC of 15,625 distribution combinations"),
    ("Unconstrained (no GP)", dict(constraint_general_pop=False), "no general-population hazard floor"),
    ("GP floor 48 mo", dict(gp_floor_start_months=48), "GP hazard floor applied from month 48"),
    ("GP floor 60 mo (base)", dict(gp_floor_start_months=60), "GP hazard floor from month 60"),
    ("GP floor 72 mo", dict(gp_floor_start_months=72), "GP hazard floor applied from month 72"),
    ("GP floor 96 mo", dict(gp_floor_start_months=96), "GP hazard floor applied from month 96"),
    ("10-year horizon", dict(time_horizon_years=10), "costs and QALYs truncated at 10 years"),
    ("20-year horizon", dict(time_horizon_years=20), "truncated at 20 years"),
    ("30-year horizon", dict(time_horizon_years=30), "truncated at 30 years"),
    ("40-year horizon (base)", dict(time_horizon_years=40), "base-case horizon"),
    ("Waning V1 (OS-survival convergence)", dict(treatment_waning=True),
     "combo OS curve converges to pembro between 5 and 20 years"),
    ("Waning V2 (hazard convergence)", "V2", "monthly mortality hazard converges from month 60 to 120"),
    ("No post-trial OS benefit", dict(no_post_trial_os_benefit=True),
     "combo OS set equal to pembro OS after month 60"),
    ("0% discount", dict(discount_rate=0.0), "undiscounted"),
    ("5% discount", dict(discount_rate=0.05), "5% annual discounting"),
]

DSA_UTILS = [
    # label, attribute, base display, (low value, low display), (high value, high display)
    ("Utility RF", "util_rf", "0.83", (0.80, "0.80"), (0.86, "0.86")),
    ("Utility LR", "util_lr", "0.64", (0.61, "0.61"), (0.67, "0.67")),
    ("Utility DM", "util_dm", "0.55", (0.52, "0.52"), (0.58, "0.58")),
    ("Pembrolizumab annual cost", "cost_keytruda_annual", "\\$220,896",
     (176_717, "\\$176,717"), (265_075, "\\$265,075")),
    ("Intismeran cost", "cost_intismeran", "\\$200,000",
     (160_000, "\\$160,000"), (240_000, "\\$240,000")),
    ("LR monthly cost", "cost_lr_monthly", "\\$3,000", (2_400, "\\$2,400"),
     (3_600, "\\$3,600")),
    ("DM monthly cost", "cost_dm_monthly", "\\$12,000", (9_600, "\\$9,600"),
     (14_400, "\\$14,400")),
    ("Discount rate", "discount_rate", "3\\%", (0.0, "0\\%"), (0.05, "5\\%")),
    # Route-specific administration: intramuscular injection (intismeran) vs the
    # intravenous infusion rate (pembrolizumab). Range spans $0 to the cost level
    # that would apply if the injection were charged like an infusion.
    ("Intismeran IM injection cost", "cost_admin_im_injection", "\\$13.91",
     (0.0, "\\$0"), (120.0, "\\$120")),
]
# Applied to both arms, so the incremental cost (and hence the ICER) cannot move.
EXPECT_NO_CHANGE = {"Pembrolizumab annual cost"}


# ---------------------------------------------------------------- builders
def build_deterministic() -> dict:
    print("== base case ==")
    CASES["base"] = run_case("base")
    base = CASES["base"]

    # Administration must reflect the stated regimen and route: 18 pembrolizumab IV
    # infusions in both arms (CPT 96365) and 9 intismeran IM injections in the
    # combination arm (CPT 96372). The implied discounted dose counts have to sit
    # just under those numbers.
    _bp = build_params()
    for label, implied, n_doses, rate in (
            ("pembrolizumab IV infusions", admin_iv_cost(_bp, "combo"),
             N_PEMBRO_DOSES, _bp.cost_admin_iv_infusion),
            ("intismeran IM injections", admin_im_cost(_bp, "combo"),
             N_INTISMERAN_DOSES, _bp.cost_admin_im_injection)):
        n_implied = implied / rate
        if not n_doses - 1.0 <= n_implied <= n_doses:
            FAILURES.append(f"{label}: implied discounted dose count {n_implied:.2f} "
                            f"does not match the stated {n_doses}")

    print("== price ladder ==")
    price_cases = {}
    for price in sorted(set(PRICE_GRID + [100_000, 300_000, 500_000])):
        cs = run_case(f"price {price}", cost_intismeran=float(price))
        # affine check: dCost moves one-for-one with price, dQALY unchanged
        check(f"price {price} dQALY invariance", cs.dqaly, base.dqaly, 1e-9)
        check(f"price {price} dCost affinity", cs.dcost,
              base.dcost + (price - DEFAULT_PRICE), 1e-6)
        price_cases[price] = cs

    print("== analytic thresholds ==")
    thresholds = {}
    for w in (100_000, 150_000):
        p_star = DEFAULT_PRICE + (w * base.dqaly - base.dcost)
        at = run_case(f"threshold {w}", cost_intismeran=p_star)
        check(f"threshold {w} ICER", at.icer, w, 1e-6)
        thresholds[w] = p_star
        CASES[f"threshold_{w}"] = at

    print("== scenarios ==")
    scenario_rows = []
    v2_cls = None
    for label, kw, definition in SCENARIOS:
        if kw == "RANK1":
            case = run_case(label, **rank1_params())
        elif kw == "V2":
            if v2_cls is None:
                v2_cls = waning_v2_cls()
            if v2_cls is None:
                continue
            case = run_case(label, model_cls=v2_cls, treatment_waning=True,
                            wane_start_month=60, wane_end_month=120)
        else:
            case = run_case(label, **kw)
        base = CASES["base"]
        if (label not in EXPECTED_EQUAL_TO_BASE
                and abs(case.dcost - base.dcost) < 1e-6
                and abs(case.dqaly - base.dqaly) < 1e-9):
            FAILURES.append(f"scenario '{label}' is numerically identical to the base "
                            f"case (unimplemented switch or wrong parametrization)")
        # A partial override is invisible to the check above: changing only the OS
        # distribution used to leave dLY bit-identical to base while cost/QALY moved
        # (CEAModelV2 read only p.os_distribution and ignored the per-arm fields).
        if (isinstance(kw, dict) and os_dists(kw) != os_dists({})
                and abs(case.dly - base.dly) < 1e-9):
            FAILURES.append(f"scenario '{label}' overrides the OS distribution "
                            f"({os_dists(kw)}) but discounted LY is unchanged from base "
                            f"-> the override was ignored")
        scenario_rows.append(("-", label, definition, case))

    print("== DSA (one-way) ==")
    dsa_rows = []

    def _dsa_row(label: str, base_disp: str, lo_val, lo_disp: str,
                 hi_val, hi_disp: str, attr: str) -> None:
        rows = {"low_case": run_case(f"dsa {label} lo", **{attr: lo_val}),
                "high_case": run_case(f"dsa {label} hi", **{attr: hi_val})}
        if label not in EXPECT_NO_CHANGE:
            b = CASES["base"]
            if (abs(rows["low_case"].icer - b.icer) < 1
                    and abs(rows["high_case"].icer - b.icer) < 1):
                FAILURES.append(f"DSA row '{label}' does not move the ICER in either "
                                f"direction -> the perturbation was ignored")
        dsa_rows.append({"param": label, "attr": attr, "base": base_disp, "low": lo_disp,
                         "high": hi_disp, **rows})

    for label, attr, base_disp, (lo, lo_disp), (hi, hi_disp) in DSA_UTILS:
        _dsa_row(label, base_disp, lo, lo_disp, hi, hi_disp, attr)
    for which, lo_f, hi_f in (("combo", 0.6, 1.4), ("pembro", 0.8, 1.2)):
        mu = getattr(build_params(), f"os_mu_{which}")
        _dsa_row(f"OS $\\mu$ {which}", f"{mu:.3f}", mu * lo_f,
                 f"-{100 - 100 * lo_f:.0f}\\%", mu * hi_f, f"+{100 * hi_f - 100:.0f}\\%",
                 f"os_mu_{which}")

    return {"base": base, "prices": price_cases, "thresholds": thresholds,
            "scenarios": scenario_rows, "dsa": dsa_rows}


# ---------------------------------------------------------------- PSA
def build_psa(extra_args: list[str] | None = None) -> dict:
    print("== PSA (delegating to psa_v2.py) ==")
    cmd = [sys.executable, os.path.join(HERE, "psa_v2.py")] + list(extra_args or [])
    proc = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"psa_v2.py failed:\n{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")

    with open(os.path.join(HERE, "output", "psa_v2_results.json")) as fh:
        summary = json.load(fh)
    with open(os.path.join(HERE, "output", "psa_v2_draws.json")) as fh:
        draws = json.load(fh)

    dq = np.array([d["dq"] for d in draws])
    dc = np.array([d["dc"] for d in draws])
    for w in (50_000, 100_000, 150_000):
        reported = summary["pce"][f"${w:,}"]
        check(f"PSA P(CE) @ {w}", reported, float(np.mean(-dc + w * dq > 0)), 1e-9)
        check(f"PSA NMB @ {w}", summary["nmb"][f"${w:,}"]["mean"],
              float(np.mean(-dc + w * dq)), 1e-2)
    check("PSA n_valid", summary["n_valid"], len(draws), 0)
    return {"summary": summary, "n_draws": len(draws)}


# ---------------------------------------------------------------- LaTeX export
def tex_esc(s: str) -> str:
    """Escape characters that are special inside a LaTeX tabular cell."""
    return s.replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%") \
            .replace("#", "\\#").replace("_", "\\_")


def _rows_scenarios(rows) -> str:
    out = []
    for _, label, _definition, c in rows:
        out.append(f"{tex_esc(label)} & {num(c.dly)} & {num(c.dqaly)} & {money(c.dcost)} & "
                   f"{money(c.icer)} & {money(c.nmb(100_000))} & {money(c.nmb(150_000))} \\\\")
    return "\n".join(out)


def _rows_price(price_cases, thresholds) -> str:
    thr = {int(round(v)): w for w, v in thresholds.items()}
    other = {100_000: 150_000, 150_000: 100_000}
    out = []
    for price in sorted(set(list(price_cases) + list(thr))):
        if price in thr:
            w = thr[price]
            # by construction NMB(W) == 0 at the analytic threshold; report the other WTP too
            nmb_other = (other[w] - w) * price_cases[DEFAULT_PRICE].dqaly
            vals = {w: 0.0, other[w]: nmb_other}
            out.append(f"{bf(money(price))} & {bf(money(w))} & "
                       f"{bf(money(vals[100_000]))} & {bf(money(vals[150_000]))} \\\\")
            continue
        c = price_cases[price]
        ic = "Dominant" if c.icer < 0 else money(c.icer)
        lbl = f"{money(price)} (base)" if price == DEFAULT_PRICE else money(price)
        out.append(f"{lbl if price else ZERO} & {ic} & "
                   f"{money(c.nmb(100_000))} & {money(c.nmb(150_000))} \\\\")
    return "\n".join(out)


def export_tex(deterministic: dict, tex_dir: str) -> None:
    os.makedirs(tex_dir, exist_ok=True)
    base = deterministic["base"]
    c, b = base.arm["combo"], base.arm["pembro"]

    def w(name, body):
        path = os.path.join(tex_dir, name)
        with open(path, "w") as fh:
            fh.write(body.rstrip() + "\n")
        print(f"  wrote {path}")

    # Table 2: base-case cost + QALY decomposition
    lines = ["\\multicolumn{4}{c}{\\textit{Cost decomposition}} \\\\"]
    for comp in c["cost_parts"]:
        lines.append(f"{comp} & {money(c['cost_parts'][comp])} & {money(b['cost_parts'][comp])} "
                     f"& {money(c['cost_parts'][comp] - b['cost_parts'][comp])} \\\\")
    lines += ["\\midrule",
              f"Total cost & {money(c['cost'])} & {money(b['cost'])} & {money(base.dcost)} \\\\",
              "\\midrule",
              "\\multicolumn{4}{c}{\\textit{QALY decomposition}} \\\\"]
    for st in ("RF", "LR", "DM"):
        lines.append(f"{st} state & {num(c['qaly_parts'][st], QDP)} & {num(b['qaly_parts'][st], QDP)} "
                     f"& {num(c['qaly_parts'][st] - b['qaly_parts'][st], QDP)} \\\\")
    lines += ["\\midrule",
              f"Total QALYs & {num(c['qaly'], QDP)} & {num(b['qaly'], QDP)} & {num(base.dqaly, QDP)} \\\\"]
    w("tab_base_body.tex", "\n".join(lines))

    # S11: cost decomposition, arms + incremental (with an explicit total row)
    s11 = []
    for comp in c["cost_parts"]:
        s11.append(f"{comp} & {money(c['cost_parts'][comp])} & {money(b['cost_parts'][comp])} "
                   f"& {money(c['cost_parts'][comp] - b['cost_parts'][comp])} \\\\")
    s11 += ["\\midrule",
            f"{bf('Total')} & {bf(money(c['cost']))} & {bf(money(b['cost']))} "
            f"& {bf(money(base.dcost))} \\\\"]
    w("tab_cost_decomp_body.tex", "\n".join(s11))

    # S12: LY + QALY decomposition
    s12 = [f"Life-years (undiscounted) & {num(c['ly_undisc'])} & {num(b['ly_undisc'])} "
           f"& {num(base.dly_undisc)} \\\\",
           f"Life-years (discounted 3\\%) & {num(c['ly_disc'])} & {num(b['ly_disc'])} "
           f"& {num(base.dly)} \\\\",
           "\\midrule"]
    for st in ("RF", "LR", "DM"):
        s12.append(f"{st} state QALYs & {num(c['qaly_parts'][st], QDP)} & "
                   f"{num(b['qaly_parts'][st], QDP)} "
                   f"& {num(c['qaly_parts'][st] - b['qaly_parts'][st], QDP)} \\\\")
    s12 += ["\\midrule",
            f"\\textbf{{Total QALYs}} & \\textbf{{{num(c['qaly'], QDP)}}} & "
            f"\\textbf{{{num(b['qaly'], QDP)}}} & \\textbf{{{num(base.dqaly, QDP)}}} \\\\"]
    w("tab_ly_qaly_body.tex", "\n".join(s12))

    # Section 3.2 paragraph and the supplement's decomposition note, generated from
    # the same values as Table 2 so the prose cannot lag the table.
    main_par, supp_par = _decomp_prose(c, b, base)
    with open(os.path.join(HERE, "tables", "para_decomp.tex"), "w") as fh:
        fh.write(main_par + "\n")
    with open(os.path.join(HERE, "tables", "para_decomp_supp.tex"), "w") as fh:
        fh.write(supp_par + "\n")
    print("  wrote tables/para_decomp.tex, tables/para_decomp_supp.tex")

    # S13 / Table 4: scenarios
    w("tab_scenarios_body.tex", _rows_scenarios(deterministic["scenarios"]))

    # Table 4 (main text): base + pricing + alternative assumptions
    by_label = {label: cs for _, label, _d, cs in deterministic["scenarios"]}
    pr = deterministic["prices"]
    main = [f"Base case (GP-constrained) & {money(base.dcost)} & {num(base.dqaly)} & "
            f"{money(base.icer)} & {money(base.nmb(150_000))} \\\\",
            "\\addlinespace",
            "\\multicolumn{5}{c}{\\textit{Pricing}} \\\\"]
    for price in (100_000, 300_000, 500_000):
        pc = pr[price]
        main.append(f"Price: {money(price)}/course & {money(pc.dcost)} & {num(pc.dqaly)} & "
                    f"{money(pc.icer)} & {money(pc.nmb(150_000))} \\\\")
    main += ["\\addlinespace",
             "\\multicolumn{5}{c}{\\textit{Alternative assumptions}} \\\\"]
    for label, disp in (("Weibull OS", "Weibull OS"),
                        ("0% discount", "Discount rate: 0\\%"),
                        ("5% discount", "Discount rate: 5\\%"),
                        ("10-year horizon", "Time horizon: 10 years"),
                        ("20-year horizon", "Time horizon: 20 years")):
        cs = by_label[label]
        main.append(f"{disp} & {money(cs.dcost)} & {num(cs.dqaly)} & "
                    f"{money(cs.icer)} & {money(cs.nmb(150_000))} \\\\")
    main += ["\\addlinespace",
             "\\multicolumn{5}{c}{\\textit{Treatment-effect persistence}} \\\\"]
    for label, disp in (("Waning V1 (OS-survival convergence)",
                         "OS-survival convergence (5\\textendash20y)"),
                        ("Waning V2 (hazard convergence)",
                         "Hazard convergence (5\\textendash10y)"),
                        ("No post-trial OS benefit", "No post-trial OS benefit")):
        cs = by_label[label]
        ic = "Dominant" if cs.icer < 0 else money(cs.icer)
        dc = f"\\textendash{money(abs(cs.dcost))}" if cs.dcost < 0 else money(cs.dcost)
        main.append(f"{disp} & {dc} & {num(cs.dqaly)} & {ic} & "
                    f"{money(cs.nmb(150_000))} \\\\")
    w("tab_main_scenarios_body.tex", "\n".join(main))

    # S16: price ladder + analytic thresholds (no bisection outputs)
    w("tab_price_body.tex", _rows_price(deterministic["prices"], deterministic["thresholds"]))

    # DSA (matches the supplement column layout: Parameter & Base & Low & High &
    # ICER range & Change in NMB)
    def _ic(c: Case) -> str:
        return "Dominant" if c.icer < 0 else f"{c.icer:,.0f}"

    dsa = []
    for row in deterministic["dsa"]:
        lo, hi = row["low_case"], row["high_case"]
        dnmb = abs(hi.nmb(150_000) - lo.nmb(150_000))
        dsa.append(f"{row['param']} & {row['base']} & {row['low']} & {row['high']} "
                   f"& {_ic(lo)}--{_ic(hi)} & {money(dnmb)} \\\\")
    w("tab_dsa_body.tex", "\n".join(dsa))


TEX_SLOTS = {
    "para_decomp": ["manuscript.tex"],
    "para_decomp_supp": ["supplementary.tex", "supplementary_blind.tex"],
    "tab_base_body": ["manuscript.tex"],
    "tab_main_scenarios_body": ["manuscript.tex"],
    "tab_cost_decomp_body": ["supplementary.tex", "supplementary_blind.tex"],
    "tab_ly_qaly_body": ["supplementary.tex", "supplementary_blind.tex"],
    "tab_scenarios_body": ["supplementary.tex", "supplementary_blind.tex"],
    "tab_dsa_body": ["supplementary.tex", "supplementary_blind.tex"],
    "tab_price_body": ["supplementary.tex", "supplementary_blind.tex"],
}


def sync_tex(frag_dir: str = "tables") -> None:
    """Inline each generated fragment between markers in the .tex sources.

    \\input{} cannot be used inside a tabular body: booktabs' \\midrule issues
    \\noalign and the following \\omit (first token of \\multicolumn) is then
    misplaced. Markers keep the pipeline as the single source without that trap.
    """
    for name, files in TEX_SLOTS.items():
        with open(os.path.join(HERE, frag_dir, f"{name}.tex")) as fh:
            body = fh.read().rstrip("\n").splitlines()
        if not body:
            raise SystemExit(f"empty fragment: {name}")
        open_m, close_m = f"% >>> pipeline {name}", f"% <<< pipeline {name}"
        legacy_in = f"\\input{{{frag_dir}/{name}.tex}}"
        for f in files:
            path = os.path.join(HERE, f)
            with open(path) as fh:
                lines = fh.read().splitlines()
            starts = [i for i, l in enumerate(lines) if l.strip() == open_m]
            if starts:
                i0 = starts[0]
                ends = [i for i, l in enumerate(lines) if l.strip() == close_m and i > i0]
                if not ends:
                    raise SystemExit(f"unclosed marker {open_m} in {f}")
                lines[i0 + 1:ends[0]] = body
            else:
                hit = [i for i, l in enumerate(lines) if l.strip() == legacy_in]
                if len(hit) != 1:
                    raise SystemExit(f"slot {name}: {len(hit)} anchors in {f} (expected 1)")
                lines[hit[0]:hit[0] + 1] = [open_m] + body + [close_m]
            with open(path, "w") as fh:
                fh.write("\n".join(lines) + "\n")
        print(f"  synced {name} -> {', '.join(files)}")


def export_dsa_json(deterministic: dict, out: str) -> None:
    """DSV data for the tornado figure (Figure 3) -- same numbers as the DSA table."""
    base = deterministic["base"]
    rows = []
    for r in deterministic["dsa"]:
        lo, hi = r["low_case"].icer, r["high_case"].icer
        rows.append({"param": r["attr"],
                     "lo_icer": "Dominant" if lo <= 0 else f"{lo:,.0f}",
                     "hi_icer": "Dominant" if hi <= 0 else f"{hi:,.0f}",
                     "lo_icer_raw": None if lo <= 0 else lo,
                     "hi_icer_raw": None if hi <= 0 else hi})
    with open(out, "w") as fh:
        json.dump({"base_icer": base.icer, "rows": rows}, fh, indent=1)
    print(f"  wrote {out}")


def _decomp_prose(c: dict, b: dict, base: Case) -> tuple[str, str]:
    """Prose quoting the Table 2 decomposition, built from the same values."""
    rf = c["qaly_parts"]["RF"] - b["qaly_parts"]["RF"]
    lr = c["qaly_parts"]["LR"] - b["qaly_parts"]["LR"]
    dm = c["qaly_parts"]["DM"] - b["qaly_parts"]["DM"]
    lrc, lrp = c["cost_parts"]["LR management"], b["cost_parts"]["LR management"]
    dmc, dmp = c["cost_parts"]["DM management"], b["cost_parts"]["DM management"]
    main = (
        f"The {money(c['cost_parts']['Intismeran'])} acquisition cost of intismeran contributed "
        f"substantially to the incremental cost. Under the corrected general-population mortality "
        f"inputs, the combination arm had {money(lrc - lrp)} higher locoregional-recurrence management "
        f"costs ({money(lrc)} vs {money(lrp)}) but {money(abs(dmc - dmp))} lower distant-metastasis "
        f"management costs ({money(dmc)} vs {money(dmp)}) than the \\mbox{{pembrolizumab}}-alone arm, "
        f"because \\mbox{{pembrolizumab}}-alone patients spend more time in the DM state over the "
        f"lifetime ({num(b['qaly_parts']['DM'], QDP)} vs {num(c['qaly_parts']['DM'], QDP)} discounted DM "
        f"QALYs). The incremental QALY gain of {num(base.dqaly, QDP)} was driven primarily by additional "
        f"time spent in the recurrence-free state ({num(rf, QDP)} RF QALYs gained), with a modest "
        f"additional {num(lr, QDP)} LR QALYs, partly offset by a {num(abs(dm), QDP)}-QALY reduction in "
        f"the DM state."
    )
    supp = (
        f"Health benefit is driven by RF state QALYs ({num(rf, QDP)} gained), with a modest additional "
        f"contribution from LR ({num(lr, QDP)}) and a small reduction in DM state QALYs "
        f"({('$-$' + num(abs(dm), QDP)) if dm < 0 else num(dm, QDP)}) "
        f"under the general-population mortality constraint."
    )
    return main, supp


def export_json(deterministic: dict, psa: dict | None, out: str) -> None:
    """JSON companion: every reported quantity at full precision."""
    def ser(c: Case) -> dict:
        return {"label": c.label, "params": c.params, "dcost": c.dcost, "dqaly": c.dqaly,
                "dly": c.dly, "dly_undisc": c.dly_undisc, "icer": c.icer,
                "nmb": {f"${w:,}": c.nmb(w) for w in WTP_GRID},
                "arm": c.arm}

    doc = {
        "base": ser(deterministic["base"]),
        "scenarios": [{"label": l, "definition": d, "result": ser(c)}
                      for _, l, d, c in deterministic["scenarios"]],
        "prices": {str(k): ser(v) for k, v in deterministic["prices"].items()},
        "thresholds_analytic": {str(k): v for k, v in deterministic["thresholds"].items()},
    }
    if psa:
        doc["psa"] = {"summary": psa["summary"], "n_draws": psa["n_draws"]}
    with open(out, "w") as fh:
        json.dump(doc, fh, indent=2)
    print(f"  wrote {out}")


def export_prose(deterministic: dict, out: str) -> None:
    base = deterministic["base"]
    c, b = base.arm["combo"], base.arm["pembro"]
    L = [
        "# Canonical prose values (do not hand-type; regenerate with results_pipeline.py)",
        "",
        f"- Combo: {c['cost']:,.0f} cost / {c['qaly']:.2f} QALY / "
        f"{c['ly_disc']:.2f} discounted LY",
        f"- Pembro: {b['cost']:,.0f} cost / {b['qaly']:.2f} QALY / "
        f"{b['ly_disc']:.2f} discounted LY",
        f"- Incremental: {base.dcost:,.0f} cost / {base.dqaly:.2f} QALY / "
        f"{base.dly:.2f} discounted LY / {base.dly_undisc:.2f} undiscounted LY",
        f"- ICER: {base.icer:,.0f} per QALY",
        f"- NMB @100K: {base.nmb(100_000):,.0f}   @150K: {base.nmb(150_000):,.0f}",
        f"- RF/LR/DM QALY (combo): {c['qaly_parts']['RF']:.2f} / {c['qaly_parts']['LR']:.2f} "
        f"/ {c['qaly_parts']['DM']:.2f}",
        f"- Incremental RF/LR/DM QALY: "
        f"{c['qaly_parts']['RF']-b['qaly_parts']['RF']:.2f} / "
        f"{c['qaly_parts']['LR']-b['qaly_parts']['LR']:.2f} / "
        f"{c['qaly_parts']['DM']-b['qaly_parts']['DM']:.2f}",
        f"- Administration cost per arm: {c['cost_parts']['Administration']:,.0f}",
        "",
        "## Analytic threshold prices",
    ]
    for wq, price in sorted(deterministic["thresholds"].items()):
        L.append(f"- {wq:,}/QALY -> {price:,.2f} per course")
    L += ["", "## Scenario ICERs"]
    for _, label, _d, cs in deterministic["scenarios"]:
        L.append(f"- {label}: {cs.icer:,.0f} (dCost {cs.dcost:,.0f}, dQALY {cs.dqaly:.2f})")
    with open(out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"  wrote {out}")


# ---------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--psa", action="store_true", help="also run PSA (slow)")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "output"))
    ap.add_argument("--tex-dir", default=os.path.join(HERE, "tables"))
    args = ap.parse_args()

    det = build_deterministic()
    psa = build_psa() if args.psa else None

    if FAILURES:
        print(f"\nFAILED {len(FAILURES)} invariant(s):")
        for f in FAILURES:
            print(f"  !! {f}")
        return 1
    print("\nAll invariants hold.")

    if not args.check_only:
        export_json(det, psa, os.path.join(args.out_dir, "canonical_results.json"))
        export_prose(det, os.path.join(args.out_dir, "canonical_values.md"))
        export_tex(det, args.tex_dir)
        sync_tex(args.tex_dir)
        export_dsa_json(det, os.path.join(args.out_dir, "dsa_canonical.json"))

    base = det["base"]
    print(f"\nbase: dCost ${base.dcost:,.0f}  dQALY {base.dqaly:.2f}  "
          f"ICER ${base.icer:,.0f}  thresholds "
          + ", ".join(f"${p:,.0f}@{w:,}" for w, p in sorted(det['thresholds'].items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
