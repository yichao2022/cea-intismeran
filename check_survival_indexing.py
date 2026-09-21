"""Guard: the model's monthly survival arrays must be indexed by month, and every
survival value printed in the supplementary tables must reproduce from them.

Regression for the cycle-indexing drift that made Supplementary Table S6 read
S(month-1), the long-term table S5 read S(month+1) and the text read S(month).

Covers: index == month for the curves read through survival_at_month(); the
published 60-month text values; the observed-vs-fitted landmark table (tab:s15),
the external plausibility table (tab:s7) and the long-term table (tab:s5) -- the
last two are hand-maintained in supplementary.tex, so this is their only guard.

Run:  python check_survival_indexing.py
"""
import hashlib
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import (ModelParams, general_pop_surv, lognorm_surv,  # noqa: E402
                   survival_at_month)
from rerun_primary import CEAModelV2  # noqa: E402

TEX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "supplementary.tex")
CURVE = {"OS": "os", "DMFS": "dmfs", "RFS": "rfs"}
ARM = {"Combo": "combo", "Pembro": "pembro"}


def table_rows(label: str) -> list:
    """Lines between \\label{<label>} and that table's \\end{table}."""
    rows, inside = [], False
    for line in open(TEX):
        if f"\\label{{{label}}}" in line:
            inside = True
            continue
        if inside:
            if line.startswith("\\end{table}"):
                break
            rows.append(line)
    assert rows, f"table {label} not found in {os.path.basename(TEX)}"
    return rows


def pct(tok: str):
    """'74.5\\%' -> 74.5; '---' / 'NR' -> None."""
    tok = tok.strip().replace("\\%", "").replace("$", "").replace("{", "").replace("}", "")
    return float(tok) if re.fullmatch(r"\d+(\.\d+)?", tok) else None


def check_tables(sur: dict) -> None:
    """Every printed survival value must come back out of the model."""
    fails = []
    checked = 0

    for line in table_rows("tab:s15"):  # observed vs fitted landmark table
        m = re.match(r"\s*(OS|DMFS|RFS)\s*&\s*(Combo|Pembro)\s*&\s*(\d+)\s*&([^&]+)&([^&]+)&", line)
        if not m:
            continue
        endpoint, arm, month, fitted = m.group(1), m.group(2), int(m.group(3)), pct(m.group(5))
        if fitted is None:
            continue
        want = 100 * survival_at_month(sur[f"{CURVE[endpoint]}_{ARM[arm]}"], month)
        checked += 1
        if abs(want - fitted) > 0.06:
            fails.append(f"tab:s15 {endpoint} {arm} {month}m: table {fitted}% vs model {want:.2f}%")

    for line in table_rows("tab:s7"):  # external plausibility, model column
        m = re.match(r"\s*(OS|DMFS|RFS)\s*&\s*(\d+)\s*&([^&]+)&", line)
        if not m:
            continue
        endpoint, year, model_val = m.group(1), int(m.group(2)), pct(m.group(3))
        if model_val is None:
            continue
        want = 100 * survival_at_month(sur[f"{CURVE[endpoint]}_pembro"], 12 * year)
        checked += 1
        if abs(want - model_val) > 0.06:
            fails.append(f"tab:s7 {endpoint} {year}y: table {model_val}% vs model {want:.2f}%")

    for line in table_rows("tab:s5"):  # long-term OS vs general population
        m = re.match(r"\s*(\d+)\s*&\s*\d+\s*&([^&]+)&([^&]+)&([^&]+)\\\\", line)
        if not m:
            continue
        year, combo, pembro, gp = int(m.group(1)), pct(m.group(2)), pct(m.group(3)), pct(m.group(4))
        for arm, printed in (("combo", combo), ("pembro", pembro)):
            if printed is None:
                continue
            want = 100 * survival_at_month(sur[f"os_{arm}"], 12 * year)
            checked += 1
            if abs(want - printed) > 0.06:
                fails.append(f"tab:s5 {year}y os_{arm}: table {printed}% vs model {want:.2f}%")
        if gp is not None:
            want_gp = 100 * float(general_pop_surv(np.array([float(year)]))[0])
            checked += 1
            if abs(want_gp - gp) > 0.06:
                fails.append(f"tab:s5 {year}y GP: table {gp}% vs life table {want_gp:.2f}%")

    assert not fails, "printed values do not reproduce from the model:\n  " + "\n  ".join(fails)
    print(f"  {checked} printed table values reproduce from the model  OK")


def analytic(month: float, mu: float, sigma: float) -> float:
    """Fitted log-normal S(month) with no GP floor -- the reference for months <= 60."""
    return float(lognorm_surv(np.array([float(month)]), mu, sigma)[0])


def main() -> int:
    p = ModelParams(constraint_general_pop=True)
    m = CEAModelV2(p)
    sur = m._survival()

    print("month | model os_pembro | analytic log-normal | diff (pp)")
    worst = 0.0
    for month in (17, 18, 23, 24, 35, 36, 47, 48, 59, 60):
        got = 100 * survival_at_month(sur["os_pembro"], month)
        want = 100 * analytic(month, p.os_mu_pembro, p.os_sigma_pembro)
        diff = got - want
        worst = max(worst, abs(diff))
        print(f"{month:>5} | {got:>15.2f} | {want:>19.2f} | {diff:>+8.3f}")

    # Landmarks at or before the GP floor start (60 months) must reproduce the
    # parametric curve exactly: no floor binds there.
    assert worst < 0.05, f"model does not reproduce S(month); worst diff {worst:.3f} pp"

    # The published text values sit on exact months.
    checks = {
        ("os_pembro", 60): 74.5,
        ("rfs_pembro", 60): 40.5,
        ("dmfs_pembro", 60): 53.8,
    }
    for (curve, month), published in checks.items():
        got = 100 * survival_at_month(sur[curve], month)
        assert abs(got - published) <= 0.15, (
            f"{curve}[{month}] = {got:.2f}%, published {published}%")
        print(f"  {curve}@{month}m = {got:.2f}%  (text {published}%)  OK")

    check_tables(sur)

    # The committed canonical results must belong to the survival table on disk.
    root = os.path.dirname(os.path.abspath(__file__))
    doc_path = os.path.join(root, "output", "canonical_results.json")
    if os.path.exists(doc_path):
        lt = os.path.join(root, "data", "life_table_surv.json")
        digest = hashlib.sha256(open(lt, "rb").read()).hexdigest()[:16]
        stamped = json.load(open(doc_path)).get("inputs", {}).get("life_table_sha256")
        assert stamped == digest, (
            f"canonical results were built from life table {stamped}, current table is "
            f"{digest} -> rerun results_pipeline.py")
        print(f"  life table fingerprint {digest} matches canonical results  OK")

    print("\nPASS: index == month for every curve read through survival_at_month()")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
