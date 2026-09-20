"""Guard: the model's monthly survival arrays must be indexed by month.

Regression for the cycle-indexing drift that made Supplementary Table S6 read
S(month-1), the long-term table S5 read S(month+1) and the text read S(month).

Run:  python check_survival_indexing.py
"""
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import ModelParams, lognorm_surv, survival_at_month  # noqa: E402
from rerun_primary import CEAModelV2  # noqa: E402


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
