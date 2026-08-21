"""
cea_intismeran/model.py
Partitioned Survival Model for Intismeran + Keytruda vs Keytruda alone
in adjuvant melanoma (Stage IIB-IV, completely resected).

Reference: KEYNOTE-942/mRNA-4157-P201 5-year follow-up (2026 ASCO)
          INTerpath-001 Phase 3 top-line (Merck 2026-08-19)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Callable, Tuple


# ── Survival distributions ──────────────────────────────────────────────────

def weibull_surv(t: np.ndarray, scale: float, shape: float) -> np.ndarray:
    """Weibull survival S(t) = exp(-(t/scale)^shape)"""
    return np.exp(-((t / scale) ** shape))


def fit_weibull_from_survrate(
    t_known: float, s_known: float, shape: float = 1.2
) -> Tuple[float, float]:
    """
    Fit Weibull assuming fixed shape (typical for cancer RFS) and
    solving for scale from S(t_known) = s_known.
    """
    scale = t_known / ((-np.log(s_known)) ** (1 / shape))
    return scale, shape


# ── Parameter container ─────────────────────────────────────────────────────

@dataclass
class ModelParams:
    """All model inputs. Fill from literature as data becomes available."""
    # ── Clinical (KEYNOTE-942 / INTerpath-001) ──
    rfs_5yr_combo: float = 0.688       # 5-year RFS, combo arm
    rfs_5yr_pembro: float = 0.491      # 5-year RFS, pembro alone
    rfs_hr: float = 0.51               # HR for RFS (95% CI 0.294-0.887)
    dmfs_hr: float = 0.411             # HR for DMFS (95% CI 0.200-0.843)
    os_5yr_combo: float = 0.922        # 5-year OS, combo (exploratory)
    os_5yr_pembro: float = 0.713       # 5-year OS, pembro alone

    # ── Costs (USD, 2026) — PLACEHOLDER, replace with literature values ──
    cost_intismeran_course: float = 200_000      # Jefferies estimate
    cost_keytruda_year: float = 190_000          # WAC, approximate
    cost_admin_per_cycle: float = 500            # infusion/admin cost
    cost_ae_combo: float = 15_000               # incremental AE cost, combo arm
    cost_dm_month: float = 12_000               # metastatic disease, monthly
    cost_sequencing: float = 5_000               # tumor sequencing + manufacturing

    # ── Utilities (EQ-5D) — PLACEHOLDER, replace with literature ──
    util_dfs: float = 0.85
    util_dm: float = 0.65
    util_disutility_ae: float = 0.05            # one-time disutility for AEs

    # ── Model settings ──
    cycle_length: float = 1 / 12                # 1 month
    time_horizon: float = 40.0                  # years
    discount_rate: float = 0.03                 # 3%/year

    # ── Derived survival parameters ──
    # (computed in __post_init__)
    rfs_scale_combo: float = field(init=False)
    rfs_scale_pembro: float = field(init=False)
    os_scale_combo: float = field(init=False)
    os_scale_pembro: float = field(init=False)
    weibull_shape: float = 1.2

    def __post_init__(self):
        # Fit Weibull scales from 5-year survival rates
        self.rfs_scale_combo, _ = fit_weibull_from_survrate(
            5.0, self.rfs_5yr_combo, self.weibull_shape
        )
        self.rfs_scale_pembro, _ = fit_weibull_from_survrate(
            5.0, self.rfs_5yr_pembro, self.weibull_shape
        )
        self.os_scale_combo, _ = fit_weibull_from_survrate(
            5.0, self.os_5yr_combo, self.weibull_shape
        )
        self.os_scale_pembro, _ = fit_weibull_from_survrate(
            5.0, self.os_5yr_pembro, self.weibull_shape
        )


# ── PSM Core ─────────────────────────────────────────────────────────────────

class PartitionedSurvivalModel:
    """
    Three-state partitioned survival model:
      DFS (disease-free survival) → DM (distant metastasis) → Death
    """

    def __init__(self, params: ModelParams):
        self.p = params
        n_cycles = int(np.ceil(params.time_horizon / params.cycle_length))
        self.t = np.arange(1, n_cycles + 1) * params.cycle_length  # years
        self.n_cycles = n_cycles

    def _surv_curves(self, arm: str) -> Tuple[np.ndarray, np.ndarray]:
        """Return PFS(t) and OS(t) for given arm."""
        if arm == "combo":
            pfs = weibull_surv(self.t, self.p.rfs_scale_combo, self.p.weibull_shape)
            os = weibull_surv(self.t, self.p.os_scale_combo, self.p.weibull_shape)
        else:
            pfs = weibull_surv(self.t, self.p.rfs_scale_pembro, self.p.weibull_shape)
            os = weibull_surv(self.t, self.p.os_scale_pembro, self.p.weibull_shape)
        return pfs, os

    def _state_proportions(self, pfs: np.ndarray, os: np.ndarray) -> dict:
        """Return dict of state proportion vectors."""
        return {
            "dfs": pfs,
            "dm": np.maximum(0, os - pfs),  # alive but progressed
            "death": 1 - os,
        }

    def run_arm(self, arm: str) -> dict:
        """Run one arm and return QALYs, costs, and state proportions."""
        pfs, os = self._surv_curves(arm)
        states = self._state_proportions(pfs, os)

        # Discount factor (annual discount, applied at cycle midpoint)
        disc = np.exp(-self.p.discount_rate * (self.t - self.p.cycle_length / 2))

        # ── QALYs ──
        qaly = (states["dfs"] * self.p.util_dfs + states["dm"] * self.p.util_dm) * disc
        total_qaly = np.sum(qaly) * self.p.cycle_length

        # ── Costs ──
        cycle_cost = np.zeros(self.n_cycles)

        if arm == "combo":
            # Intismeran 1mg q3w × 9 doses (~6 months) + Keytruda 400mg q6w × 9 (~1 year)
            # Simplification: drug cost spread over 12 months
            treatment_months = 12
            drug_cost = (self.p.cost_intismeran_course + self.p.cost_keytruda_year)
            monthly_tx = drug_cost / treatment_months
            # Treatment cost only in the first year
            first_year_mask = self.t <= 1.0
            cycle_cost[first_year_mask] = monthly_tx * self.p.cycle_length * 12

            # Sequencing cost (one-time)
            cycle_cost[0] += self.p.cost_sequencing

            # AE cost (one-time, first cycle)
            cycle_cost[0] += self.p.cost_ae_combo

        else:
            # Keytruda alone, ~1 year
            treatment_months = 12
            monthly_tx = self.p.cost_keytruda_year / treatment_months
            first_year_mask = self.t <= 1.0
            cycle_cost[first_year_mask] = monthly_tx * self.p.cycle_length * 12

        # Admin cost per cycle (both arms)
        cycle_cost += self.p.cost_admin_per_cycle

        # Metastatic disease cost (only patients in DM state pay this)
        dm_cost = states["dm"] * self.p.cost_dm_month * 12 * self.p.cycle_length

        total_cost = np.sum((cycle_cost + dm_cost) * disc)

        return {
            "arm": arm,
            "total_qaly": total_qaly,
            "total_cost": total_cost,
            "states": states,
            "pfs": pfs,
            "os": os,
        }

    def run(self) -> dict:
        """Run both arms and return ICER."""
        combo = self.run_arm("combo")
        pembro = self.run_arm("pembro")

        delta_cost = combo["total_cost"] - pembro["total_cost"]
        delta_qaly = combo["total_qaly"] - pembro["total_qaly"]
        icer = delta_cost / delta_qaly if delta_qaly > 0 else np.inf

        return {
            "combo": combo,
            "pembro": pembro,
            "delta_cost": delta_cost,
            "delta_qaly": delta_qaly,
            "icer": icer,
        }


# ── Base case runner ─────────────────────────────────────────────────────────

def base_case() -> dict:
    """Run base case with default parameters."""
    params = ModelParams()
    model = PartitionedSurvivalModel(params)
    result = model.run()
    return result


def print_results(result: dict):
    """Pretty-print base case results."""
    r = result
    print(f"{'='*60}")
    print(f"  CEA: Intismeran + Keytruda vs Keytruda (Adjuvant Melanoma)")
    print(f"  Reference: KEYNOTE-942 5-year follow-up (2026 ASCO)")
    print(f"{'='*60}")
    for arm_name in ["pembro", "combo"]:
        arm = r[arm_name]
        print(f"\n  {arm['arm'].upper()}:")
        print(f"    Total QALYs:  {arm['total_qaly']:>8.2f}")
        print(f"    Total Cost:  ${arm['total_cost']:>10,.0f}")

    print(f"\n  {'─'*50}")
    print(f"  Incremental Cost:  ${r['delta_cost']:>10,.0f}")
    print(f"  Incremental QALYs: {r['delta_qaly']:>10.2f}")
    print(f"  ICER:             ${r['icer']:>10,.0f} / QALY")
    print(f"  {'─'*50}")
    print(f"\n  ⚠  COSTS AND UTILITIES ARE PLACEHOLDERS.")
    print(f"     Replace with literature values before publication.\n")


if __name__ == "__main__":
    result = base_case()
    print_results(result)