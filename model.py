"""
cea_intismeran/model.py
Partitioned Survival Model for Intismeran + Keytruda vs Keytruda alone
in adjuvant melanoma (Stage IIB-IV, completely resected).

Reference: KEYNOTE-942/mRNA-4157-P201 5-year follow-up (2026 ASCO)
          INTerpath-001 Phase 3 top-line (Merck 2026-08-19)

Cost sources:
  - Keytruda WAC: $24,544/q6w dose (GoodRx, May 2026; Keytruda.com)
  - Intismeran: ~$200K/course (Jefferies analyst estimate)
  - Sequencing: <$1,000 (WES, industry 2025)
  - Admin: CMS Physician Fee Schedule
  - Metastatic melanoma: published CEA literature

Utility sources:
  - DFS: 0.83 (Bensimon 2019, KEYNOTE-054 EQ-5D)
  - DM: 0.65 (melanoma CEA literature)
  - AE disutility: 0.05 (Bensimon 2019)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Tuple, Optional


# ── Survival distributions ──────────────────────────────────────────────────

def weibull_surv(t: np.ndarray, scale: float, shape: float) -> np.ndarray:
    """Weibull survival S(t) = exp(-(t/scale)^shape)"""
    return np.exp(-((t / scale) ** shape))


def fit_weibull_from_survrate(
    t_known: float, s_known: float, shape: float = 1.2
) -> Tuple[float, float]:
    """Fit Weibull scale from survival rate at known timepoint."""
    scale = t_known / ((-np.log(s_known)) ** (1 / shape))
    return scale, shape


# ── Parameter container ─────────────────────────────────────────────────────

@dataclass
class ModelParams:
    """Base-case model inputs. Replace with literature values as available."""

    # ── Clinical (KEYNOTE-942 / INTerpath-001) ──
    rfs_5yr_combo: float = 0.688       # 5-year RFS, combo arm
    rfs_5yr_pembro: float = 0.491      # 5-year RFS, pembro alone
    rfs_hr: float = 0.51               # HR for RFS (95% CI 0.294-0.887)
    dmfs_hr: float = 0.411             # HR for DMFS (95% CI 0.200-0.843)
    os_5yr_combo: float = 0.922        # 5-year OS, combo (exploratory)
    os_5yr_pembro: float = 0.713       # 5-year OS, pembro alone

    # ── Costs (USD, 2026) ──
    # Keytruda 400mg q6w × 9 doses: $24,544/dose (GoodRx May 2026)
    cost_keytruda_year: float = 220_896
    # Intismeran: Jefferies estimate ~$200K/course
    cost_intismeran_course: float = 200_000
    # Whole-exome sequencing for neoantigen identification
    cost_sequencing: float = 1_000
    # IV infusion admin per cycle (CMS Physician Fee Schedule)
    cost_admin_per_cycle: float = 500
    # Incremental AE management cost for combo arm (published CEA)
    cost_ae_combo: float = 15_000
    # Monthly metastatic melanoma treatment cost (literature)
    cost_dm_month: float = 12_000

    # ── Utilities (EQ-5D) ──
    util_dfs: float = 0.83              # Bensimon 2019, KEYNOTE-054
    util_dm: float = 0.65               # metastatic melanoma literature
    util_disutility_ae: float = 0.05    # one-time AE disutility

    # ── Model settings ──
    cycle_length: float = 1 / 12        # 1 month
    time_horizon: float = 40.0          # years (lifetime)
    discount_rate: float = 0.03         # 3%/year

    # ── Derived survival parameters ──
    rfs_scale_combo: float = field(init=False)
    rfs_scale_pembro: float = field(init=False)
    os_scale_combo: float = field(init=False)
    os_scale_pembro: float = field(init=False)
    weibull_shape: float = 1.2          # ponytail: fixed; estimate from KM digitization

    def __post_init__(self):
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


def sample_psa_params(n: int = 1000, seed: int = 42) -> list:
    """Generate PSA parameter draws from assumed distributions."""
    rng = np.random.default_rng(seed)
    draws = []

    for _ in range(n):
        p = ModelParams()

        # Clinical: log-normal on HRs
        log_hr_rfs = np.log(0.51)
        se_log_hr_rfs = (np.log(0.887) - np.log(0.294)) / (2 * 1.96)
        p.rfs_hr = np.exp(rng.normal(log_hr_rfs, se_log_hr_rfs))

        # Costs: gamma distribution (mean, SE ~20% of mean)
        p.cost_keytruda_year = rng.gamma(25, 220_896 / 25)
        p.cost_intismeran_course = rng.gamma(25, 200_000 / 25)
        p.cost_dm_month = rng.gamma(25, 12_000 / 25)

        # Utilities: beta distribution on logit scale
        def sample_beta(mean, se):
            alpha = mean * (mean * (1 - mean) / se**2 - 1)
            beta = (1 - mean) * (mean * (1 - mean) / se**2 - 1)
            return rng.beta(max(alpha, 1), max(beta, 1))

        p.util_dfs = sample_beta(0.83, 0.03)
        p.util_dm = sample_beta(0.65, 0.04)

        # Survival rates: beta on 5-year rates
        p.rfs_5yr_combo = sample_beta(0.688, 0.04)
        p.rfs_5yr_pembro = sample_beta(0.491, 0.05)
        p.os_5yr_combo = sample_beta(0.922, 0.03)
        p.os_5yr_pembro = sample_beta(0.713, 0.06)

        draws.append(p)

    return draws


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
        if arm == "combo":
            pfs = weibull_surv(self.t, self.p.rfs_scale_combo, self.p.weibull_shape)
            os = weibull_surv(self.t, self.p.os_scale_combo, self.p.weibull_shape)
        else:
            pfs = weibull_surv(self.t, self.p.rfs_scale_pembro, self.p.weibull_shape)
            os = weibull_surv(self.t, self.p.os_scale_pembro, self.p.weibull_shape)
        return pfs, os

    def _state_proportions(self, pfs: np.ndarray, os: np.ndarray) -> dict:
        return {
            "dfs": pfs,
            "dm": np.maximum(0, os - pfs),
            "death": 1 - os,
        }

    def run_arm(self, arm: str) -> dict:
        pfs, os = self._surv_curves(arm)
        states = self._state_proportions(pfs, os)
        disc = np.exp(-self.p.discount_rate * (self.t - self.p.cycle_length / 2))

        # QALYs
        qaly = (states["dfs"] * self.p.util_dfs + states["dm"] * self.p.util_dm) * disc
        total_qaly = np.sum(qaly) * self.p.cycle_length

        # Costs
        cycle_cost = np.zeros(self.n_cycles)

        if arm == "combo":
            treatment_months = 12
            total_drug = self.p.cost_intismeran_course + self.p.cost_keytruda_year
            monthly_tx = total_drug / treatment_months
            first_year = self.t <= 1.0
            cycle_cost[first_year] = monthly_tx * self.p.cycle_length * 12
            cycle_cost[0] += self.p.cost_sequencing
            cycle_cost[0] += self.p.cost_ae_combo
        else:
            monthly_tx = self.p.cost_keytruda_year / 12
            first_year = self.t <= 1.0
            cycle_cost[first_year] = monthly_tx * self.p.cycle_length * 12

        cycle_cost += self.p.cost_admin_per_cycle
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


# ── Results ──────────────────────────────────────────────────────────────────

def base_case() -> dict:
    return PartitionedSurvivalModel(ModelParams()).run()


def run_psa(n: int = 1000) -> dict:
    """Run PSA and return ICER distribution."""
    params_list = sample_psa_params(n)
    icers = []
    costs = {"combo": [], "pembro": []}
    qalys = {"combo": [], "pembro": []}

    for p in params_list:
        model = PartitionedSurvivalModel(p)
        r = model.run()
        costs["combo"].append(r["combo"]["total_cost"])
        costs["pembro"].append(r["pembro"]["total_cost"])
        qalys["combo"].append(r["combo"]["total_qaly"])
        qalys["pembro"].append(r["pembro"]["total_qaly"])
        icers.append(r["icer"])

    return {
        "icers": np.array(icers),
        "costs": costs,
        "qalys": qalys,
        "mean_icer": np.mean(icers),
        "median_icer": np.median(icers),
        "ci_icer": np.percentile(icers, [2.5, 97.5]),
        "p_cost_effective_100k": np.mean(np.array(icers) <= 100_000),
        "p_cost_effective_150k": np.mean(np.array(icers) <= 150_000),
    }


def print_results(result: dict):
    r = result
    print(f"{'='*65}")
    print(f"  CEA: Intismeran + Keytruda vs Keytruda (Adjuvant Melanoma)")
    print(f"  Reference: KEYNOTE-942 5-year follow-up (2026 ASCO)")
    print(f"  Costs: literature-based | Utilities: Bensimon 2019")
    print(f"{'='*65}")
    for arm_name in ["pembro", "combo"]:
        arm = r[arm_name]
        print(f"\n  {arm['arm'].upper()}:")
        print(f"    Total QALYs:  {arm['total_qaly']:>8.2f}")
        print(f"    Total Cost:  ${arm['total_cost']:>10,.0f}")

    print(f"\n  {'─'*55}")
    print(f"  Incremental Cost:  ${r['delta_cost']:>10,.0f}")
    print(f"  Incremental QALYs: {r['delta_qaly']:>10.2f}")
    print(f"  ICER:             ${r['icer']:>10,.0f} / QALY")
    print(f"  {'─'*55}\n")


def print_psa_results(psa: dict):
    print(f"  PSA Results ({len(psa['icers'])} iterations):")
    print(f"  Mean ICER:     ${psa['mean_icer']:>10,.0f}")
    print(f"  Median ICER:   ${psa['median_icer']:>10,.0f}")
    print(f"  95% CI:        ${psa['ci_icer'][0]:>10,.0f} — ${psa['ci_icer'][1]:>10,.0f}")
    print(f"  P(CE @ $100K):  {psa['p_cost_effective_100k']:.1%}")
    print(f"  P(CE @ $150K):  {psa['p_cost_effective_150k']:.1%}")


if __name__ == "__main__":
    print("── Base Case ──")
    result = base_case()
    print_results(result)

    print("── PSA (1000 iterations) ──")
    psa = run_psa(1000)
    print_psa_results(psa)

    # Save to repo for paper
    np.save("/Users/cary/cea-intismeran/output/psa_icers.npy", psa["icers"])
    print("\n  PSA results saved to output/")