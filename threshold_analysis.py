"""Threshold analysis: ICER at intismeran prices + value-based pricing.
GP-constrained model (final)."""
import sys, json
sys.path.insert(0, "/Users/cary/cea-intismeran")
from model import ModelParams
from rerun_primary import CEAModelV2

def icer_at_price(price: float) -> tuple:
    p = ModelParams(constraint_general_pop=True, cost_intismeran=price)
    r = CEAModelV2(p).run()
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    return dc, dq, dc / dq

def find_price_for_icer(target, low=0, high=5_000_000, tol=100):
    for _ in range(60):
        mid = (low + high) / 2
        _, _, icer = icer_at_price(mid)
        if abs(icer - target) < tol:
            break
        if icer < target:
            low = mid
        else:
            high = mid
    return mid, icer

if __name__ == "__main__":
    # Base case
    _, _, icer_base = icer_at_price(200_000)
    print(f"Base case ($200K/course): ICER = ${icer_base:,.0f}/QALY\n")
    # Price ladder
    print("Price → ICER")
    for price in [100_000, 200_000, 300_000, 400_000, 500_000, 600_000, 800_000, 1_000_000]:
        _, _, icer = icer_at_price(price)
        print(f"  ${price:>7,}/course → ${icer:>8,.0f}/QALY")
    # Value-based pricing
    print("\nValue-based pricing:")
    for target in [100_000, 150_000]:
        price, icer = find_price_for_icer(target)
        print(f"  ICER=${target:>6,}/QALY → price=${price:>8,.0f}/course (actual ICER ${icer:,.0f})")