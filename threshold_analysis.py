"""
Threshold analysis: ICER at intismeran prices $200K-$2M + value-based pricing.
Uses current CEAModel API.
"""
import sys
sys.path.insert(0, "/Users/cary/cea-intismeran")
from model import ModelParams, CEAModel


def icer_at_price(price: float) -> tuple:
    p = ModelParams(cost_intismeran=price)
    r = CEAModel(p).run()
    d_cost = r["combo"]["cost"] - r["pembro"]["cost"]
    d_qaly = r["combo"]["qaly"] - r["pembro"]["qaly"]
    return d_cost, d_qaly, d_cost / d_qaly


def find_price_for_icer(target_icer, price_low=0, price_high=20_000_000, tol=500):
    """Binary search for intismeran price giving target ICER."""
    mid = None
    for _ in range(60):
        mid = (price_low + price_high) / 2
        _, _, icer = icer_at_price(mid)
        if abs(icer - target_icer) < tol:
            break
        if icer < target_icer:
            price_low = mid
        else:
            price_high = mid
    return mid, icer


if __name__ == "__main__":
    # Base case
    _, _, icer_base = icer_at_price(200_000)
    print(f"Base case ($200K/course): ICER = ${icer_base:,.0f}/QALY\n")

    # Price points within $200K-$2M range
    print("Intismeran price → ICER")
    for price in [200_000, 300_000, 500_000, 750_000, 1_000_000, 1_500_000, 2_000_000]:
        _, _, icer = icer_at_price(price)
        print(f"  ${price:>7,}/course → ${icer:>8,.0f}/QALY")

    # Value-based pricing
    print("\nValue-based pricing:")
    for target in [50_000, 100_000, 150_000, 200_000]:
        price, icer = find_price_for_icer(target)
        print(f"  ICER=${target:>6,}/QALY → price=${price:>8,.0f}/course (actual ICER ${icer:,.0f})")