"""
Threshold analysis: find the intismeran price that yields a target ICER.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/cary/cea-intismeran")
from model import ModelParams, PartitionedSurvivalModel

def find_price_for_icer(target_icer, price_low=0, price_high=500_000, tol=100):
    """Binary search for intismeran price that gives target ICER."""
    for _ in range(50):
        price_mid = (price_low + price_high) / 2
        p = ModelParams(cost_intismeran_course=price_mid)
        icer = PartitionedSurvivalModel(p).run()["icer"]
        if abs(icer - target_icer) < tol:
            break
        if icer < target_icer:
            price_low = price_mid
        else:
            price_high = price_mid
    return price_mid, icer

# At $150K/QALY threshold
price_150, icer_150 = find_price_for_icer(150_000)
print(f"ICER = $150,000/QALY → Intismeran price = ${price_150:,.0f}/course")
print(f"  (actual ICER at this price: ${icer_150:,.0f})")

# At $100K/QALY threshold
price_100, icer_100 = find_price_for_icer(100_000)
print(f"ICER = $100,000/QALY → Intismeran price = ${price_100:,.0f}/course")
print(f"  (actual ICER at this price: ${icer_100:,.0f})")

# At $200K/QALY threshold
price_200, icer_200 = find_price_for_icer(200_000)
print(f"ICER = $200,000/QALY → Intismeran price = ${price_200:,.0f}/course")
print(f"  (actual ICER at this price: ${icer_200:,.0f})")

# Price points for reference
print()
for test_price in [100_000, 150_000, 200_000, 250_000, 300_000]:
    p = ModelParams(cost_intismeran_course=test_price)
    icer = PartitionedSurvivalModel(p).run()["icer"]
    print(f"Price ${test_price:>6,}/course → ICER ${icer:>8,.0f}/QALY")