import numpy as np
d = np.load("/tmp/figures/psa_v2_draws.npz")
dq, dc = d["dq"], d["dc"]

print("n =", len(dq))
print(f"Incr cost:  mean={dc.mean():,.0f}  median={np.median(dc):,.0f}  95%CI=[{np.percentile(dc,2.5):,.0f}, {np.percentile(dc,97.5):,.0f}]")
print(f"Incr QALY:  mean={dq.mean():.2f}  median={np.median(dq):.2f}  95%CI=[{np.percentile(dq,2.5):.2f}, {np.percentile(dq,97.5):.2f}]")
print(f"NMB@100K: mean={(-dc+100_000*dq).mean():,.0f}")
print(f"NMB@150K: mean={(-dc+150_000*dq).mean():,.0f}")
mask = dc > 0
icers = dc[mask]/dq[mask]
print(f"ICER (cost-increasing only, n={mask.sum()}): mean={icers.mean():,.0f} median={np.median(icers):,.0f}")
print(f"  cost-saving draws: {(~mask).sum()} ({np.mean(~mask)*100:.1f}%)")