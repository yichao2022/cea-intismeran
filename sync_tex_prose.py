"""Sync the hand-typed numbers that the pipeline does not own.

results_pipeline.py rewrites the generated tables (tab_base, tab_price, tab_dsa,
tab_scenarios, tab_survival_combos); everything else -- the PSA summary table,
the abstract's rounded thresholds, the one-way DSA prose and the Supplementary
PSA Results block -- is typed by hand and drifts silently on every rerun.

Every pair asserts a single match, so a changed anchor fails loudly instead of
skipping. Re-running on an already-synced manuscript is a no-op (0 matches).

Usage:  python sync_tex_prose.py
"""
import sys

MANUSCRIPT = [
    # abstract: rounded threshold prices and the worst-case scenario ICER
    ("The deterministic threshold prices were approximately \\$468K and \\$701K.",
     "The deterministic threshold prices were approximately \\$443K and \\$646K."),
    ("(\\$104,591/QALY)", "(\\$104,971/QALY)"),
    # highlights: break-even per-course price
    ("up to about \\$468,000 per course", "up to about \\$443,000 per course"),
    # abstract / highlights / PSA section: P(CE) at $50K and cost-saving share
    ("from 61.1\\% at \\$50,000/QALY", "from 63.9\\% at \\$50,000/QALY"),
    ("15.8\\% of simulations were cost-saving", "16.8\\% of simulations were cost-saving"),
    ("and 15.8\\% were cost-saving", "and 16.8\\% were cost-saving"),
    # tab:psa -- mean / median / 95% CI
    ("Incremental cost  & \\$164,003 & \\$200,679 & $-$\\$677,115 to \\$612,313",
     "Incremental cost  & \\$128,915 & \\$166,863 & $-$\\$671,440 to \\$524,694"),
    ("Incremental QALYs & 4.95      & 4.78      & 2.87--7.90",
     "Incremental QALYs & 4.36      & 4.18      & 2.56--7.10"),
    # one-way DSA prose: ranking AND values both moved with the rerun
    ("the largest NMB changes arose from the discount rate (\\$37,614--\\$48,727/QALY "
     "across 0\\%--5\\%), the combination-arm OS parameter (\\$28,283--\\$42,814/QALY "
     "across its $-40\\%$ to $+40\\%$ range), and the \\mbox{pembrolizumab}-arm OS "
     "parameter (\\$58,503/QALY to dominant across its $-20\\%$ to $+20\\%$ range)",
     "the largest NMB changes arose from the combination-arm OS parameter "
     "(\\$13,725--\\$40,486/QALY across its $-40\\%$ to $+40\\%$ range), the "
     "\\mbox{pembrolizumab}-arm OS parameter (\\$57,365/QALY to dominant across its "
     "$-20\\%$ to $+20\\%$ range), and the discount rate (\\$33,050--\\$47,535/QALY "
     "across 0\\%--5\\%)"),
    ("produced ICERs of \\$33,816--\\$51,008/QALY", "produced ICERs of \\$30,102--\\$49,849/QALY"),
    ("the maximum finite ICER was \\$58,503/QALY", "the maximum finite ICER was \\$57,365/QALY"),
    # discussion: hazard-convergence ICER
    ("yields the highest ICER of \\$104,591/QALY", "yields the highest ICER of \\$104,971/QALY"),
]

SUPPLEMENT = [
    # PSA Results block + CE plane caption
    ("\\item Mean $\\Delta$QALY: 4.95 (95\\% CI 2.87--7.90).",
     "\\item Mean $\\Delta$QALY: 4.36 (95\\% CI 2.56--7.10)."),
    ("\\item Mean $\\Delta$Cost: \\$164,003 (95\\% CI \\textendash\\$677,115 to \\$612,313).",
     "\\item Mean $\\Delta$Cost: \\$128,915 (95\\% CI \\textendash\\$671,440 to \\$524,694)."),
    ("Probability cost-effective: 61.1\\% at \\$50K", "Probability cost-effective: 63.9\\% at \\$50K"),
    ("15.8\\% of draws were cost-saving", "16.8\\% of draws were cost-saving"),
    ("15.8\\% of iterations fall in the southeast quadrant",
     "16.8\\% of iterations fall in the southeast quadrant"),
]

changed = 0
unresolved = []
for path, pairs in (("manuscript.tex", MANUSCRIPT), ("supplementary.tex", SUPPLEMENT)):
    text = open(path).read()
    for old, new in pairs:
        if new in text:
            continue  # already synced
        n = text.count(old)
        if n == 0:
            # Neither the old nor the new string is present: the anchor drifted (the
            # sentence was reworded or regenerated). Failing here is the point -- a
            # silently skipped pair is how a stale hand-typed number survives a rerun.
            unresolved.append(f"{path}: {old[:70]!r}")
            continue
        assert n == 1, f"{path}: anchor matched {n} times: {old[:60]}"
        text = text.replace(old, new)
        changed += 1
        print(f"{path}: {old[:58]}... -> {new[:44]}...")
    open(path, "w").write(text)

if unresolved:
    print("DRIFTED anchors (fix the pair or the sentence):", *unresolved, sep="\n  ")
    sys.exit(1)
print(f"{changed} replacement(s)")
