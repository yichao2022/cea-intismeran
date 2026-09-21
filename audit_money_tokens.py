"""Audit: every money amount and key rate in the .tex files must be a value the
canonical JSONs actually produce. Flags tokens that no canonical output contains.

Usage: python audit_money_tokens.py [--decimals]
  --decimals  also scan plain decimal tokens (x.xx / x.x%) outside money strings
"""
import json
import os
import re
import sys

canon = json.load(open("output/canonical_results.json"))
dec = json.load(open("output/decomposition_v2.json"))
dsa = json.load(open("output/dsa_canonical.json"))
psa = json.load(open("output/psa_v2_results.json"))

# every canonical output the manuscript draws on — a token missing here is hand-typed
CANON_FILES = ["output/canonical_results.json", "output/decomposition_v2.json",
               "output/dsa_canonical.json", "output/psa_v2_results.json"]


def walk(o, out):
    if isinstance(o, dict):
        for v in o.values():
            walk(v, out)
    elif isinstance(o, list):
        for v in o:
            walk(v, out)
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        out.append(float(o))
    elif isinstance(o, str):
        out.append(o)


vals = []
for doc in (canon, dec, dsa, psa):
    walk(doc, vals)
nums = [v for v in vals if isinstance(v, float)]

known = set()
for v in nums:
    for fmt in (f"{v:,.0f}", f"{v:.2f}", f"{v:.1f}", f"{v:.3f}",
                f"{100 * v:.1f}", f"{100 * v:.2f}"):
        known.add(fmt)

MONEY = re.compile(r"(?<![\d.,])\d{1,3}(?:,\d{3})+(?![\d])")
DECIMAL = re.compile(r"(?<![\d.,])\d+\.\d{1,3}(?![\d])")

# Numbers reach the .tex from three shapes the walker above misses: string values,
# dict keys, and already-formatted strings. Scan the raw JSON text so a token that
# really is in a canonical output never gets flagged (2026-09-21: 58/58 were this).
for _f in CANON_FILES:
    known.update(MONEY.findall(open(_f).read()))

flag = {}
for f in ("manuscript.tex", "supplementary.tex", "cover_letter_pharmacoeconomics.tex",
          "title_page.tex"):
    for i, line in enumerate(open(f), 1):
        if line.lstrip().startswith("%"):
            continue
        for m in MONEY.finditer(line):
            tok = m.group(0)
            if tok not in known:
                flag.setdefault(tok, []).append(f"{f}:{i}")
        if "--decimals" in sys.argv:
            for m in DECIMAL.finditer(line):
                tok = m.group(0)
                if tok in known:
                    continue
                ctx = line[max(0, m.start() - 12):m.start()]
                if any(k in ctx for k in ("p=", "p =", "n=", "=", "CI", "SE", "e.g.", "vs", "\\cite")):
                    continue
                flag.setdefault(tok, []).append(f"{f}:{i}")

print(f"canonical scalars: {len(nums)} | flagged tokens: {len(flag)}")
for tok, where in sorted(flag.items(), key=lambda kv: -len(kv[1])):
    print(f"  {tok:<12} x{len(where):<3} {where[0]}")
