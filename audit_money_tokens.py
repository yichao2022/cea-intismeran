"""Audit: every money amount and key rate in the .tex files must be a value the
canonical JSONs actually produce. Flags tokens that no canonical output contains.

Usage: python audit_money_tokens.py [--decimals]
  --decimals  also scan plain decimal tokens (x.xx / x.x%) outside money strings
"""
import json
import glob
import os
import re
import sys

canon = json.load(open("output/canonical_results.json"))
dec = json.load(open("output/decomposition_v2.json"))
dsa = json.load(open("output/dsa_canonical.json"))
psa = json.load(open("output/psa_v2_results.json"))

# every canonical output the manuscript draws on — a token missing here is hand-typed
CANON_FILES = ["output/canonical_results.json", "output/decomposition_v2.json",
               "output/dsa_canonical.json", "output/psa_v2_results.json",
               "output/printed_values.json"]


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
                f"{100 * v:.1f}", f"{100 * v:.2f}",
                # magnitudes too: a .tex prints "$-$671,440" as two tokens, so the
                # signed value never matches the printed one
                f"{abs(v):,.0f}"):
        known.add(fmt)

MONEY = re.compile(r"(?<![\d.,])\d{1,3}(?:,\d{3})+(?![\d])")
DECIMAL = re.compile(r"(?<![\d.,])\d+\.\d{1,3}(?![\d])")

# Numbers reach the .tex from three shapes the walker above misses: string values,
# dict keys, and already-formatted strings. Scan the raw JSON text so a token that
# really is in a canonical output never gets flagged (2026-09-21: 58/58 were this).
for _f in CANON_FILES:
    known.update(MONEY.findall(open(_f).read()))

# Declared inputs: amounts the manuscript legitimately states as parameters that no
# output JSON reproduces -- drug prices, administration and management costs, the life
# table. They are read from where they are declared, and reported as their own class so
# an unattributed amount is not buried among them.
INPUT_FILES = sorted(glob.glob("*.py")) + sorted(glob.glob("data/*.json"))
inputs = set()
for _f in INPUT_FILES:
    for _m in re.finditer(r"\d[\d_]*(?:\.\d+)?", open(_f).read()):
        _v = float(_m.group(0).replace("_", ""))
        if _v >= 1000:
            inputs.add(f"{_v:,.0f}")

# Amounts that are not model outputs: literature values, external price lists, and
# counts the money regex happens to match. Each entry names its provenance, and the
# audit prints them as their own class -- a value that needs an entry here but has no
# source is a manuscript problem, not an allowlist entry.
EXTERNAL = {
    "75,206": "prior exploratory analysis ICER (Discussion comparison)",
    "75,000": "published interim-data ICER (Discussion comparison)",
    "1,137": "INTerpath-001 enrolment (literature)",
    "443,000": "break-even price rounded from the analytic threshold 443,167",
    "12,336": "pembrolizumab WAC per 200 mg dose (external price list)",
    "38,050": "published pembrolizumab CEA ICER, lower bound (literature)",
    "78,400": "published pembrolizumab CEA ICER, upper bound (literature)",
    "4,000": "iteration count in a figure caption, not an amount",
}

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

def klass(tok: str) -> str:
    if tok in inputs:
        return "input"
    if tok in EXTERNAL:
        return "external"
    return "UNATTRIBUTED"


counts = {k: sum(1 for t in flag if klass(t) == k) for k in ("input", "external", "UNATTRIBUTED")}
print(f"canonical scalars: {len(nums)} | declared-input amounts: {len(inputs)} | "
      f"flagged tokens: {len(flag)} ({counts['input']} declared inputs, "
      f"{counts['external']} external, {counts['UNATTRIBUTED']} unattributed)")
for tok, where in sorted(flag.items(), key=lambda kv: -len(kv[1])):
    k = klass(tok)
    note = f"  # {EXTERNAL[tok]}" if k == "external" else ""
    print(f"  {tok:<12} x{len(where):<3} {k:<12} {where[0]}{note}")
if counts["UNATTRIBUTED"]:
    sys.exit(1)  # an amount nobody can source is a failure, not a warning
