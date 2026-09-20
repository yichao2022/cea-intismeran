"""Sync numbers in the .tex files to the freshly regenerated canonical results.

Compares the committed output/canonical_results.json (git HEAD) with the new one,
builds old-formatted -> new-formatted string pairs, reports every occurrence in
manuscript.tex / supplementary.tex, and (with --write) applies the ones that are
unambiguous.

    python sync_tex_numbers.py            # report only
    python sync_tex_numbers.py --write    # apply
"""
import json
import subprocess
import sys

TEX_FILES = ["manuscript.tex", "supplementary.tex"]
JSON_PATH = "output/canonical_results.json"


def load_old() -> dict:
    raw = subprocess.run(["git", "show", f"HEAD:{JSON_PATH}"],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(raw)


def walk(node, path=""):
    if isinstance(node, dict):
        for k, v in node.items():
            yield from walk(v, f"{path}/{k}")
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from walk(v, f"{path}[{i}]")
    elif isinstance(node, (int, float)) and not isinstance(node, bool):
        yield path, float(node)


def formats(v: float) -> list[str]:
    out = []
    if abs(v) >= 1000:
        out.append(f"{v:,.0f}")
    out += [f"{v:.2f}", f"{v:.1f}", f"{v:.0f}"]
    return out


def main() -> int:
    old = dict(walk(load_old(), "old"))
    new = dict(walk(json.load(open(JSON_PATH)), "new"))
    mapping = {}
    for path, ov in old.items():
        nv = new.get(path)
        if nv is None or abs(nv - ov) < 1e-9:
            continue
        for o, n in zip(formats(ov), formats(nv)):
            if o != n:
                mapping.setdefault(o, set()).add(n)

    ambiguous = {o: ns for o, ns in mapping.items() if len(ns) > 1}
    clean = {o: next(iter(ns)) for o, ns in mapping.items() if len(ns) == 1}
    print(f"可映射的旧值 {len(clean)} 个，歧义 {len(ambiguous)} 个\n")

    total = 0
    for fname in TEX_FILES:
        text = open(fname).read()
        hits = []
        for o, n in sorted(clean.items(), key=lambda kv: -len(kv[0])):
            c = text.count(o)
            if c:
                hits.append((o, n, c))
        if not hits:
            continue
        print(f"=== {fname} ===")
        for o, n, c in hits:
            print(f"  {o:>12} -> {n:<12} x{c}")
            total += c
        if "--write" in sys.argv:
            for o, n, _ in hits:
                text = text.replace(o, n)
            open(fname, "w").write(text)
    print(f"\n共 {total} 处")
    if ambiguous:
        print("\n歧义（需人工判断）:")
        for o, ns in sorted(ambiguous.items()):
            print(f"  {o} -> {sorted(ns)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
