#!/usr/bin/env python3
"""Build blinded copies of the manuscript and supplement from the live sources.

Anonymisation only: every other character is untouched. The sources in the
repository stay non-blinded (single source of truth); this script derives the
blinded PDFs on demand into a scratch directory and copies the results to the
repository root.

    python3 make_blind.py           # writes manuscript_blind.pdf,
                                    # supplementary_blind.pdf, manuscript_full_blind.pdf
"""
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
BUILD = pathlib.Path("/tmp/blind_build")
REPLACEMENTS = {
    # (pattern, replacement, required) — required=False tolerates patterns that a
    # future revision removes from the source (e.g. the Declarations blocks were
    # dropped once they moved to the title page).
    "manuscript.tex": [
        # author block -> anonymised
        (r"""\author{Yichao Jin, Ph.D. \\
  School of Economic, Political and Policy Sciences, \\
  University of Texas at Dallas, Richardson, TX, USA \\
  \texttt{Yichao.Jin@UTDallas.edu} \\
  ORCID: 0009-0003-7667-5143}""",
         r"\author{Blinded for peer review}", True),
        # repository URL (appears in the data- and code-availability statements)
        (r"\url{https://github.com/yichao2022/cea-intismeran}",
         r"\texttt{[repository URL withheld for double-blind review]}", False),
        (r"\textbf{Author contributions:} Yichao Jin conceived and designed the study, "
         r"developed the model, conducted the analysis, and wrote the manuscript.",
         r"\textbf{Author contributions:} Withheld for double-blind review.", False),
    ],
    "supplementary.tex": [
        ("% Yichao Jin, PhD", "% Blinded for peer review", True),
        (r"\author{Yichao Jin, Ph.D.}", r"\author{Blinded for peer review}", True),
        (r"\url{https://github.com/yichao2022/cea-intismeran} (tag \texttt{v1.0-submission})",
         r"\texttt{[repository URL withheld for double-blind review]} "
         r"(tag \texttt{v1.0-submission})", True),
    ],
}
# patterns that must not survive into the blinded PDFs
LEAKS = ["yichao", "UTDallas", "utdallas", "0009-0003-7667-5143",
         "University of Texas at Dallas"]


def build() -> None:
    if BUILD.exists():
        shutil.rmtree(BUILD)
    BUILD.mkdir(parents=True)
    for png in ROOT.glob("*.png"):
        shutil.copy(png, BUILD / png.name)
    for extra in ("refs.bib", "provenance_table.tex"):   # \input / \bibliography deps
        shutil.copy(ROOT / extra, BUILD / extra)

    for name, pairs in REPLACEMENTS.items():
        text = (ROOT / name).read_text()
        for old, new, required in pairs:
            n = text.count(old)
            if n == 0:
                if required:
                    sys.exit(f"FAIL: pattern not found in {name}: {old[:60]}")
                print(f"  {name}: pattern absent, skipped -> {old[:50]}")
                continue
            text = text.replace(old, new)
            print(f"  {name}: {n} anonymised -> {new[:60]}")
        (BUILD / name).write_text(text)

    for f in ("manuscript", "supplementary"):
        for cmd in (["pdflatex", "-interaction=nonstopmode", f"{f}.tex"],
                    ["bibtex", f],
                    ["pdflatex", "-interaction=nonstopmode", f"{f}.tex"],
                    ["pdflatex", "-interaction=nonstopmode", f"{f}.tex"]):
            subprocess.run(cmd, cwd=BUILD, capture_output=True)
        errs = len(re.findall(r"^!", (BUILD / f"{f}.log").read_text(), re.M))
        print(f"  {f}: {errs} LaTeX errors")

    subprocess.run(["pdfunite", "manuscript.pdf", "supplementary.pdf", "manuscript_full.pdf"],
                   cwd=BUILD, capture_output=True)
    for src, dst in (("manuscript.pdf", "manuscript_blind.pdf"),
                     ("supplementary.pdf", "supplementary_blind.pdf"),
                     ("manuscript_full.pdf", "manuscript_full_blind.pdf")):
        shutil.copy(BUILD / src, ROOT / dst)
        txt = subprocess.run(["pdftotext", str(ROOT / dst), "-"],
                             capture_output=True, text=True).stdout
        hits = {p: txt.lower().count(p.lower()) for p in LEAKS}
        bad = {k: v for k, v in hits.items() if v}
        print(f"  {dst}: {'LEAKS ' + str(bad) if bad else 'no identifying strings'}")
        if bad:
            sys.exit(1)


if __name__ == "__main__":
    build()
    print("blinded build complete")
