#!/usr/bin/env python3
"""Produce the blinded submission files from the live sources.

The canonical paths (manuscript.pdf, supplementary.pdf, manuscript_full.pdf) hold the
blinded builds, because those are the files that get linked and uploaded. No
non-blinded copies are kept: the sources in the repository are already non-blinded,
so the plain build can be regenerated at any time with pdflatex.

Anonymisation is only: the author block, the repository URLs, the author-contributions
statement and the supplement's author line. Everything else is copied unchanged.

    python3 make_blind.py
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
        (r"""\author{Yichao Jin, Ph.D. \\
  School of Economic, Political and Policy Sciences, \\
  University of Texas at Dallas, Richardson, TX, USA \\
  \texttt{Yichao.Jin@UTDallas.edu} \\
  ORCID: 0009-0003-7667-5143}""",
         r"\author{Blinded for peer review}", True),
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
LEAKS = ["yichao", "UTDallas", "utdallas", "0009-0003-7667-5143",
         "University of Texas at Dallas"]


def main() -> None:
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
        if errs:
            sys.exit(1)
    subprocess.run(["pdfunite", "manuscript.pdf", "supplementary.pdf", "manuscript_full.pdf"],
                   cwd=BUILD, capture_output=True)

    for name in ("manuscript.pdf", "supplementary.pdf", "manuscript_full.pdf"):
        shutil.copy(BUILD / name, ROOT / name)
        txt = subprocess.run(["pdftotext", str(ROOT / name), "-"],
                             capture_output=True, text=True).stdout
        bad = {p: txt.lower().count(p.lower()) for p in LEAKS}
        bad = {k: v for k, v in bad.items() if v}
        ok = "Blinded for peer review" in txt
        print(f"  {name}: {'LEAKS ' + str(bad) if bad else 'no identifying strings'}"
              f"{', anonymised' if ok else ', NOT ANONYMISED'}")
        if bad or not ok:
            sys.exit(1)
    print("done: canonical PDFs blinded")


if __name__ == "__main__":
    main()
