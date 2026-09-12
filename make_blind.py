#!/usr/bin/env python3
"""Produce the blinded submission files from the live sources.

The canonical paths (manuscript.pdf, supplementary.pdf, manuscript_full.pdf) are the
blinded versions, because those are the files that get linked and uploaded; the
non-blinded builds stay beside them as *_unblinded.pdf for the authors' records.
Sources in the repository remain non-blinded.

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
BLIND = pathlib.Path("/tmp/blind_build")
PLAIN = pathlib.Path("/tmp/plain_build")
# (pattern, replacement, required)

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


def prepare(dest: pathlib.Path, anonymise: bool) -> None:
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    for png in ROOT.glob("*.png"):
        shutil.copy(png, dest / png.name)
    for extra in ("refs.bib", "provenance_table.tex"):   # \input / \bibliography deps
        shutil.copy(ROOT / extra, dest / extra)
    for name, pairs in REPLACEMENTS.items():
        text = (ROOT / name).read_text()
        if anonymise:
            for old, new, required in pairs:
                n = text.count(old)
                if n == 0:
                    if required:
                        sys.exit(f"FAIL: pattern not found in {name}: {old[:60]}")
                    print(f"  {name}: pattern absent, skipped -> {old[:50]}")
                    continue
                text = text.replace(old, new)
                print(f"  {name}: {n} anonymised -> {new[:60]}")
        (dest / name).write_text(text)


def compile_in(where: pathlib.Path) -> None:
    for f in ("manuscript", "supplementary"):
        for cmd in (["pdflatex", "-interaction=nonstopmode", f"{f}.tex"],
                    ["bibtex", f],
                    ["pdflatex", "-interaction=nonstopmode", f"{f}.tex"],
                    ["pdflatex", "-interaction=nonstopmode", f"{f}.tex"]):
            subprocess.run(cmd, cwd=where, capture_output=True)
        errs = len(re.findall(r"^!", (where / f"{f}.log").read_text(), re.M))
        print(f"  {where.name}/{f}: {errs} LaTeX errors")
        if errs:
            sys.exit(1)
    subprocess.run(["pdfunite", "manuscript.pdf", "supplementary.pdf", "manuscript_full.pdf"],
                   cwd=where, capture_output=True)


def text_of(pdf: pathlib.Path) -> str:
    return subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True).stdout


def main() -> None:
    print("building blinded:"); prepare(BLIND, anonymise=True); compile_in(BLIND)
    print("building non-blinded:"); prepare(PLAIN, anonymise=False); compile_in(PLAIN)

    outputs = [("manuscript.pdf", "manuscript.pdf"), ("supplementary.pdf", "supplementary.pdf"),
               ("manuscript_full.pdf", "manuscript_full.pdf")]
    for src, dst in outputs:
        shutil.copy(BLIND / src, ROOT / dst)
    for src, dst in (("manuscript.pdf", "manuscript_unblinded.pdf"),
                     ("supplementary.pdf", "supplementary_unblinded.pdf"),
                     ("manuscript_full.pdf", "manuscript_full_unblinded.pdf")):
        shutil.copy(PLAIN / src, ROOT / dst)

    # the canonical files must be anonymous, the *_unblinded copies must not
    for src, _ in outputs:
        txt = text_of(ROOT / src)
        bad = {p: txt.lower().count(p.lower()) for p in LEAKS}
        bad = {k: v for k, v in bad.items() if v}
        print(f"  {src}: {'LEAKS ' + str(bad) if bad else 'no identifying strings'}")
        if bad:
            sys.exit(1)
        if "Blinded for peer review" not in txt:
            sys.exit(f"FAIL: {src} is not anonymised")
    plain = text_of(ROOT / "manuscript_unblinded.pdf")
    if "Yichao Jin" not in plain:
        sys.exit("FAIL: manuscript_unblinded.pdf lost the author block")
    print("  manuscript_unblinded.pdf: author block present (as intended)")


if __name__ == "__main__":
    main()
    print("done: canonical PDFs blinded, *_unblinded.pdf kept for the record")
