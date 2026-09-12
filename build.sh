#!/bin/bash
# Single build entry point.
#
# The sources are blinded, so a plain compile produces the submission files. Only the
# title page and the cover letter carry author identity. The script fails if an author
# string ever reappears in the manuscript or supplement.
set -e
cd "$(dirname "$0")"

for f in manuscript supplementary title_page cover_letter_pharmacoeconomics; do
  pdflatex -interaction=nonstopmode "$f.tex" > /dev/null 2>&1
  bibtex "$f" > /dev/null 2>&1 || true
  pdflatex -interaction=nonstopmode "$f.tex" > /dev/null 2>&1
  pdflatex -interaction=nonstopmode "$f.tex" > /dev/null 2>&1
  errs=$(grep -c '^!' "$f.log" || true)
  printf '%-36s %s errors, %s pages\n' "$f" "$errs" \
    "$(pdfinfo "$f.pdf" | awk '/^Pages/{print $2}')"
  [ "$errs" = "0" ] || { echo "FAIL: $f has LaTeX errors"; exit 1; }
done

pdfunite manuscript.pdf supplementary.pdf manuscript_full.pdf
printf '%-36s %s pages\n' manuscript_full.pdf "$(pdfinfo manuscript_full.pdf | awk '/^Pages/{print $2}')"

echo
for f in manuscript.pdf supplementary.pdf manuscript_full.pdf; do
  hits=$(pdftotext "$f" - | grep -ci 'yichao\|utdallas\|0009-0003' || true)
  printf '%-36s %s\n' "$f" "$([ "$hits" = "0" ] && echo 'blinded' || echo "LEAK ($hits)")"
  [ "$hits" = "0" ] || { echo "FAIL: $f leaks author identity"; exit 1; }
done
for f in title_page.pdf cover_letter_pharmacoeconomics.pdf; do
  printf '%-36s %s\n' "$f" "$(pdftotext "$f" - | grep -qi 'yichao' && echo 'AUTHOR DETAILS (as intended)' || echo 'MISSING AUTHOR DETAILS')"
done
