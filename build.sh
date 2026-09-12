#!/bin/bash
# Single build entry point.
#
# manuscript.pdf / supplementary.pdf / manuscript_full.pdf must stay blinded, so this
# never leaves a plain pdflatex output at those paths: it compiles the documents that
# are not blinded (title page, cover letter) and then lets make_blind.py rebuild the
# three blinded PDFs from the sources.
set -e
cd "$(dirname "$0")"

for f in title_page cover_letter_pharmacoeconomics; do
  pdflatex -interaction=nonstopmode "$f.tex" > /dev/null 2>&1
  pdflatex -interaction=nonstopmode "$f.tex" > /dev/null 2>&1
  errs=$(grep -c '^!' "$f.log" || true)
  printf '%-36s %s errors, %s pages\n' "$f" "$errs" \
    "$(pdfinfo "$f.pdf" | awk '/^Pages/{print $2}')"
done

python3 make_blind.py

echo
for f in manuscript.pdf supplementary.pdf manuscript_full.pdf title_page.pdf \
         cover_letter_pharmacoeconomics.pdf; do
  printf '%-36s %s pages  %s\n' "$f" "$(pdfinfo "$f" | awk '/^Pages/{print $2}')" \
    "$(pdftotext "$f" - | grep -qE 'Yichao|yichao' && echo 'AUTHOR DETAILS' || echo 'blinded')"
done
