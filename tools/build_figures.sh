#!/bin/sh
# Author: Xu Chen, University of Washington
# Build every labs/*/figures/*.tex (TikZ, standalone) to .pdf (vector), .svg and .png (for README/notebooks).
# Needs a LaTeX installation (not part of the Docker image; the built figures are committed).
#   usage: tools/build_figures.sh [lab-filter]
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
for tex in "$ROOT"/labs/*${1}*/figures/*.tex; do
  dir=$(dirname "$tex"); base=$(basename "$tex" .tex)
  (cd "$dir" && pdflatex -interaction=nonstopmode -halt-on-error "$base.tex" >/dev/null) || { echo "FAILED $tex (see $dir/$base.log)"; exit 1; }
  pdftocairo -svg "$dir/$base.pdf" "$dir/$base.svg"
  pdftocairo -png -r 220 -singlefile "$dir/$base.pdf" "$dir/$base"
  rm -f "$dir/$base.aux" "$dir/$base.log"
  echo "built ${tex#$ROOT/}"
done
