#!/usr/bin/env bash
# Render the Lean sources to Markdown with mdgen and to HTML with pandoc.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mkdir -p docs/generated
lake build ReproductiveNumbersProofs
lake exe mdgen ReproductiveNumbersProofs docs/generated

# One combined document, in reading order.
COMBINED="docs/generated/ReproductiveNumbersProofs.md"
{
  cat docs/preamble.md
  for f in Restriction RankOne BlockTriangular TwoByTwo Examples; do
    printf '\n\n'
    cat "docs/generated/${f}.md"
  done
} > "$COMBINED"

pandoc "$COMBINED" \
  --standalone \
  --toc \
  --mathjax \
  --metadata title="ReproductiveNumbers.jl: machine-checked mathematics" \
  --css style.css \
  --output docs/ReproductiveNumbersProofs.html

test -s docs/ReproductiveNumbersProofs.html
printf 'Generated %s\n' docs/ReproductiveNumbersProofs.html
