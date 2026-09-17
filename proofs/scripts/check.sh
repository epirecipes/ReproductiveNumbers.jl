#!/usr/bin/env bash
# Build every proof, refuse placeholders, and render the documentation.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if grep -R -n -E '^[[:space:]]*(sorry|admit|axiom|unsafe)([[:space:]]|$)|:=[[:space:]]*(sorry|by[[:space:]]+sorry)' \
  --include='*.lean' ReproductiveNumbersProofs ReproductiveNumbersProofs.lean; then
  echo "Forbidden proof placeholder or declaration found." >&2
  exit 1
fi

lake build ReproductiveNumbersProofs
./scripts/build-docs.sh
