# Formal proofs for ReproductiveNumbers.jl

Lean 4 / Mathlib formalisation of the matrix algebra used by the package. See
`docs/preamble.md` for the list of results, and `docs/ReproductiveNumbersProofs.html` for
the rendered version.

```sh
cd proofs
lake exe cache get      # download Mathlib build artefacts (first time only)
./scripts/check.sh      # build the proofs, reject placeholders, render docs
```

Rendering needs `pandoc`; `mdgen` is pulled in as a Lake dependency.
