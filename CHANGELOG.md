# Changelog

## v0.1.0 (unreleased)

  - Initial release: next-generation matrices (`T`, `Σ`, `K_L`, `K`) from ModelingToolkit
    systems and Catalyst reaction networks following Diekmann, Heesterbeek & Roberts (2010).
  - Closed-form `R₀` for `1 × 1`, rank-one and `2 × 2` irreducible blocks; numeric evaluation
    otherwise. Type reproduction numbers, small-domain matrices, decomposition validation.
  - Quarto vignettes (Julia engine) and Lean 4 / Mathlib proofs of the underlying linear algebra.
