# Changelog

## v0.1.0 (unreleased)

  - Initial release: next-generation matrices (`T`, `Σ`, `K_L`, `K`) from ModelingToolkit
    systems and Catalyst reaction networks following Diekmann, Heesterbeek & Roberts (2010).
  - Closed-form `R₀` for `1 × 1`, rank-one and `2 × 2` irreducible blocks; numeric evaluation
    otherwise. Type reproduction numbers, small-domain matrices, decomposition validation.
  - Quarto vignettes (Julia engine) and Lean 4 / Mathlib proofs of the underlying linear algebra.
  - Transmission classification is recorded in the result and selectable by name; warnings
    for misfiled loss terms; polynomial infection-free steady states with Nemo; general
    small-domain factorisation and automatic fallback to `K_S`.
  - `suggest_infected`, `mean_sojourn_times`, `sensitivities`, `elasticities`,
    `perron_vectors`, `effective_reproduction_number`, a ForwardDiff constructor from
    functions, and a Latexify extension.
  - Effective reproduction number `R_t` (closed form, at a state or time, along solutions,
    as an observed variable, with time-varying rates), size-aware algorithm routing with
    retries, and `abbreviate`/`expand_definitions` for readable matrices.
  - Steady-state regression tests (SteadyStateDiffEq, NonlinearSolve), models of sections
    4.2 and 4.3 of Diekmann et al., eight further vignettes, Lean proofs of the rank-one
    converse and the block-triangular spectrum.
