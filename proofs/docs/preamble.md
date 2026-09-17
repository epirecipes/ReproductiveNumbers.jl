# Overview

These Lean 4 files, checked against Mathlib, formalise the linear algebra behind
`ReproductiveNumbers.jl`. The Julia package builds a next-generation matrix by splitting the
linearised infected subsystem of a compartmental model into transmissions `T` and
transitions `Σ` and forming `K_L = -T Σ⁻¹` and `K = Eᵀ K_L E` (Diekmann, Heesterbeek and
Roberts, *J. R. Soc. Interface* 2010). It then returns `R₀` as the spectral radius of `K`,
symbolically when a closed form exists. The proofs below cover exactly the steps the
software relies on:

 1. **Restriction** (`Restriction.lean`): `K` and `K_L` have the same non-zero eigenvalues,
    so restricting to the states-at-infection loses nothing; and `det(T + Σ) = 0` exactly
    when `1` is an eigenvalue of `K_L`, which is the `R₀ = 1` threshold.
 2. **Rank one** (`RankOne.lean`): when the states-at-infection are entered in fixed
    proportions the next-generation matrix has rank one and `R₀` is its trace.
 3. **Two by two** (`TwoByTwo.lean`): the closed form
    `R₀ = (tr K + √((a - d)² + 4bc)) / 2` used for `2 × 2` blocks is a real eigenvalue, it
    dominates every other real eigenvalue, it reduces to `√(bc)` for vector–host and two-sex
    models, and `R₀ > 1` iff `tr K > 2` or the characteristic polynomial is negative at `1`.
 4. **Examples** (`Examples.lean`): the SIR, SEIR, two-latent-stage and two-type models used
    in the package's tests and vignettes, worked through explicitly.

The documentation is generated from the Lean sources with
[mdgen](https://github.com/Seasawher/mdgen) and rendered with pandoc; every displayed
declaration type-checks, and `scripts/check.sh` refuses `sorry`, `admit` and new axioms.
