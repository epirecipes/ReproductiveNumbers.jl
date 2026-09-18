# Formal proofs

The linear algebra that the package relies on is formalised in Lean 4 against Mathlib in
the `proofs/` directory of the repository. The rendered report, generated from the Lean
sources with [mdgen](https://github.com/Seasawher/mdgen) and pandoc, is bundled with this
site: **[ReproductiveNumbers.jl: machine-checked mathematics](proofs/ReproductiveNumbersProofs.html)**.

Every displayed declaration type-checks against Mathlib v4.30.0, and `proofs/scripts/check.sh`
rejects `sorry`, `admit` and new axioms.

## What is proved

**`Restriction.lean`**

  - `mul_comm_eigenvector`, `mul_comm_hasEigenvalue_iff`: for matrices `A : m × n` and
    `B : n × m`, `A B` and `B A` have the same non-zero eigenvalues. Applied with
    `A = Eᵀ K_L` and `B = E` this shows that the next-generation matrix `K = Eᵀ K_L E` and
    the large-domain matrix `K_L` share their non-zero spectrum, so both give `R₀`.
  - `det_charmatrix_mul_comm`: the same statement for roots of the characteristic
    polynomials, via `det(1 − AB) = det(1 − BA)`.
  - `det_jacobian_eq_det_ngm_sub_one`: `det(T + Σ) = det(−TΣ⁻¹ − 1) · det(−Σ)`, so the
    infection-free steady state loses stability exactly when `1` is an eigenvalue of `K_L`.

**`RankOne.lean`**

  - `eigenvalue_vecMulVec`: every eigenvalue of a rank-one matrix `u vᵀ` is `0` or its trace.
    This is the closed form used when the states-at-infection are entered in fixed
    proportions.
  - `minors_vecMulVec`: all `2 × 2` minors of `u vᵀ` vanish, the criterion the software
    tests.

**`TwoByTwo.lean`**

  - `disc_nonneg`, `spectralRadius₂_charpoly`, `det_charmatrix_spectralRadius₂`: the
    closed form `(a + d + √((a − d)² + 4bc)) / 2` is real and is an eigenvalue of
    `!![a, b; c, d]`.
  - `le_spectralRadius₂`, `abs_le_spectralRadius₂`: it dominates every real eigenvalue in
    absolute value, so it is the spectral radius of a non-negative `2 × 2` matrix.
  - `spectralRadius₂_zero_diag`: with a zero diagonal it reduces to `√(bc)` (vector–host and
    two-sex models).
  - `one_lt_spectralRadius₂_iff`: `R₀ > 1` iff `trace K > 2` or the characteristic
    polynomial is negative at `1`.

**`Examples.lean`**

  - The SIR, SEIR, two-latent-stage (section 2.1) and two-type (section 4.1) models of the
    paper and the vignettes, with explicit inverses of `Σ`, `K_L`, `K`, and `R₀`.

## Building the proofs

```sh
cd proofs
lake exe cache get      # Mathlib build artefacts, first time only
./scripts/check.sh      # build, reject placeholders, render docs/ReproductiveNumbersProofs.html
```
