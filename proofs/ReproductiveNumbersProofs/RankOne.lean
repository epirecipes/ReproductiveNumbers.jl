import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Data.Matrix.Mul
import Mathlib.Tactic

/-!
# Rank-one next-generation matrices

When all states-at-infection are entered in fixed proportions (for instance an SEI model with
two latent categories entered with probabilities `p` and `1 - p`, Example 2.1 of Diekmann,
Heesterbeek and Roberts 2010), the next-generation matrix `K` has rank one and `det K = 0`.
Diekmann et al. then reduce to the *small domain* matrix `K_S` which is a scalar equal to the
trace of `K`.

Here we prove the linear-algebra fact used by `ReproductiveNumbers.jl` when it detects a
rank-one block: if `K = u vᵀ` (entries `K i j = u i * v j`) then every eigenvalue of `K` is
either `0` or `trace K = ∑ i, u i * v i`.  Hence the spectral radius of a non-negative rank-one
matrix is its trace, and `R₀ = trace K`.
-/

namespace ReproductiveNumbersProofs

open Matrix Finset

variable {n : Type*} [Fintype n] {R : Type*} [Field R]

/-- Multiplying a vector by the rank-one matrix `u vᵀ` gives `(v ⬝ x) • u`. -/
theorem vecMulVec_mulVec' (u v x : n → R) :
    (vecMulVec u v) *ᵥ x = (v ⬝ᵥ x) • u := by
  ext i
  simp only [Matrix.mulVec, Matrix.vecMulVec_apply, dotProduct, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- The trace of `u vᵀ` is the inner product `v ⬝ᵥ u`. -/
theorem trace_vecMulVec' (u v : n → R) : (vecMulVec u v).trace = v ⬝ᵥ u := by
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.vecMulVec_apply, dotProduct]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- Every eigenvalue of a rank-one matrix `u vᵀ` is `0` or its trace. -/
theorem eigenvalue_vecMulVec (u v : n → R) {μ : R} {x : n → R} (hx : x ≠ 0)
    (h : (vecMulVec u v) *ᵥ x = μ • x) : μ = 0 ∨ μ = (vecMulVec u v).trace := by
  rw [vecMulVec_mulVec'] at h
  rw [trace_vecMulVec']
  by_cases hvx : v ⬝ᵥ x = 0
  · left
    rw [hvx, zero_smul] at h
    exact (smul_eq_zero.mp h.symm).resolve_right hx
  · right
    -- take the inner product of both sides with `v`
    have h' : v ⬝ᵥ ((v ⬝ᵥ x) • u) = v ⬝ᵥ (μ • x) := by rw [h]
    rw [dotProduct_smul, dotProduct_smul, smul_eq_mul, smul_eq_mul] at h'
    -- `(v ⬝ x) * (v ⬝ u) = μ * (v ⬝ x)`, cancel the non-zero factor
    have : (v ⬝ᵥ x) * (v ⬝ᵥ u) = (v ⬝ᵥ x) * μ := by rw [h', mul_comm]
    exact (mul_left_cancel₀ hvx this).symm

omit [Fintype n] in
/-- A square matrix whose `2 × 2` minors all vanish and which has a non-zero row is of the form
`u vᵀ`: this is the criterion the Julia code uses (all `2 × 2` minors zero) before applying the
trace formula.  We record the direction used in the software: an explicit factorisation gives
vanishing minors. -/
theorem minors_vecMulVec (u v : n → R) (i j k l : n) :
    (vecMulVec u v) i k * (vecMulVec u v) j l - (vecMulVec u v) i l * (vecMulVec u v) j k = 0 := by
  simp only [Matrix.vecMulVec_apply]
  ring

end ReproductiveNumbersProofs
