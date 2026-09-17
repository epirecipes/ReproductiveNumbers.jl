import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

/-!
# Restricting the next-generation matrix to states-at-infection

Diekmann, Heesterbeek and Roberts (2010) define the next-generation matrix with
*large domain* `K_L = -T Σ⁻¹` and the next-generation matrix `K = -Eᵀ T Σ⁻¹ E`, where the
columns of `E` are the unit vectors of the states-at-infection (the non-zero rows of `T`).
They state that `K` and `K_L` have the same non-zero eigenvalues, so both give `R₀`.

Because every row of `T` outside the states-at-infection vanishes, `T = E Eᵀ T`, hence
`K_L = E · (Eᵀ (-T Σ⁻¹))` and `K = (Eᵀ (-T Σ⁻¹)) · E`.  Writing `B = E` and
`A = Eᵀ (-T Σ⁻¹)` this is the classical fact that `A B` and `B A` share their non-zero
eigenvalues, which is what we prove here, in eigenvector form and in determinant form.
-/

namespace ReproductiveNumbersProofs

open Matrix

variable {m n : Type*} [Fintype m] [Fintype n]
variable {R : Type*} [Field R]

/-- If `x` is an eigenvector of `A * B` for a non-zero eigenvalue `μ`, then `B *ᵥ x` is a
non-zero eigenvector of `B * A` for the same eigenvalue. -/
theorem mul_comm_eigenvector (A : Matrix m n R) (B : Matrix n m R) {μ : R} (hμ : μ ≠ 0)
    {x : m → R} (hx : x ≠ 0) (h : (A * B) *ᵥ x = μ • x) :
    B *ᵥ x ≠ 0 ∧ (B * A) *ᵥ (B *ᵥ x) = μ • (B *ᵥ x) := by
  refine ⟨?_, ?_⟩
  · intro hBx
    apply hx
    have h' : A *ᵥ (B *ᵥ x) = μ • x := by rw [Matrix.mulVec_mulVec]; exact h
    rw [hBx, Matrix.mulVec_zero] at h'
    exact (smul_eq_zero.mp h'.symm).resolve_left hμ
  · calc (B * A) *ᵥ (B *ᵥ x) = B *ᵥ ((A * B) *ᵥ x) := by
          simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc]
      _ = B *ᵥ (μ • x) := by rw [h]
      _ = μ • (B *ᵥ x) := Matrix.mulVec_smul _ _ _

/-- The non-zero eigenvalues of `A * B` and `B * A` coincide (eigenvector formulation). -/
theorem mul_comm_hasEigenvalue_iff (A : Matrix m n R) (B : Matrix n m R) {μ : R} (hμ : μ ≠ 0) :
    (∃ x : m → R, x ≠ 0 ∧ (A * B) *ᵥ x = μ • x) ↔
      (∃ y : n → R, y ≠ 0 ∧ (B * A) *ᵥ y = μ • y) := by
  constructor
  · rintro ⟨x, hx, h⟩
    exact ⟨B *ᵥ x, mul_comm_eigenvector A B hμ hx h⟩
  · rintro ⟨y, hy, h⟩
    exact ⟨A *ᵥ y, mul_comm_eigenvector B A hμ hy h⟩

variable [DecidableEq m] [DecidableEq n]

/-- Determinant formulation: for `μ ≠ 0`, `μ` is a root of the characteristic polynomial of
`A * B` iff it is a root of the characteristic polynomial of `B * A`. -/
theorem det_charmatrix_mul_comm (A : Matrix m n R) (B : Matrix n m R) {μ : R} (hμ : μ ≠ 0) :
    det (μ • (1 : Matrix m m R) - A * B) = 0 ↔ det (μ • (1 : Matrix n n R) - B * A) = 0 := by
  have h1 : μ • (1 : Matrix m m R) - A * B = μ • (1 - (μ⁻¹ • A) * B) := by
    rw [smul_sub, ← Matrix.smul_mul, smul_smul, mul_inv_cancel₀ hμ, one_smul]
  have h2 : (1 : Matrix n n R) - B * (μ⁻¹ • A) = μ⁻¹ • (μ • (1 : Matrix n n R) - B * A) := by
    rw [Matrix.mul_smul, smul_sub, smul_smul, inv_mul_cancel₀ hμ, one_smul]
  rw [h1, Matrix.det_smul, Matrix.det_one_sub_mul_comm, h2, Matrix.det_smul]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h | h
    · exact absurd h (pow_ne_zero _ hμ)
    · rcases mul_eq_zero.mp h with h | h
      · exact absurd h (pow_ne_zero _ (inv_ne_zero hμ))
      · exact h
  · intro h
    rw [h]; ring

/-- The Jacobian of the infected subsystem is `T + Σ = (K_L - 1) · (-Σ)` with `K_L = -T Σ⁻¹`.
Hence `det (T + Σ) = 0` exactly when `1` is an eigenvalue of `K_L`: the threshold `R₀ = 1`
coincides with the loss of stability of the infection-free steady state.  (`Σ` is a reserved
token in Lean, so the transition matrix is written `Sg` in the code.) -/
theorem det_jacobian_eq_det_ngm_sub_one (T Sg : Matrix n n R) (hS : IsUnit Sg.det) :
    det (T + Sg) = det (-(T * Sg⁻¹) - 1) * det (-Sg) := by
  have hfac : (-(T * Sg⁻¹) - 1) * (-Sg) = T + Sg := by
    rw [sub_mul, neg_mul, Matrix.mul_neg, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hS,
      Matrix.mul_one, neg_neg, one_mul, sub_neg_eq_add]
  rw [← Matrix.det_mul, hfac]

end ReproductiveNumbersProofs
