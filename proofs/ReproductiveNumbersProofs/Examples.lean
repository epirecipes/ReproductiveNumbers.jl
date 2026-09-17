import ReproductiveNumbersProofs.Restriction
import ReproductiveNumbersProofs.RankOne
import ReproductiveNumbersProofs.TwoByTwo

/-!
# Worked examples

The next-generation matrices that `ReproductiveNumbers.jl` derives for the models used in
its vignettes and tests, checked symbolically in Lean.  Each example fixes the transmission
matrix `T` and the transition matrix `Σ` (written `Sg` because `Σ` is a Lean token), exhibits
the inverse of `Σ`, and computes `K_L = -T Σ⁻¹` and `K = Eᵀ K_L E`.
-/

namespace ReproductiveNumbersProofs

open Matrix

noncomputable section

/-! ## SIR

For `İ = β S I / N - γ I` linearised at `S = N`: `T = [β]`, `Σ = [-γ]`, `R₀ = β / γ`. -/
section SIR

variable (β γ : ℝ)

def sirT : Matrix (Fin 1) (Fin 1) ℝ := !![β]
def sirSg : Matrix (Fin 1) (Fin 1) ℝ := !![-γ]

theorem sirSg_inv (hγ : γ ≠ 0) : (sirSg γ)⁻¹ = !![-1 / γ] := by
  apply Matrix.inv_eq_right_inv
  ext i j
  fin_cases i; fin_cases j
  simp [sirSg, Matrix.mul_apply]
  field_simp

theorem sir_ngm (hγ : γ ≠ 0) : -(sirT β * (sirSg γ)⁻¹) = !![β / γ] := by
  rw [sirSg_inv γ hγ]
  ext i j
  fin_cases i; fin_cases j
  simp [sirT, Matrix.mul_apply]
  ring

end SIR

/-! ## SEIR

Infected states `(E, I)`, `Ė = β S I / N - σ E`, `İ = σ E - γ I`.  Only `E` is a
state-at-infection, so `K` is the scalar `β / γ` even though `K_L` is `2 × 2`. -/
section SEIR

variable (β σ γ : ℝ)

def seirT : Matrix (Fin 2) (Fin 2) ℝ := !![0, β; 0, 0]
def seirSg : Matrix (Fin 2) (Fin 2) ℝ := !![-σ, 0; σ, -γ]
def seirSgInv : Matrix (Fin 2) (Fin 2) ℝ := !![-1 / σ, 0; -1 / γ, -1 / γ]
/-- Selection matrix of the single state-at-infection `E`. -/
def seirE : Matrix (Fin 2) (Fin 1) ℝ := !![1; 0]

theorem seirSg_mul_inv (hσ : σ ≠ 0) (hγ : γ ≠ 0) : seirSg σ γ * seirSgInv σ γ = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [seirSg, seirSgInv, Matrix.mul_apply, Fin.sum_univ_two] <;> field_simp <;> ring

theorem seirSg_inv (hσ : σ ≠ 0) (hγ : γ ≠ 0) : (seirSg σ γ)⁻¹ = seirSgInv σ γ :=
  Matrix.inv_eq_right_inv (seirSg_mul_inv σ γ hσ hγ)

/-- The next-generation matrix with large domain of the SEIR model. -/
theorem seir_KL (hσ : σ ≠ 0) (hγ : γ ≠ 0) :
    -(seirT β * (seirSg σ γ)⁻¹) = !![β / γ, β / γ; 0, 0] := by
  rw [seirSg_inv σ γ hσ hγ]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [seirT, seirSgInv] <;> ring

/-- Restricting to the state-at-infection gives the scalar `β / γ`. -/
theorem seir_K :
    (seirE)ᵀ * !![β / γ, β / γ; 0, 0] * seirE = !![β / γ] := by
  ext i j
  fin_cases i; fin_cases j
  simp [seirE, Matrix.mul_apply, Fin.sum_univ_two]

/-- `K_L` factorises through `E`, which is what makes `K` and `K_L` share their non-zero
spectrum (`mul_comm_hasEigenvalue_iff`). -/
theorem seir_KL_factor :
    seirE * ((seirE)ᵀ * !![β / γ, β / γ; 0, 0]) = !![β / γ, β / γ; 0, 0] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [seirE, Matrix.mul_apply, Fin.sum_univ_two]

/-- The non-zero eigenvalues of `K_L` are those of `K = [β / γ]`. -/
theorem seir_eigen_iff {μ : ℝ} (hμ : μ ≠ 0) :
    (∃ x : Fin 2 → ℝ, x ≠ 0 ∧ !![β / γ, β / γ; 0, 0] *ᵥ x = μ • x) ↔
      (∃ y : Fin 1 → ℝ, y ≠ 0 ∧ !![β / γ] *ᵥ y = μ • y) := by
  have hKL : !![β / γ, β / γ; 0, 0] = seirE * ((seirE)ᵀ * !![β / γ, β / γ; 0, 0]) :=
    (seir_KL_factor β γ).symm
  have hK : !![β / γ] = ((seirE)ᵀ * !![β / γ, β / γ; 0, 0]) * seirE := (seir_K β γ).symm
  conv_lhs => rw [hKL]
  conv_rhs => rw [hK]
  exact (mul_comm_hasEigenvalue_iff ((seirE)ᵀ * !![β / γ, β / γ; 0, 0]) seirE hμ).symm

/-- The characteristic polynomial of `K_L` is `μ (μ - β / γ)`. -/
theorem seir_charpoly (μ : ℝ) :
    Matrix.det (μ • (1 : Matrix (Fin 2) (Fin 2) ℝ) - !![β / γ, β / γ; 0, 0]) = μ * (μ - β / γ) := by
  rw [Matrix.det_fin_two]
  simp [Matrix.one_fin_two]
  ring

end SEIR

/-! ## SEI with two latent categories (Diekmann et al. 2010, Example 2.1)

New infections enter `E₁` with probability `p` and `E₂` with probability `1 - p`; the
transmission matrix has rank one and `R₀ = trace K = p ν₁ β / ((ν₁ + μ)(γ + μ)) +
(1 - p) ν₂ β / ((ν₂ + μ)(γ + μ))`. -/
section TwoLatent

variable (β p ν₁ ν₂ γ μ : ℝ)

def latT : Matrix (Fin 3) (Fin 3) ℝ := !![0, 0, β * p; 0, 0, β * (1 - p); 0, 0, 0]
def latSg : Matrix (Fin 3) (Fin 3) ℝ :=
  !![-(ν₁ + μ), 0, 0; 0, -(ν₂ + μ), 0; ν₁, ν₂, -(γ + μ)]
def latSgInv : Matrix (Fin 3) (Fin 3) ℝ :=
  !![-1 / (ν₁ + μ), 0, 0;
     0, -1 / (ν₂ + μ), 0;
     -ν₁ / ((ν₁ + μ) * (γ + μ)), -ν₂ / ((ν₂ + μ) * (γ + μ)), -1 / (γ + μ)]

theorem latSg_mul_inv (h₁ : ν₁ + μ ≠ 0) (h₂ : ν₂ + μ ≠ 0) (h₃ : γ + μ ≠ 0) :
    latSg ν₁ ν₂ γ μ * latSgInv ν₁ ν₂ γ μ = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [latSg, latSgInv, Matrix.mul_apply, Fin.sum_univ_three] <;> field_simp <;> ring

theorem latSg_inv (h₁ : ν₁ + μ ≠ 0) (h₂ : ν₂ + μ ≠ 0) (h₃ : γ + μ ≠ 0) :
    (latSg ν₁ ν₂ γ μ)⁻¹ = latSgInv ν₁ ν₂ γ μ :=
  Matrix.inv_eq_right_inv (latSg_mul_inv ν₁ ν₂ γ μ h₁ h₂ h₃)

/-- The column and row of the rank-one factorisation `K_L = u vᵀ`. -/
def latU : Fin 3 → ℝ := ![β * p, β * (1 - p), 0]
def latV : Fin 3 → ℝ :=
  ![ν₁ / ((ν₁ + μ) * (γ + μ)), ν₂ / ((ν₂ + μ) * (γ + μ)), 1 / (γ + μ)]

/-- `K_L = -T Σ⁻¹` is the rank-one matrix `u vᵀ`. -/
theorem lat_KL (h₁ : ν₁ + μ ≠ 0) (h₂ : ν₂ + μ ≠ 0) (h₃ : γ + μ ≠ 0) :
    -(latT β p * (latSg ν₁ ν₂ γ μ)⁻¹) = vecMulVec (latU β p) (latV ν₁ ν₂ γ μ) := by
  rw [latSg_inv ν₁ ν₂ γ μ h₁ h₂ h₃]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [latT, latSgInv, latU, latV, Matrix.vecMulVec_apply] <;> field_simp

/-- Hence every non-zero eigenvalue of `K_L`, in particular `R₀`, equals the trace
`p ν₁ β / ((ν₁ + μ)(γ + μ)) + (1 - p) ν₂ β / ((ν₂ + μ)(γ + μ))`. -/
theorem lat_R0 (h₁ : ν₁ + μ ≠ 0) (h₂ : ν₂ + μ ≠ 0) (h₃ : γ + μ ≠ 0) {ρ : ℝ} (hρ : ρ ≠ 0)
    {x : Fin 3 → ℝ} (hx : x ≠ 0)
    (h : (-(latT β p * (latSg ν₁ ν₂ γ μ)⁻¹)) *ᵥ x = ρ • x) :
    ρ = β * p * ν₁ / ((ν₁ + μ) * (γ + μ)) + β * (1 - p) * ν₂ / ((ν₂ + μ) * (γ + μ)) := by
  rw [lat_KL β p ν₁ ν₂ γ μ h₁ h₂ h₃] at h
  rcases eigenvalue_vecMulVec _ _ hx h with h0 | htr
  · exact absurd h0 hρ
  · rw [htr, trace_vecMulVec']
    simp [latU, latV, dotProduct, Fin.sum_univ_three]
    ring

end TwoLatent

/-! ## Two-sex or vector–host transmission (Diekmann et al. 2010, Example 4.1)

With no within-type transmission `K = !![0, k₁₂; k₂₁, 0]` and `R₀ = √(k₁₂ k₂₁)`. -/
section TwoType

theorem two_type_R0 (k₁₂ k₂₁ : ℝ) (h : 0 ≤ k₁₂ * k₂₁) :
    spectralRadius₂ 0 k₁₂ k₂₁ 0 = Real.sqrt (k₁₂ * k₂₁) :=
  spectralRadius₂_zero_diag k₁₂ k₂₁ h

/-- With the explicit entries of Example 4.1. -/
theorem sti_R0 (β₁ β₂ ν₁ ν₂ γ₁ γ₂ μ : ℝ)
    (h : 0 ≤ (β₁ * ν₂ / ((ν₂ + μ) * (γ₁ + μ))) * (β₂ * ν₁ / ((ν₁ + μ) * (γ₂ + μ)))) :
    spectralRadius₂ 0 (β₁ * ν₂ / ((ν₂ + μ) * (γ₁ + μ))) (β₂ * ν₁ / ((ν₁ + μ) * (γ₂ + μ))) 0 =
      Real.sqrt ((β₁ * ν₂ / ((ν₂ + μ) * (γ₁ + μ))) * (β₂ * ν₁ / ((ν₁ + μ) * (γ₂ + μ)))) :=
  spectralRadius₂_zero_diag _ _ h

end TwoType

end

end ReproductiveNumbersProofs
