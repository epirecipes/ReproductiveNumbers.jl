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

/-! ## SEIR with demography

Infected states `(E, I)`, `Ė = β S I / N - (σ + μ) E`, `İ = σ E - (γ + μ) I`, linearised at
`S = N`, exactly as in the package's second vignette. Only `E` is a state-at-infection, so
`K` is the scalar `β σ / ((σ + μ)(γ + μ))` even though `K_L` is `2 × 2`. -/
section SEIR

variable (β σ γ μ : ℝ)

def seirT : Matrix (Fin 2) (Fin 2) ℝ := !![0, β; 0, 0]
def seirSg : Matrix (Fin 2) (Fin 2) ℝ := !![-(σ + μ), 0; σ, -(γ + μ)]
def seirSgInv : Matrix (Fin 2) (Fin 2) ℝ :=
  !![-1 / (σ + μ), 0; -σ / ((σ + μ) * (γ + μ)), -1 / (γ + μ)]
/-- Selection matrix of the single state-at-infection `E`. -/
def seirE : Matrix (Fin 2) (Fin 1) ℝ := !![1; 0]
/-- The next-generation matrix with large domain, as computed by the package. -/
def seirKL : Matrix (Fin 2) (Fin 2) ℝ :=
  !![β * σ / ((σ + μ) * (γ + μ)), β / (γ + μ); 0, 0]
/-- The next-generation matrix, `R₀` as a `1 × 1` matrix. -/
def seirK : Matrix (Fin 1) (Fin 1) ℝ := !![β * σ / ((σ + μ) * (γ + μ))]

theorem seirSg_mul_inv (hσ : σ + μ ≠ 0) (hγ : γ + μ ≠ 0) :
    seirSg σ γ μ * seirSgInv σ γ μ = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [seirSg, seirSgInv, Matrix.mul_apply, Fin.sum_univ_two] <;> field_simp <;> ring

theorem seirSg_inv (hσ : σ + μ ≠ 0) (hγ : γ + μ ≠ 0) :
    (seirSg σ γ μ)⁻¹ = seirSgInv σ γ μ :=
  Matrix.inv_eq_right_inv (seirSg_mul_inv σ γ μ hσ hγ)

/-- `K_L = -T Σ⁻¹` for the SEIR model with demography. -/
theorem seir_KL (hσ : σ + μ ≠ 0) (hγ : γ + μ ≠ 0) :
    -(seirT β * (seirSg σ γ μ)⁻¹) = seirKL β σ γ μ := by
  rw [seirSg_inv σ γ μ hσ hγ]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [seirT, seirSgInv, seirKL] <;> ring

/-- Restricting to the state-at-infection gives the scalar `β σ / ((σ + μ)(γ + μ))`. -/
theorem seir_K : (seirE)ᵀ * seirKL β σ γ μ * seirE = seirK β σ γ μ := by
  ext i j
  fin_cases i; fin_cases j
  simp [seirE, seirKL, seirK, Matrix.mul_apply, Fin.sum_univ_two]

/-- `K_L` factorises through `E`, which is what makes `K` and `K_L` share their non-zero
spectrum (`mul_comm_hasEigenvalue_iff`). -/
theorem seir_KL_factor : seirE * ((seirE)ᵀ * seirKL β σ γ μ) = seirKL β σ γ μ := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [seirE, seirKL, Matrix.mul_apply, Fin.sum_univ_two]

/-- The non-zero eigenvalues of `K_L` are those of `K`. -/
theorem seir_eigen_iff {ν : ℝ} (hν : ν ≠ 0) :
    (∃ x : Fin 2 → ℝ, x ≠ 0 ∧ seirKL β σ γ μ *ᵥ x = ν • x) ↔
      (∃ y : Fin 1 → ℝ, y ≠ 0 ∧ seirK β σ γ μ *ᵥ y = ν • y) := by
  have hKL : seirKL β σ γ μ = seirE * ((seirE)ᵀ * seirKL β σ γ μ) :=
    (seir_KL_factor β σ γ μ).symm
  have hK : seirK β σ γ μ = ((seirE)ᵀ * seirKL β σ γ μ) * seirE := (seir_K β σ γ μ).symm
  conv_lhs => rw [hKL]
  conv_rhs => rw [hK]
  exact (mul_comm_hasEigenvalue_iff ((seirE)ᵀ * seirKL β σ γ μ) seirE hν).symm

/-- The characteristic polynomial of `K_L` is `ν (ν - R₀)`. -/
theorem seir_charpoly (ν : ℝ) :
    Matrix.det (ν • (1 : Matrix (Fin 2) (Fin 2) ℝ) - seirKL β σ γ μ) =
      ν * (ν - β * σ / ((σ + μ) * (γ + μ))) := by
  rw [Matrix.det_fin_two]
  simp [Matrix.one_fin_two, seirKL]
  ring

end SEIR

/-! ## Ross–Macdonald (vector-borne transmission)

Infected states `(I_H, I_V)`: humans are infected by mosquitoes at rate `a b I_V / N_H` per
susceptible human, mosquitoes by humans at rate `a c I_H / N_H` per susceptible mosquito;
humans recover at rate `γ` and mosquitoes die at rate `μ_V`. Linearised at `S_H = N_H`,
`S_V = N_V`. Both states are states-at-infection and `K` has a zero diagonal, so
`R₀ = √(K₁₂ K₂₁)`, as in the package's fourth vignette. (Identifiers are ASCII: `NH`, `NV`,
`gam`, `muV`.) -/
section RossMacdonald

variable (a b c NH NV gam muV : ℝ)

def rmT : Matrix (Fin 2) (Fin 2) ℝ := !![0, a * b; a * c * NV / NH, 0]
def rmSg : Matrix (Fin 2) (Fin 2) ℝ := !![-gam, 0; 0, -muV]
def rmSgInv : Matrix (Fin 2) (Fin 2) ℝ := !![-1 / gam, 0; 0, -1 / muV]
/-- The next-generation matrix: `K₁₂ = a b / μ_V`, `K₂₁ = a c N_V / (N_H γ)`. -/
def rmK : Matrix (Fin 2) (Fin 2) ℝ := !![0, a * b / muV; a * c * NV / (NH * gam), 0]

theorem rmSg_mul_inv (hg : gam ≠ 0) (hm : muV ≠ 0) : rmSg gam muV * rmSgInv gam muV = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rmSg, rmSgInv, Matrix.mul_apply, Fin.sum_univ_two] <;> field_simp

theorem rmSg_inv (hg : gam ≠ 0) (hm : muV ≠ 0) : (rmSg gam muV)⁻¹ = rmSgInv gam muV :=
  Matrix.inv_eq_right_inv (rmSg_mul_inv gam muV hg hm)

/-- `K = -T Σ⁻¹`; both rows of `T` are non-zero so `K_L = K`. -/
theorem rm_K (hg : gam ≠ 0) (hm : muV ≠ 0) :
    -(rmT a b c NH NV * (rmSg gam muV)⁻¹) = rmK a b c NH NV gam muV := by
  rw [rmSg_inv gam muV hg hm]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [rmT, rmSgInv, rmK] <;> ring

/-- `R₀ = √(a² b c N_V / (N_H γ μ_V))`. -/
theorem rm_R0 (hbc : 0 ≤ (a * b / muV) * (a * c * NV / (NH * gam))) :
    spectralRadius₂ 0 (a * b / muV) (a * c * NV / (NH * gam)) 0 =
      Real.sqrt (a ^ 2 * b * c * NV / (NH * gam * muV)) := by
  rw [spectralRadius₂_zero_diag _ _ hbc]
  congr 1
  ring

end RossMacdonald

/-! ## SEI with two latent categories (Diekmann et al. 2010, section 2.1)

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

/-! ## Two-sex or vector–host transmission (Diekmann et al. 2010, section 4.1)

With no within-type transmission `K = !![0, k₁₂; k₂₁, 0]` and `R₀ = √(k₁₂ k₂₁)`. -/
section TwoType

theorem two_type_R0 (k₁₂ k₂₁ : ℝ) (h : 0 ≤ k₁₂ * k₂₁) :
    spectralRadius₂ 0 k₁₂ k₂₁ 0 = Real.sqrt (k₁₂ * k₂₁) :=
  spectralRadius₂_zero_diag k₁₂ k₂₁ h

/-- With the explicit entries of section 4.1. -/
theorem sti_R0 (β₁ β₂ ν₁ ν₂ γ₁ γ₂ μ : ℝ)
    (h : 0 ≤ (β₁ * ν₂ / ((ν₂ + μ) * (γ₁ + μ))) * (β₂ * ν₁ / ((ν₁ + μ) * (γ₂ + μ)))) :
    spectralRadius₂ 0 (β₁ * ν₂ / ((ν₂ + μ) * (γ₁ + μ))) (β₂ * ν₁ / ((ν₁ + μ) * (γ₂ + μ))) 0 =
      Real.sqrt ((β₁ * ν₂ / ((ν₂ + μ) * (γ₁ + μ))) * (β₂ * ν₁ / ((ν₁ + μ) * (γ₂ + μ)))) :=
  spectralRadius₂_zero_diag _ _ h

end TwoType

end

end ReproductiveNumbersProofs
