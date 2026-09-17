import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Data.Real.Sqrt
import Mathlib.Tactic

/-!
# The spectral radius of a non-negative `2 × 2` matrix

Next-generation matrices with two states-at-infection are ubiquitous: two host types
(Example 2.2), vector–host transmission, or heterosexual transmission (Example 4.1) in
Diekmann, Heesterbeek and Roberts (2010).  For such a matrix
`K = !![a, b; c, d]` with non-negative entries the paper's equation (2.12) gives

`R₀ = (trace K + √(trace K ^ 2 - 4 det K)) / 2`.

`ReproductiveNumbers.jl` uses exactly this formula for every irreducible `2 × 2` block.  We
verify that it is legitimate:

* the discriminant `trace² - 4 det = (a - d)² + 4 b c` is non-negative, so the square root is
  real;
* the formula is a root of the characteristic polynomial, i.e. an eigenvalue;
* it dominates every real eigenvalue in absolute value, so it is the spectral radius;
* when the diagonal vanishes (vector–host or two-sex models) it reduces to `√(b c)`;
* the threshold criterion `R₀ > 1 ↔ trace K > 2 ∨ 1 - trace K + det K < 0`.
-/

namespace ReproductiveNumbersProofs

open Matrix

/-- The discriminant of the characteristic polynomial of `!![a, b; c, d]`. -/
noncomputable def disc (a b c d : ℝ) : ℝ := (a + d) ^ 2 - 4 * (a * d - b * c)

/-- Equation (2.12) of Diekmann et al.: the dominant eigenvalue of a `2 × 2` matrix. -/
noncomputable def spectralRadius₂ (a b c d : ℝ) : ℝ :=
  ((a + d) + Real.sqrt (disc a b c d)) / 2

theorem disc_eq (a b c d : ℝ) : disc a b c d = (a - d) ^ 2 + 4 * (b * c) := by
  unfold disc; ring

/-- The discriminant is non-negative whenever the off-diagonal entries have the same sign. -/
theorem disc_nonneg (a b c d : ℝ) (hbc : 0 ≤ b * c) : 0 ≤ disc a b c d := by
  rw [disc_eq]
  positivity

theorem sq_sqrt_disc (a b c d : ℝ) (hbc : 0 ≤ b * c) :
    Real.sqrt (disc a b c d) ^ 2 = disc a b c d :=
  Real.sq_sqrt (disc_nonneg a b c d hbc)

/-- `spectralRadius₂` is a root of the characteristic polynomial `λ² - (a + d) λ + (a d - b c)`. -/
theorem spectralRadius₂_charpoly (a b c d : ℝ) (hbc : 0 ≤ b * c) :
    spectralRadius₂ a b c d ^ 2 - (a + d) * spectralRadius₂ a b c d + (a * d - b * c) = 0 := by
  have hs := sq_sqrt_disc a b c d hbc
  unfold spectralRadius₂
  unfold disc at hs ⊢
  nlinarith [hs]

/-- `spectralRadius₂ a b c d` is an eigenvalue of `!![a, b; c, d]`. -/
theorem det_charmatrix_spectralRadius₂ (a b c d : ℝ) (hbc : 0 ≤ b * c) :
    Matrix.det (spectralRadius₂ a b c d • (1 : Matrix (Fin 2) (Fin 2) ℝ) - !![a, b; c, d]) = 0 := by
  have h := spectralRadius₂_charpoly a b c d hbc
  rw [Matrix.det_fin_two, Matrix.one_fin_two]
  simp
  linear_combination h

/-- Both bounds `-s ≤ 2 μ - (a + d) ≤ s` for a root `μ`, where `s = √disc`. -/
theorem root_bounds (a b c d : ℝ) (hbc : 0 ≤ b * c) {μ : ℝ}
    (hμ : μ ^ 2 - (a + d) * μ + (a * d - b * c) = 0) :
    -Real.sqrt (disc a b c d) ≤ 2 * μ - (a + d) ∧
      2 * μ - (a + d) ≤ Real.sqrt (disc a b c d) := by
  have hs := sq_sqrt_disc a b c d hbc
  have hnn := Real.sqrt_nonneg (disc a b c d)
  have hsq : (2 * μ - (a + d)) ^ 2 = Real.sqrt (disc a b c d) ^ 2 := by
    rw [hs]; unfold disc; linear_combination 4 * hμ
  have h1 := sq_le_sq.mp (le_of_eq hsq)
  rw [abs_of_nonneg hnn] at h1
  exact abs_le.mp h1

/-- Every real eigenvalue `μ` of `!![a, b; c, d]` satisfies `μ ≤ spectralRadius₂ a b c d`. -/
theorem le_spectralRadius₂ (a b c d : ℝ) (hbc : 0 ≤ b * c) {μ : ℝ}
    (hμ : μ ^ 2 - (a + d) * μ + (a * d - b * c) = 0) : μ ≤ spectralRadius₂ a b c d := by
  have h := (root_bounds a b c d hbc hμ).2
  unfold spectralRadius₂
  linarith

/-- For non-negative entries the dominant eigenvalue is non-negative. -/
theorem spectralRadius₂_nonneg (a b c d : ℝ) (ha : 0 ≤ a) (hd : 0 ≤ d) (hbc : 0 ≤ b * c) :
    0 ≤ spectralRadius₂ a b c d := by
  unfold spectralRadius₂
  have := Real.sqrt_nonneg (disc a b c d)
  linarith

/-- For non-negative entries every real eigenvalue `μ` satisfies `|μ| ≤ spectralRadius₂`, so
`spectralRadius₂` is the spectral radius. -/
theorem abs_le_spectralRadius₂ (a b c d : ℝ) (ha : 0 ≤ a) (hd : 0 ≤ d) (hbc : 0 ≤ b * c) {μ : ℝ}
    (hμ : μ ^ 2 - (a + d) * μ + (a * d - b * c) = 0) : |μ| ≤ spectralRadius₂ a b c d := by
  have h := root_bounds a b c d hbc hμ
  have hnn := Real.sqrt_nonneg (disc a b c d)
  rw [abs_le]
  unfold spectralRadius₂
  constructor <;> linarith [h.1, h.2]

/-- Vector–host and two-sex models have a zero diagonal, and then `R₀ = √(b c)`: one
"generation" is a host-to-vector-to-host cycle counted as two infections. -/
theorem spectralRadius₂_zero_diag (b c : ℝ) (hbc : 0 ≤ b * c) :
    spectralRadius₂ 0 b c 0 = Real.sqrt (b * c) := by
  unfold spectralRadius₂ disc
  have h4 : (0 + 0) ^ 2 - 4 * (0 * 0 - b * c) = 2 ^ 2 * (b * c) := by ring
  rw [h4, Real.sqrt_mul (by positivity), Real.sqrt_sq (by norm_num)]
  ring

/-- The threshold criterion for a non-negative `2 × 2` next-generation matrix:
`R₀ > 1` iff `trace K > 2` or the characteristic polynomial is negative at `1`. -/
theorem one_lt_spectralRadius₂_iff (a b c d : ℝ) (hbc : 0 ≤ b * c) :
    1 < spectralRadius₂ a b c d ↔ 2 < a + d ∨ 1 - (a + d) + (a * d - b * c) < 0 := by
  have hs := sq_sqrt_disc a b c d hbc
  have hnn := Real.sqrt_nonneg (disc a b c d)
  unfold spectralRadius₂
  unfold disc at hs hnn ⊢
  set s := Real.sqrt ((a + d) ^ 2 - 4 * (a * d - b * c)) with hsdef
  constructor
  · intro h
    by_cases htr : 2 < a + d
    · exact Or.inl htr
    · right
      push Not at htr
      -- `2 - (a + d) < s` with both sides non-negative, so square
      have h1 : 2 - (a + d) < s := by linarith
      nlinarith [hs, h1, htr, hnn]
  · rintro (h | h)
    · linarith
    · by_cases htr : 2 < a + d
      · linarith
      · push Not at htr
        -- `(2 - (a + d))² < s²` and `s ≥ 0`, so `2 - (a + d) < s`
        have hlt : (2 - (a + d)) ^ 2 < s ^ 2 := by rw [hs]; nlinarith
        have : 2 - (a + d) < s := by
          by_contra hcon
          push Not at hcon
          nlinarith [hcon, hnn]
        linarith

end ReproductiveNumbersProofs
