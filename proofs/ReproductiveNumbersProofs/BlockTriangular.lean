import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Data.Matrix.Block

/-!
# Block-triangular next-generation matrices

`ReproductiveNumbers.jl` permutes a next-generation matrix `K` into block upper triangular
form using the strongly connected components of its non-zero pattern, computes the dominant
eigenvalue of each irreducible diagonal block in closed form, and returns the maximum. This
is legitimate because the characteristic polynomial of a block-triangular matrix is the
product of the characteristic polynomials of its diagonal blocks, so its spectrum is the
union of the spectra of the blocks; for a non-negative matrix the spectral radius is
therefore the largest of the block spectral radii. Multi-strain models with no interaction
between strains are the typical case.

We prove the two-block statement; the general case follows by induction on the number of
blocks.
-/

namespace ReproductiveNumbersProofs

open Matrix

variable {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
variable {R : Type*} [CommRing R]

/-- `μ • 1 - fromBlocks A B 0 D` is again block upper triangular. -/
theorem charmatrix_fromBlocks (A : Matrix m m R) (B : Matrix m n R) (D : Matrix n n R)
    (μ : R) :
    μ • (1 : Matrix (m ⊕ n) (m ⊕ n) R) - fromBlocks A B 0 D =
      fromBlocks (μ • (1 : Matrix m m R) - A) (-B) 0 (μ • (1 : Matrix n n R) - D) := by
  rw [← fromBlocks_one, fromBlocks_smul, sub_eq_add_neg, fromBlocks_neg, fromBlocks_add]
  simp [sub_eq_add_neg]

/-- The characteristic polynomial of a block upper triangular matrix, evaluated at `μ`, is
the product of those of its diagonal blocks. -/
theorem det_charmatrix_fromBlocks (A : Matrix m m R) (B : Matrix m n R) (D : Matrix n n R)
    (μ : R) :
    det (μ • (1 : Matrix (m ⊕ n) (m ⊕ n) R) - fromBlocks A B 0 D) =
      det (μ • (1 : Matrix m m R) - A) * det (μ • (1 : Matrix n n R) - D) := by
  rw [charmatrix_fromBlocks, det_fromBlocks_zero₂₁]

variable {F : Type*} [Field F]

/-- `μ` is an eigenvalue (root of the characteristic polynomial) of a block upper triangular
matrix iff it is an eigenvalue of one of the diagonal blocks. Hence the spectral radius of a
non-negative block-triangular matrix is the maximum of the block spectral radii, which is
the rule used by `ReproductiveNumbers.spectral_radius`. -/
theorem det_charmatrix_fromBlocks_eq_zero_iff (A : Matrix m m F) (B : Matrix m n F)
    (D : Matrix n n F) (μ : F) :
    det (μ • (1 : Matrix (m ⊕ n) (m ⊕ n) F) - fromBlocks A B 0 D) = 0 ↔
      det (μ • (1 : Matrix m m F) - A) = 0 ∨ det (μ • (1 : Matrix n n F) - D) = 0 := by
  rw [det_charmatrix_fromBlocks, mul_eq_zero]

end ReproductiveNumbersProofs
