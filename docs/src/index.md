# ReproductiveNumbers.jl

*Next-generation matrices and reproduction numbers for compartmental epidemic models
written with [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) or
[Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/), computed symbolically where a
closed form exists and numerically otherwise.*

The construction follows Diekmann, Heesterbeek and Roberts (2010), *The construction of
next-generation matrices for compartmental epidemic models*, J. R. Soc. Interface
7:873–885 ([doi:10.1098/rsif.2009.0386](https://doi.org/10.1098/rsif.2009.0386)).
See [Mathematical background](@ref) for the recipe and [Tutorial](@ref) for a worked
introduction.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/epirecipes/ReproductiveNumbers.jl")
```

## A first example

```@example index
using ReproductiveNumbers, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters β σ γ μ N
@variables S(t) E(t) I(t) R(t)
eqs = [D(S) ~ μ * N - β * S * I / N - μ * S,
    D(E) ~ β * S * I / N - (σ + μ) * E,
    D(I) ~ σ * E - (γ + μ) * I,
    D(R) ~ γ * I - μ * R]
seir = complete(System(eqs, t; name = :seir))

ngm = next_generation_matrix(seir, [E, I])
```

```@example index
basic_reproduction_number(ngm)
```

```@example index
basic_reproduction_number(ngm, Dict(β => 0.5, σ => 0.25, γ => 0.2, μ => 0.01, N => 1e3))
```

## What the package provides

  - [`next_generation_matrix`](@ref) splits the linearised infected subsystem into a
    transmission matrix `T` and a transition matrix `Σ`, and builds the next-generation
    matrix with large domain `K_L = -TΣ⁻¹` and its restriction `K` to the
    states-at-infection. The infection-free steady state is found automatically when it is
    unique and can be supplied otherwise. The classification of terms (or, for Catalyst,
    reactions) as transmissions can be overridden, see
    [Choosing what counts as a transmission](@ref).
  - [`basic_reproduction_number`](@ref) returns `R₀` in closed form when the
    next-generation matrix decomposes into irreducible blocks that are `1 × 1`, rank one or
    `2 × 2`, and otherwise throws a [`NoClosedFormError`](@ref); with parameter values it
    always returns a number.
  - [`type_reproduction_number`](@ref), [`small_domain_matrix`](@ref),
    [`characteristic_polynomial`](@ref), [`spectral_radius`](@ref),
    [`validate_decomposition`](@ref) and [`evaluate`](@ref) for further analysis.

## Where to go next

  - The [Vignettes](@ref) walk from the SIR model to type reproduction numbers, each a
    self-contained Quarto document.
  - The linear algebra behind the closed forms is machine-checked in Lean 4, see
    [Formal proofs](@ref).
  - The [API reference](@ref) documents every exported function and the internal helpers.

## Citing

Please cite the underlying method (Diekmann, Heesterbeek & Roberts 2010) and, if the
software was useful, the repository (`CITATION.cff`).
