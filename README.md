# ReproductiveNumbers.jl

[![CI](https://github.com/epirecipes/ReproductiveNumbers.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/epirecipes/ReproductiveNumbers.jl/actions/workflows/CI.yml)
[![Lean proofs](https://github.com/epirecipes/ReproductiveNumbers.jl/actions/workflows/Proofs.yml/badge.svg)](https://github.com/epirecipes/ReproductiveNumbers.jl/actions/workflows/Proofs.yml)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

Next-generation matrices and reproduction numbers for compartmental epidemic models
written with [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) or
[Catalyst.jl](https://github.com/SciML/Catalyst.jl), computed symbolically where a closed
form exists and numerically otherwise.

The construction follows Diekmann, Heesterbeek and Roberts (2010), *The construction of
next-generation matrices for compartmental epidemic models*, J. R. Soc. Interface
7:873–885 ([doi:10.1098/rsif.2009.0386](https://doi.org/10.1098/rsif.2009.0386)):

 1. take the infected subsystem and linearise it at the infection-free steady state;
 2. split the Jacobian into transmissions `T` (new infections) and transitions `Σ`
    (progression, recovery, death), so that `-Σ⁻¹` holds expected sojourn times;
 3. form the next-generation matrix with large domain `K_L = -T Σ⁻¹` and its restriction
    `K = Eᵀ K_L E` to the states-at-infection (the non-zero rows of `T`);
 4. `R₀` is the spectral radius of `K`, which equals that of `K_L`.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/epirecipes/ReproductiveNumbers.jl")
```

## Quick start

```julia
using ReproductiveNumbers, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters β σ γ μ N
@variables S(t) E(t) I(t) R(t)
eqs = [D(S) ~ μ * N - β * S * I / N - μ * S,
    D(E) ~ β * S * I / N - (σ + μ) * E,
    D(I) ~ σ * E - (γ + μ) * I,
    D(R) ~ γ * I - μ * R]
seir = complete(System(eqs, t; name = :seir))

ngm = next_generation_matrix(seir, [E, I])   # infection-free steady state found automatically
ngm.T, ngm.Σ, ngm.K_L, ngm.K                  # symbolic matrices
basic_reproduction_number(ngm)                 # (β*σ) / ((γ + μ)*(μ + σ))
basic_reproduction_number(ngm, Dict(β => 0.5, σ => 0.25, γ => 0.2, μ => 0.01, N => 1e3))
```

With Catalyst, transmission reactions are identified from the stoichiometry (a reaction is a
transmission if it increases the number of infected individuals):

```julia
using Catalyst
rn = @reaction_network seir begin
    @parameters β σ γ μ N
    β / N, S + I --> E + I
    σ, E --> I
    γ, I --> R
    μ * N, 0 --> S
    μ, (S, E, I, R) --> 0
end
basic_reproduction_number(rn, [:E, :I])
```

## What you get

  - `next_generation_matrix(sys, infected; equilibrium, transmission)` builds a
    `NextGenerationMatrix` holding `T`, `Σ`, `K_L`, `E`, `K` and the states-at-infection.
    The equilibrium is computed automatically when it is unique (models with demography) and
    must be given otherwise (`equilibrium = Dict(S => N)`); any state can be supplied to get
    a reproduction number at that state. The split into transmissions and transitions can be
    overridden with a predicate on terms (or reactions) or an explicit vector of
    new-infection rates, because that choice changes `R₀` but not its threshold.
  - `basic_reproduction_number(ngm)` returns a closed form when one exists: `K` is reduced
    to irreducible blocks; each block is handled if it is `1 × 1`, rank one (spectral radius
    = trace) or `2 × 2` (quadratic formula); the result is the symbolic maximum over blocks.
    Otherwise a `NoClosedFormError` is thrown and `basic_reproduction_number(ngm, p)` gives
    the value for parameter values `p` (a `Dict`, pairs, or an `ODEProblem`).
  - `type_reproduction_number(ngm, types)` for targeted control (Roberts & Heesterbeek 2003),
    `small_domain_matrix`, `characteristic_polynomial`, `spectral_radius`,
    `irreducible_blocks`, `validate_decomposition`, `evaluate`, `disease_free_equilibrium`.

## Vignettes

Quarto vignettes (Julia engine), from simple to complex, live in
[`vignettes/`](vignettes/): SIR, SEIR with demography, two latent stages (rank one),
vector-borne transmission and the choice of generation, Catalyst reaction networks,
heterogeneous mixing and age structure, and type reproduction numbers.

## Formal proofs

The matrix algebra the package relies on is formalised in Lean 4 with Mathlib under
[`proofs/`](proofs/): the restricted and large-domain matrices share their non-zero
spectrum, `det(T + Σ) = 0` iff `1` is an eigenvalue of `K_L`, rank-one matrices have their
trace as only non-zero eigenvalue, the `2 × 2` closed form and its threshold criterion, and
the worked examples. The proofs are rendered to HTML with mdgen and pandoc
(`proofs/scripts/check.sh`).

## Tests

```julia
using Pkg;
Pkg.test("ReproductiveNumbers");
```

Unit tests cover the symbolic helpers, spectral radius and decomposition; integration tests
run the whole pipeline on ModelingToolkit and Catalyst models; regression tests reproduce
closed forms from the literature and check the threshold property against the Jacobian and
against simulated epidemics.

## References

  - Diekmann O, Heesterbeek JAP, Roberts MG (2010). The construction of next-generation
    matrices for compartmental epidemic models. *J R Soc Interface* 7:873–885.
  - van den Driessche P, Watmough J (2002). Reproduction numbers and sub-threshold endemic
    equilibria for compartmental models of disease transmission. *Math Biosci* 180:29–48.
  - Roberts MG, Heesterbeek JAP (2003). A new method for estimating the effort required to
    control an infectious disease. *Proc R Soc Lond B* 270:1359–1364.
