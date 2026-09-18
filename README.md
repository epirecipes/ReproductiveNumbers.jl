# ReproductiveNumbers.jl

[![Stable docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://epirecipes.github.io/ReproductiveNumbers.jl/stable/)
[![Dev docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://epirecipes.github.io/ReproductiveNumbers.jl/dev/)
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
    `NextGenerationMatrix` holding `T`, `Σ`, `K_L`, `E`, `K`, the states-at-infection and a
    record of how transmissions were identified (`transmission_method`). The equilibrium is
    computed automatically when it is unique (linear, or polynomial with Nemo loaded) and
    must be given otherwise (`equilibrium = Dict(S => N)`); any state can be supplied to get
    a reproduction number at that state. Transmissions are identified by named strategies
    (`:auto`, `:uninfected_dependence`, `:nonlinear_in_infected`, `:stoichiometry` for
    Catalyst), a predicate, or an explicit vector of new-infection rates; misfiled loss terms
    trigger a warning. `suggest_infected(sys)` proposes the infected compartments.
  - `basic_reproduction_number(ngm)` returns a closed form when one exists: `K` is reduced
    to irreducible blocks; each block is handled if it is `1 × 1`, rank one (spectral radius
    = trace) or `2 × 2` (quadratic formula); the result is the symbolic maximum over blocks,
    and the small-domain matrix `K_S` is tried when `K` has none. Otherwise a
    `NoClosedFormError` is thrown and `basic_reproduction_number(ngm, p)` gives the value
    for parameter values `p` (a `Dict`, pairs, or an `ODEProblem`).
  - `effective_reproduction_number(sys, infected)` gives the effective reproduction number
    `R_t` in closed form as a function of the state (`R₀ S(t)/N` for the SIR model),
    evaluates it at a state or a time of a solution, or along a whole trajectory; and
    `add_effective_reproduction_number(sys, infected)` adds `Rt(t)` to the system as an
    observed variable so that `sol[Rt]` comes straight from the solver;
    `type_reproduction_number(ngm, types)` for targeted control (Roberts & Heesterbeek 2003),
    `sensitivities` and `elasticities` of `R₀`, `perron_vectors` (age/type distribution of
    new infections and reproductive values), `mean_sojourn_times`, `small_domain_matrix`,
    `characteristic_polynomial`, `spectral_radius`, `irreducible_blocks`,
    `validate_decomposition`, `evaluate`, `disease_free_equilibrium`.
  - Models written as plain Julia functions: `next_generation_matrix(F, V, x₀, p)` differentiates
    van den Driessche & Watmough's `𝓕` and `𝒱` with ForwardDiff.
  - With Latexify loaded, `latexify(ngm)` renders the decomposition.

## Documentation and vignettes

The documentation site ([https://epirecipes.github.io/ReproductiveNumbers.jl/dev/](https://epirecipes.github.io/ReproductiveNumbers.jl/dev/)) has a
mathematical background page, a tutorial, a guide to choosing what counts as a
transmission, the API reference, the rendered vignettes and the rendered Lean proofs. It
is built by the `Documentation` workflow on every push to `main` (and previews for pull
requests) with Documenter.jl and deployed to the `gh-pages` branch.

Fifteen Quarto vignettes (Julia engine), from simple to complex, live in
[`vignettes/`](vignettes/): SIR, SEIR with demography, two latent stages (rank one),
vector-borne transmission and the choice of generation, Catalyst reaction networks,
heterogeneous mixing, type reproduction numbers, vertical and sexual transmission, bovine
viral diarrhoea, staged progression, multiple strains, waning immunity checked against
steady states, within-host dynamics, an age-structured contact-matrix model, and models as
plain functions with leaky versus all-or-nothing vaccines.

To build the site locally:

```sh
(cd vignettes && julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()' && quarto render)
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl     # site in docs/build/
```

**One-time setup after the first push:** the `Documentation` workflow deploys with the
repository's `GITHUB_TOKEN`; if GitHub Pages is not enabled automatically, set *Settings →
Pages → Source* to *Deploy from a branch*, branch `gh-pages`, folder `/ (root)`.

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

Unit tests cover the symbolic helpers, spectral radius, decomposition and analysis
functions; integration tests run the whole pipeline on ModelingToolkit and Catalyst models;
regression tests reproduce closed forms from the literature (Diekmann et al. sections 2.1,
2.2, 4.1, 4.2, 4.3; van den Driessche & Watmough's treatment and staged-progression
models), check the threshold property against the Jacobian and against simulated
epidemics, and check `R₀` against the steady states of the full non-linear model computed
with SteadyStateDiffEq and NonlinearSolve over parameter sweeps.

## How it is implemented: Symbolics.jl and automatic differentiation

The package is a thin layer of linear algebra on top of the SciML symbolic stack.

**Getting the equations.** A ModelingToolkit `System` is queried for its unknowns and
right-hand sides (`equations`, or `full_equations` after `mtkcompile` when there are
algebraic or observed variables). A Catalyst `ReactionSystem` is converted with
`Catalyst.ode_model`, but the split into transmissions and transitions is done at the level
of reactions, using `netstoichmat` and `oderatelaw`, before any ODE is formed.

**Splitting terms.** Each infected right-hand side is broken into its top-level additive
terms with SymbolicUtils (`isadd`, `arguments`); products are deliberately not expanded, so
a term such as `β (N − I) I / N` is classified as a whole. Whether a term depends on an
uninfected state uses `Symbolics.get_variables`; whether it is non-linear in the infected
states uses a symbolic Hessian, `Symbolics.jacobian` applied twice.

**Linearisation.** `T` and `Σ` are symbolic Jacobians, `Symbolics.jacobian(F, x)` and
`Symbolics.jacobian(G, x)`, evaluated at the infection-free steady state with
`Symbolics.substitute`. The steady state itself comes from `Symbolics.symbolic_linear_solve`
on the uninfected subsystem, or from `Symbolics.symbolic_solve` (which needs Nemo) when
that subsystem is polynomial.

**Symbolic matrix algebra.** `K_L = -TΣ⁻¹` uses Symbolics' `inv` for `Matrix{Num}`;
`K`, `K_S`, `-Σ⁻¹`, the type-reproduction-number matrix and the characteristic polynomial
are ordinary matrix products, inverses and determinants over `Num`. Results are tidied by
walking the expression tree with `SymbolicUtils.Rewriters.Postwalk` to cancel paired minus
signs (Symbolics canonicalises `−(γ + μ)` as `−γ − μ`, so the sign has to be recovered from
products and quotients) and by `Symbolics.simplify_fractions` to cancel common factors.

**Deciding whether an expression is zero.** The block decomposition, the rank test and the
de-duplication of block spectral radii all hinge on recognising identically zero expressions.
A candidate is first evaluated at fixed pseudo-random points (a non-zero value is a proof
of non-vanishing); if it is a rational function of the parameters, `simplify_fractions`
followed by `Symbolics.expand` of the numerator decides exactly; otherwise (square roots,
`max`) agreement at the random points is accepted as overwhelming evidence. This is the only
non-formal step in the symbolic pipeline and can be disabled.

**Closed forms.** The dominant eigenvalue is assembled from the irreducible blocks of `K`
(Tarjan's algorithm on the pattern of structurally non-zero entries) using formulae that
are proved in the Lean files: `1 × 1` entries, the trace for rank-one blocks, and
`(a + d + √((a − d)² + 4bc))/2` for `2 × 2` blocks, joined by a symbolic `max`. Sensitivities
and elasticities of the result are `Symbolics.derivative` with respect to each parameter.

**Numbers.** Numeric evaluation substitutes parameter values (`Symbolics.substitute`),
reading them from a `Dict`, from pairs, or through SymbolicIndexingInterface (`getp`,
`getu`) from an `ODEProblem` or a solution; calls that Symbolics leaves unevaluated
(`sqrt(6.0)`, `max(2.0, 3.0)`) are folded by a small recursive evaluator. The numeric path
then uses LinearAlgebra (`-T / Σ`, `eigvals`, `eigen`).

**Automatic differentiation.** Two kinds are used. The symbolic Jacobians above are
symbolic differentiation by Symbolics.jl, which is what makes closed forms possible. For
models written as plain Julia functions `F(x, p)` and `V(x, p)` rather than as symbolic
systems, `next_generation_matrix(F, V, x₀, p)` obtains `T` and `Σ` as `ForwardDiff.jacobian`
of `F` and `-V` at the infection-free state, so any model that can be evaluated with dual
numbers can be analysed numerically without being rewritten symbolically. The numeric
elasticities in the contact-matrix vignette use central finite differences of the numeric
eigenvalue, which is the pragmatic choice when no closed form exists.

**Extensions.** Catalyst and Latexify support live in package extensions
(`ext/`), so neither is loaded unless the user loads it.

## References

  - Diekmann O, Heesterbeek JAP, Roberts MG (2010). The construction of next-generation
    matrices for compartmental epidemic models. *J R Soc Interface* 7:873–885.
  - van den Driessche P, Watmough J (2002). Reproduction numbers and sub-threshold endemic
    equilibria for compartmental models of disease transmission. *Math Biosci* 180:29–48.
  - Roberts MG, Heesterbeek JAP (2003). A new method for estimating the effort required to
    control an infectious disease. *Proc R Soc Lond B* 270:1359–1364.
