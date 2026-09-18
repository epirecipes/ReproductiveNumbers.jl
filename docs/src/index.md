# ReproductiveNumbers.jl

*Next-generation matrices, the basic reproduction number `R₀` and the effective
reproduction number `R_t` for compartmental epidemic models written with
[ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) or
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

The effective reproduction number `R_t` is the same construction linearised at the
current state instead of the infection-free steady state, so it is a function of the
uninfected compartments:

```@example index
effective_reproduction_number(seir, [E, I])
```

It can be evaluated along a simulated trajectory, or added to the system as an observed
variable `Rt(t)` that the solver returns directly. Here `R_t` starts at `R₀` and falls
below one as susceptibles are depleted, which is when the epidemic peaks:

```@example index
using OrdinaryDiffEqTsit5, Plots
seir_Rt, Rt = add_effective_reproduction_number(seir, [E, I])
prob = ODEProblem(seir_Rt, [S => 999.0, E => 0.0, I => 1.0, R => 0.0,
                            β => 0.5, σ => 0.25, γ => 0.2, μ => 0.01, N => 1e3], (0.0, 120.0))
sol = solve(prob, Tsit5(); saveat = 1.0)
plot(sol; idxs = [Rt], label = "R_t", xlabel = "time (days)", ylabel = "R_t", linewidth = 2)
hline!([1.0]; label = "threshold", linestyle = :dash, color = :black)
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
    `2 × 2` (or when the small-domain matrix does), and otherwise throws a
    [`NoClosedFormError`](@ref); with parameter values it always returns a number.
  - [`effective_reproduction_number`](@ref) returns `R_t` as a closed-form function of the
    state (also with time-varying rates), evaluates it at a state or at a time of a
    solution, or along a whole trajectory; [`add_effective_reproduction_number`](@ref)
    adds `Rt(t)` to a system as an observed variable.
  - [`abbreviate`](@ref) rewrites the matrices in terms of named sojourn times, transition
    probabilities and transmission rates, or of your own definitions, so that they read the
    way the paper derives them.
  - [`type_reproduction_number`](@ref) for targeted control, [`sensitivities`](@ref) and
    [`elasticities`](@ref) of `R₀`, [`perron_vectors`](@ref), [`mean_sojourn_times`](@ref),
    [`small_domain_matrix`](@ref), [`characteristic_polynomial`](@ref),
    [`suggest_infected`](@ref), [`validate_decomposition`](@ref) and [`evaluate`](@ref)
    for further analysis, and a ForwardDiff-based constructor for models written as plain
    Julia functions.

## Where to go next

  - The [Vignettes](@ref) walk from the SIR model to type reproduction numbers, each a
    self-contained Quarto document.
  - The linear algebra behind the closed forms is machine-checked in Lean 4, see
    [Formal proofs](@ref).
  - The [API reference](@ref) documents every exported function and the internal helpers.

## Citing

Please cite the underlying method (Diekmann, Heesterbeek & Roberts 2010) and, if the
software was useful, the repository (`CITATION.cff`).
