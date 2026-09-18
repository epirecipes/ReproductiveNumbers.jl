# Tutorial

This tutorial builds a next-generation matrix step by step, first from a ModelingToolkit
system and then from a Catalyst reaction network, and shows what to do when no closed
form exists. The [Vignettes](@ref) go further on each topic.

## A ModelingToolkit model

We use an SEIR model with births and deaths, so that the infection-free steady state is
unique.

```@example tutorial
using ReproductiveNumbers, ModelingToolkit, Symbolics, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters β σ γ μ N
@variables S(t) E(t) I(t) R(t)
eqs = [D(S) ~ μ * N - β * S * I / N - μ * S,
    D(E) ~ β * S * I / N - (σ + μ) * E,
    D(I) ~ σ * E - (γ + μ) * I,
    D(R) ~ γ * I - μ * R]
seir = complete(System(eqs, t; name = :seir))
nothing # hide
```

### Which states are infected?

You tell the package which unknowns are infected compartments, as symbolic variables or as
`Symbol`s. Their order fixes the rows and columns of `T`, `Σ` and `K_L`.

```@example tutorial
sub = infected_subsystem(seir, [E, I])
sub.f_infected
```

### The infection-free steady state

Setting the infected states to zero and solving the remaining equations:

```@example tutorial
disease_free_equilibrium(seir, [E, I])
```

If the steady state is not unique (an SIR model without demography, say), pass it with the
`equilibrium` keyword. Any state is accepted, which gives a reproduction number at that
state rather than the basic one.

### Transmissions, transitions and the next-generation matrices

```@example tutorial
ngm = next_generation_matrix(seir, [E, I])
```

The object holds everything:

```@example tutorial
T, Σ = transmission_transition_matrices(ngm)
T
```

```@example tutorial
Σ
```

```@example tutorial
ngm.K_L
```

Only `E` is a state-at-infection, so `K` is `1 × 1` and `R₀` is its single entry:

```@example tutorial
states_at_infection(ngm), ngm.K
```

```@example tutorial
R0 = basic_reproduction_number(ngm)
```

### Numbers

Give parameter values as a `Dict`, pairs, or an `ODEProblem` built from the same system:

```@example tutorial
p = Dict(β => 0.5, σ => 0.25, γ => 0.2, μ => 0.01, N => 1000.0)
basic_reproduction_number(ngm, p)
```

[`evaluate`](@ref) returns a numeric [`NextGenerationMatrix`](@ref) with the same fields:

```@example tutorial
num = evaluate(ngm, p)
num.K_L
```

```@example tutorial
validate_decomposition(num)
```

`R₀ > 1` exactly when the Jacobian `T + Σ` has an eigenvalue with positive real part:

```@example tutorial
(basic_reproduction_number(num) > 1, maximum(real, eigvals(num.T + num.Σ)) > 0)
```

## A Catalyst model

With a reaction network the classification is read from the stoichiometry: a reaction is a
transmission if it increases the number of infected individuals.

```@example tutorial
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

## When there is no closed form

Three or more interacting types generally have no closed-form dominant eigenvalue.

```@example tutorial
@parameters q γ₃ C[1:3, 1:3] Nᵢ[1:3]
@variables Sᵢ(t)[1:3] Iᵢ(t)[1:3]
Cm, Nv, Sv, Iv = collect(C), collect(Nᵢ), collect(Sᵢ), collect(Iᵢ)
λ = [q * sum(Cm[i, j] * Iv[j] / Nv[j] for j in 1:3) for i in 1:3]
age = complete(System(
    [[D(Sv[i]) ~ -λ[i] * Sv[i] for i in 1:3];
     [D(Iv[i]) ~ λ[i] * Sv[i] - γ₃ * Iv[i] for i in 1:3]],
    t;
    name = :age))
ngm_age = next_generation_matrix(age, Iv; equilibrium = Dict(Sv[i] => Nv[i] for i in 1:3))
try
    basic_reproduction_number(ngm_age)
catch err
    err
end
```

```@example tutorial
contacts = [10.0 3.0 1.0; 3.0 8.0 2.0; 1.0 2.0 4.0]
vals = Dict{Any, Any}(Cm[i, j] => contacts[i, j] for i in 1:3, j in 1:3)
for i in 1:3
    vals[Nv[i]] = 100.0
end
vals[q] = 0.02
vals[γ₃] = 0.2
basic_reproduction_number(ngm_age, vals)
```

## Building from `T` and `Σ` directly

If you already have the matrices (from a paper, say), skip the model:

```@example tutorial
@parameters b s g
ngm_direct = next_generation_matrix([0 b; 0 0], [-s 0; s -g]; infected = [:E, :I])
basic_reproduction_number(ngm_direct)
```
