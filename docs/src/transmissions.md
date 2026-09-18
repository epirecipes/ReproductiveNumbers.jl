# Choosing what counts as a transmission

The split of the linearised infected subsystem into transmissions `T` and transitions `Σ`
is a modelling decision. The package has a default for each kind of model and several
ways to override it.

## Defaults

**ModelingToolkit systems.** Each infected equation is split into its top-level additive
terms (products are not expanded, so `β (N − I) I / N` is one term). A term is a
transmission if it depends on an uninfected state (`β S I / N`), or if it is non-linear in
the infected states (`β (N − I) I / N` after `S` has been eliminated). Terms that are linear
in the infected states alone (`σ E`, `−γ I`, `p μ I`) are transitions. See
[`ReproductiveNumbers.default_is_transmission`](@ref).

**Catalyst networks.** A reaction is a transmission if its net stoichiometry increases the
total number of individuals in the infected species. `S + I → E + I` and `I → I + J`
(vertical transmission) are transmissions; `E → I` and `I → R` are not.

The defaults reproduce the choices made in Diekmann et al. (2010) for their examples, with
one important exception: vertical transmission written as an ODE term (`p μ I` in the
equation for infected juveniles) is linear in the infected states and involves no
susceptible class, so the ODE heuristic files it under transitions. Write such models as
reaction networks, or override.

## Overriding for ModelingToolkit systems

The `transmission` keyword of [`next_generation_matrix`](@ref) accepts:

  - a `Function` called on each additive term (a `Num`) returning `Bool`;
  - a `Vector` giving the new-infection rate `Fᵢ` of every infected compartment explicitly
    (van den Driessche & Watmough's `𝓕`); the transitions are the remainder.

```@example transmissions
using ReproductiveNumbers, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters β μ ν γ p N
@variables S(t) J(t) I(t)
eqs = [D(S) ~ μ * N - p * μ * I - β * S * I / N - μ * S,
    D(J) ~ p * μ * I - (ν + μ) * J,
    D(I) ~ β * S * I / N + ν * J - (γ + μ) * I]
vertical = complete(System(eqs, t; name = :vertical))

# default: births of infected juveniles are transitions, only I is a state-at-infection
states_at_infection(next_generation_matrix(vertical, [J, I]))
```

```@example transmissions
ngm = next_generation_matrix(vertical, [J, I];
    transmission = [p * μ * I, β * S * I / N])
states_at_infection(ngm), ngm.T
```

```@example transmissions
basic_reproduction_number(ngm)
```

## Overriding for Catalyst networks

For a `ReactionSystem`, `transmission` may be a predicate on `Reaction`s, a
`Vector{Bool}` with one flag per reaction, or a vector of reaction indices.

```@example transmissions
using Catalyst
malaria = @reaction_network malaria begin
    @parameters a b c γ μᵥ N_H N_V
    a * b / N_H, S_H + I_V --> I_H + I_V
    γ, I_H --> S_H
    a * c / N_H, S_V + I_H --> I_V + I_H
    μᵥ * N_V, 0 --> S_V
    μᵥ, (S_V, I_V) --> 0
end
@unpack S_H, I_H, S_V, I_V, N_H, N_V = malaria
eq = Dict(S_H => N_H, S_V => N_V)
R0_two_generations = basic_reproduction_number(malaria, [I_H, I_V]; equilibrium = eq)
```

```@example transmissions
R0_host = basic_reproduction_number(
    malaria, [I_H, I_V]; equilibrium = eq, transmission = [1])
```

Both are legitimate: they agree on whether `R₀ > 1` but the host-to-host version is the
square of the two-generation version.

## Sanity checks

Whatever the choice, `T` must be non-negative and `Σ` must have non-negative off-diagonal
and non-positive diagonal entries. Terms that *remove* individuals from an infected
compartment but involve an uninfected state (density-dependent death `-d (S + I) I`, say)
are misfiled by the default rule; the package warns when a transmission term or an entry
of `T` is manifestly negative:

```@example transmissions
@parameters d
dd = complete(System([D(S) ~ -β * S * I / N,
                      D(I) ~ β * S * I / N - γ * I - d * (S + I) * I], t; name = :dd))
ngm_dd = next_generation_matrix(dd, [I]; equilibrium = Dict(S => N),
                                transmission = [β * S * I / N])
basic_reproduction_number(ngm_dd)
```

[`validate_decomposition`](@ref) checks the sign conventions at given parameter values;
the regression tests use it together with the threshold property. Finally,
[`suggest_infected`](@ref) proposes the infected compartments when you are unsure:

```@example transmissions
suggest_infected(dd)
```
