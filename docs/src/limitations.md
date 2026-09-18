# Limitations and failure modes

This page collects the situations in which the package cannot do what you might expect,
and what to do instead.

## Which terms are transmissions

The default term classification (`:auto`) treats a term as a transmission if it involves an
uninfected state or is non-linear in the infected states. Two kinds of model defeat it:

- **Vertical transmission written as an ODE term** (`p μ I` feeding infected juveniles) is
  linear in the infected states and involves no susceptible class, so it is filed under
  transitions. Write the model as a reaction network (the stoichiometry rule recognises
  `I --> I + J`), or give the transmission terms explicitly.
- **Loss terms that involve uninfected states**, such as density-dependent death
  `-d (S + I) I` or predation on infecteds, are classified as transmissions and would give
  a negative entry of `T`. The package warns (`warn = true`) when a transmission term or
  an entry of `T` is manifestly negative; reassign such terms with an explicit
  `transmission` vector or a predicate.

Whatever strategy is used is recorded in the result and printed with it
([`transmission_method`](@ref)). For numeric values, [`validate_decomposition`](@ref)
checks the sign conventions.

## Infection-free steady states

[`disease_free_equilibrium`](@ref) solves the uninfected subsystem at `x = 0`. It handles
linear subsystems directly, and polynomial ones (logistic host growth, say) through
`Symbolics.symbolic_solve` when Nemo is loaded (and Groebner for several coupled unknowns);
among polynomial solutions it accepts exactly one candidate in which no uninfected state
vanishes and otherwise lists the candidates. Steady states that are not unique (an SIR
model without demography, for which every `S` is a steady state) or not polynomial must
be supplied with the `equilibrium` keyword.

## Closed forms

`R₀` is returned in closed form only when every irreducible block of `K` is `1 × 1`, rank
one, or `2 × 2`, or when the small-domain matrix `K_S` has that property. General
`3 × 3` and larger irreducible blocks have no closed-form dominant eigenvalue, and the
symbolic cubic/quartic solver of Symbolics cannot be applied to matrices with symbolic
entries in the current release, so a [`NoClosedFormError`](@ref) is thrown; use parameter
values, or [`characteristic_polynomial`](@ref).

## The rank test is partly probabilistic

Whether a block has rank one is decided from its `2 × 2` minors. A minor that is a rational
function of the parameters is simplified exactly. Minors involving other functions
(square roots from a nested closed form, say) are tested at random parameter values with a
fixed seed: agreement at four random points is overwhelming, but not formal, evidence of
identity. Pass `numeric_rank_check = false` to [`basic_reproduction_number`](@ref) or
[`spectral_radius`](@ref) to disable the probabilistic step (some rank-one matrices are
then reported as having no closed form).

## Models the package does not handle

- Non-autonomous systems (explicit dependence on time) are rejected; seasonally forced
  models need a different definition of `R₀` (Bacaër and Guernaoui 2006).
- Delay, stochastic, spatial and discrete-time models are out of scope.
- Reaction networks with coupled non-reaction equations, constant or boundary species, or
  hybrid noise/jump components fall back to the ODE route (with a warning) or fail to
  convert; ModelingToolkit parameter dependencies are not resolved before linearisation.
- Composed ModelingToolkit systems with several unknowns of the same name must be referred
  to by their symbolic variables, not by `Symbol`.

## Numerical caveats

Numeric next-generation matrices use dense linear algebra (`-T / Σ`, `eigvals`), which is
adequate for compartmental models with up to a few hundred compartments. The type
reproduction number is undefined, and an error is raised, when the infection can persist
among the types outside the chosen set (`ρ((I - P)K) ≥ 1`).
