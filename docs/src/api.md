# API reference

```@meta
CurrentModule = ReproductiveNumbers
```

```@docs
ReproductiveNumbers
```

## Building next-generation matrices

```@docs
next_generation_matrix
NextGenerationMatrix
transmission_transition_matrices
transmission_method
states_at_infection
infected_states
small_domain_matrix
mean_sojourn_times
evaluate
validate_decomposition
suggest_infected
```

## Reproduction numbers

```@docs
basic_reproduction_number
effective_reproduction_number
add_effective_reproduction_number
type_reproduction_number
NoClosedFormError
```

## Readable forms

```@docs
abbreviate
expand_definitions
```

## Sensitivity and structure of the next-generation matrix

```@docs
sensitivities
elasticities
perron_vectors
```

## The infected subsystem

```@docs
infected_subsystem
disease_free_equilibrium
default_is_transmission
depends_on_uninfected
nonlinear_in_infected
split_terms
linearise
ode_right_hand_sides
resolve_states
assemble
```

## Spectral radius

```@docs
spectral_radius
spectral_radius_2x2
irreducible_blocks
isrankone
characteristic_polynomial
```

## Routing

```@docs
symbolic_inverse
first_success
expression_size
TIDY_SIZE_LIMIT
LAPLACE_INVERSE_LIMIT
```

## Symbolic helpers

```@docs
symbolic_isequal
symbolic_iszero
tidy
additive_terms
symbolic_variables
depends_on
manifestly_negative
to_number
substitution_map
```

## Catalyst extension

Loading Catalyst adds methods of [`next_generation_matrix`](@ref),
[`infected_subsystem`](@ref) and [`disease_free_equilibrium`](@ref) for `ReactionSystem`s.
A reaction is classified as a transmission when its net stoichiometry increases the number
of infected individuals; see [Choosing what counts as a transmission](@ref) for the
overrides.

```@docs
reaction_is_transmission
```

## Latexify extension

Loading Latexify adds `latexify(ngm)` and a `text/latex` `show` method, so a
[`NextGenerationMatrix`](@ref) renders as an aligned display of `T`, `Σ`, `K` and `R₀`
in notebooks and Quarto documents.
