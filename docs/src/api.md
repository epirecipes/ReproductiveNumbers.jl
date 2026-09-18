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
states_at_infection
infected_states
small_domain_matrix
evaluate
validate_decomposition
```

## Reproduction numbers

```@docs
basic_reproduction_number
type_reproduction_number
NoClosedFormError
```

## The infected subsystem

```@docs
infected_subsystem
disease_free_equilibrium
default_is_transmission
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

## Symbolic helpers

```@docs
symbolic_isequal
symbolic_iszero
tidy
additive_terms
symbolic_variables
depends_on
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
