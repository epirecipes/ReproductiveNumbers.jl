"""
    ReproductiveNumbers

Next-generation matrices and reproduction numbers for compartmental epidemic models
defined with [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) or
[Catalyst.jl](https://github.com/SciML/Catalyst.jl).

The construction follows Diekmann, Heesterbeek and Roberts (2010), *The construction of
next-generation matrices for compartmental epidemic models*, J. R. Soc. Interface
7:873–885 ([https://doi.org/10.1098/rsif.2009.0386](https://doi.org/10.1098/rsif.2009.0386)). The linearised infected subsystem
`ẋ = (T + Σ) x` is split into a transmission part `T` (new infections) and a transition
part `Σ` (progression, recovery, death). The next-generation matrix with large domain is
`K_L = -T Σ⁻¹`, the next-generation matrix `K = -Eᵀ T Σ⁻¹ E` restricts `K_L` to the
states-at-infection, and `R₀` is the spectral radius of either.

See [`next_generation_matrix`](@ref) and [`basic_reproduction_number`](@ref).
"""
module ReproductiveNumbers

using DocStringExtensions
using LinearAlgebra
using Random: Xoshiro
using ForwardDiff
using SciMLBase: AbstractTimeseriesSolution
using Symbolics
using Symbolics: Num, unwrap, wrap, get_variables, jacobian, substitute
using SymbolicUtils
using SymbolicUtils: iscall, operation, arguments
using SymbolicUtils.Rewriters: Postwalk
using ModelingToolkit
using ModelingToolkit: AbstractSystem, full_equations, equations, unknowns, parameters,
                       get_iv,
                       isdifferential, mtkcompile
using SymbolicIndexingInterface: is_parameter, is_variable, is_independent_variable, getp,
                                 getu

export NextGenerationMatrix, NoClosedFormError
export next_generation_matrix, basic_reproduction_number, type_reproduction_number
export effective_reproduction_number, add_effective_reproduction_number
export disease_free_equilibrium, transmission_transition_matrices, infected_subsystem
export suggest_infected, mean_sojourn_times
export spectral_radius, characteristic_polynomial, irreducible_blocks, perron_vectors
export sensitivities, elasticities
export states_at_infection, infected_states, small_domain_matrix, evaluate,
       validate_decomposition, transmission_method
export symbolic_isequal, reaction_is_transmission

include("utils.jl")
include("routing.jl")
include("types.jl")
include("decomposition.jl")
include("ngm.jl")
include("spectral.jl")
include("r0.jl")
include("analysis.jl")

end
