"""
$(TYPEDEF)

Thrown when [`spectral_radius`](@ref) cannot produce a closed-form expression for the
dominant eigenvalue of a symbolic matrix. Evaluate numerically instead with
[`basic_reproduction_number`](@ref)`(ngm, parameter_values)`, or inspect
[`characteristic_polynomial`](@ref).
"""
struct NoClosedFormError <: Exception
    msg::String
end
Base.showerror(io::IO, e::NoClosedFormError) = print(io, "NoClosedFormError: ", e.msg)

"""
$(TYPEDEF)

The decomposition of the linearised infected subsystem of a compartmental model into
transmissions and transitions, together with the next-generation matrices built from it
(Diekmann, Heesterbeek & Roberts 2010).

The matrices are symbolic (`Matrix{Num}`) when built from a model and numeric
(`Matrix{Float64}`) after [`evaluate`](@ref). The fields `F`, `G` and `method` record how
the split into transmissions and transitions was made, because that choice determines the
value of `R₀`.

# Fields
$(TYPEDFIELDS)
"""
struct NextGenerationMatrix{M <: AbstractMatrix, E <: AbstractMatrix{Int}}
    "Infected state variables; row/column order of `T`, `Σ` and `K_L`."
    infected::Vector{Num}
    "Uninfected state variables."
    uninfected::Vector{Num}
    "The state at which the system was linearised (normally the infection-free steady state)."
    equilibrium::Dict{Num, Any}
    "Transmission matrix: `T[i, j]` is the rate at which infected individuals in state `j` produce new infections in state `i`."
    T::M
    "Transition matrix: `Σ[i, j]` is the rate of all other transitions from `j` to `i`; `-Σ⁻¹` gives the expected times spent in each state."
    Σ::M
    "Next-generation matrix with large domain, `K_L = -T Σ⁻¹`."
    K_L::M
    "Selection matrix whose columns are the unit vectors of the states-at-infection."
    E::E
    "Next-generation matrix, `K = Eᵀ K_L E`, indexed by the states-at-infection."
    K::M
    "Which infected states are states-at-infection (non-zero rows of `T`)."
    states_at_infection::Vector{Num}
    "New-infection rate of each infected compartment before linearisation (van den Driessche & Watmough's `𝓕`); empty when built directly from matrices."
    F::Vector{Num}
    "Remaining (transition) rate of each infected compartment before linearisation; empty when built directly from matrices."
    G::Vector{Num}
    "How transmissions were identified: `:auto`, `:uninfected_dependence`, `:nonlinear_in_infected`, `:stoichiometry`, `:predicate`, `:explicit`, `:matrices` or `:functions`."
    method::Symbol
end

"""
$(TYPEDSIGNATURES)

The infected state variables of the decomposition.
"""
infected_states(ngm::NextGenerationMatrix) = ngm.infected

"""
$(TYPEDSIGNATURES)

The states-at-infection: infected states an individual can occupy immediately after being
infected, i.e. the non-zero rows of the transmission matrix `T`.
"""
states_at_infection(ngm::NextGenerationMatrix) = ngm.states_at_infection

"""
$(TYPEDSIGNATURES)

The transmission (`T`, van den Driessche & Watmough's `F`) and transition (`Σ`, van den
Driessche & Watmough's `-V`) matrices of the decomposition.
"""
transmission_transition_matrices(ngm::NextGenerationMatrix) = (ngm.T, ngm.Σ)

"""
$(TYPEDSIGNATURES)

The strategy that was used to classify terms or reactions as transmissions (see the
`transmission` keyword of [`next_generation_matrix`](@ref)), and the resulting
new-infection rates `F` and transition rates `G` of every infected compartment, as a
named tuple `(method, F, G)`. `F` and `G` are empty when the decomposition was built
directly from matrices or functions.
"""
transmission_method(ngm::NextGenerationMatrix) = (method = ngm.method, F = ngm.F, G = ngm.G)

"""
$(TYPEDSIGNATURES)

`true` if the matrices of `ngm` are numeric.
"""
isnumeric(ngm::NextGenerationMatrix) = eltype(ngm.K) <: Number && !(eltype(ngm.K) <: Num)

function Base.show(io::IO, ::MIME"text/plain", ngm::NextGenerationMatrix)
    kind = isnumeric(ngm) ? "numeric" : "symbolic"
    println(io, "NextGenerationMatrix ($kind) with ", length(ngm.infected),
        " infected state(s) and ", length(ngm.states_at_infection),
        " state(s)-at-infection")
    println(io, "  infected:            ", ngm.infected)
    println(io, "  states-at-infection: ", ngm.states_at_infection)
    println(io, "  transmissions identified by: ", repr(ngm.method))
    if !isempty(ngm.F)
        for (x, f) in zip(ngm.infected, ngm.F)
            println(io, "    new infections in ", x, ": ", f)
        end
    end
    println(io, "  T = ")
    Base.print_matrix(io, ngm.T, "      ")
    println(io)
    println(io, "  Σ = ")
    Base.print_matrix(io, ngm.Σ, "      ")
    println(io)
    println(io, "  K = ")
    Base.print_matrix(io, ngm.K, "      ")
end
