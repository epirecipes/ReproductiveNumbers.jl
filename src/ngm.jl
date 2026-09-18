# Building next-generation matrices.

"""
$(TYPEDSIGNATURES)

Given `T` and `Σ`, build the next-generation matrix with large domain `K_L = -T Σ⁻¹`, the
selection matrix `E` of states-at-infection (non-zero rows of `T`) and the next-generation
matrix `K = Eᵀ K_L E`.
"""
function assemble(infected, uninfected, equilibrium, T::AbstractMatrix, Σ::AbstractMatrix)
    n = length(infected)
    size(T) == (n, n) && size(Σ) == (n, n) ||
        throw(DimensionMismatch("T and Σ must be $(n)×$(n)"))
    K_L = _large_domain(T, Σ)
    rows = [i for i in 1:n if any(j -> !_structural_zero(T[i, j]), 1:n)]
    E = zeros(Int, n, length(rows))
    for (k, i) in enumerate(rows)
        E[i, k] = 1
    end
    K = _restrict(K_L, E)
    return NextGenerationMatrix(Vector{Num}(infected), Vector{Num}(uninfected),
        Dict{Num, Any}(equilibrium), T, Σ, K_L, E, K, infected[rows])
end

_structural_zero(x::Num) = symbolic_iszero(x; numeric = false)
_structural_zero(x::Number) = iszero(x)

function _large_domain(T::AbstractMatrix{Num}, Σ::AbstractMatrix{Num})
    Σinv = try
        inv(Σ)
    catch err
        throw(ArgumentError("the transition matrix Σ could not be inverted symbolically " *
                            "($(sprint(showerror, err))); is every infected state eventually left?"))
    end
    return tidy(-T * Σinv)
end
function _large_domain(T::AbstractMatrix{<:Real}, Σ::AbstractMatrix{<:Real})
    return -T / Σ
end

_restrict(K_L::AbstractMatrix{Num}, E) = tidy(Num.(E' * K_L * E); fractions = false)
_restrict(K_L::AbstractMatrix{<:Real}, E) = E' * K_L * E

"""
    next_generation_matrix(sys, infected; equilibrium = nothing, transmission = nothing)

Construct the [`NextGenerationMatrix`](@ref) of the ModelingToolkit system `sys`
(a `System` built from ordinary differential equations, or a Catalyst `ReactionSystem`
when Catalyst is loaded).

# Arguments

  - `infected`: the infected state variables, as symbolic variables (`[E, I]`, `[sys.E, sys.I]`)
    or as `Symbol`s (`[:E, :I]`). Their order fixes the row and column order of `T`, `Σ`
    and `K_L`.
  - `equilibrium`: the state at which to linearise, as a `Dict`/pairs mapping unknowns to
    values or expressions (`Dict(S => N)`). Infected states default to zero. When
    omitted, [`disease_free_equilibrium`](@ref) is used. Any other state may be given to
    obtain a reproduction number at that state (for instance an effective reproduction
    number with partial immunity).
  - `transmission`: how to decide which terms of the infected equations are new
    infections. `nothing` uses [`default_is_transmission`](@ref) (terms involving an
    uninfected state, or non-linear in the infected states). A `Function` is called on
    each additive term (a `Num`) and must return `Bool`. A `Vector` of expressions gives
    the new-infection rate `F_i` of each infected compartment explicitly (van den
    Driessche & Watmough's `𝓕`); the transitions are then the remainder of each equation.
    For Catalyst models, a `Function` is instead called on each `Reaction`, and a
    `Vector{Bool}` or vector of reaction indices marks the transmission reactions; by
    default a reaction is a transmission if it increases the total number of infected
    individuals.

Different choices of what counts as a transmission give different next-generation matrices
and different values of `R₀`, but they all share the threshold property `R₀ > 1` if and
only if the infection-free steady state is unstable (Diekmann et al. 2010, section 2).
"""
function next_generation_matrix(sys::AbstractSystem, infected;
        equilibrium = nothing, transmission = nothing)
    sub = infected_subsystem(sys, infected)
    x, y = sub.infected, sub.uninfected
    eq = _equilibrium(sys, infected, equilibrium, x, y)
    F, G = _split(sub.f_infected, x, y, transmission)
    T, Σ = linearise(F, G, x, eq)
    return assemble(x, y, eq, T, Σ)
end

function _equilibrium(sys, infected, equilibrium, x, y)
    if equilibrium === nothing
        return disease_free_equilibrium(sys, infected)
    end
    eq = Dict{Num, Any}()
    for v in x
        eq[v] = 0
    end
    for (k, v) in substitution_map(equilibrium)
        i = _findsym(k, x)
        j = _findsym(k, y)
        if i !== nothing
            eq[x[i]] = v
        elseif j !== nothing
            eq[y[j]] = v
        else
            throw(ArgumentError("`$(k)` in `equilibrium` is not an unknown of the system"))
        end
    end
    return eq
end

function _split(f, x, y, transmission)
    if transmission === nothing
        return split_terms(f, term -> default_is_transmission(term, x, y))
    elseif transmission isa Function
        return split_terms(f, transmission)
    elseif transmission isa AbstractVector
        length(transmission) == length(f) ||
            throw(ArgumentError("`transmission` must give one new-infection rate per infected state ($(length(f)))"))
        F = Num.(transmission)
        return F, f .- F
    else
        throw(ArgumentError("`transmission` must be `nothing`, a predicate on terms, or a vector of new-infection rates"))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a [`NextGenerationMatrix`](@ref) directly from a transmission matrix `T` and a
transition matrix `Σ` (or, equivalently, van den Driessche & Watmough's `F` and `V` with
`Σ = -V`). `infected` names the rows/columns, as symbolic variables or `Symbol`s.
"""
function next_generation_matrix(T::AbstractMatrix, Σ::AbstractMatrix;
        infected = [Symbolics.variable(:x, i) for i in 1:size(T, 1)])
    x = Num[v isa Symbol ? Symbolics.variable(v) : Num(v) for v in infected]
    symbolic = !(_isnumeric(T) && _isnumeric(Σ))
    return assemble(x, Num[], Dict{Num, Any}(v => 0 for v in x), _promote(T, symbolic),
        _promote(Σ, symbolic))
end
# `Num <: Real`, so a symbolic matrix must be recognised before the numeric fallback.
function _isnumeric(A::AbstractMatrix)
    eltype(A) <: Real && !(eltype(A) <: Num) && !any(x -> x isa Num, A)
end
_promote(A, symbolic::Bool) = symbolic ? Matrix{Num}(Num.(A)) : Matrix{Float64}(A)

"""
$(TYPEDSIGNATURES)

The next-generation matrix with small domain `K_S = -R Σ⁻¹ C` obtained from a rank
factorisation `T = C R` (Diekmann et al. 2010, section 2.3). When `T` has rank one, so
that all states-at-infection are entered in fixed proportions, `K_S` is a scalar equal to
`R₀`. Only the rank-one factorisation is attempted; if `T` is not rank one, `K` itself is
returned.
"""
function small_domain_matrix(ngm::NextGenerationMatrix)
    T = ngm.T
    n = size(T, 1)
    nz = [i for i in 1:n if !all(_structural_zero, T[i, :])]
    isempty(nz) && return ngm.K
    Rrow = T[first(nz), :]
    j = findfirst(x -> !_structural_zero(x), Rrow)
    Cvec = [zero(eltype(T)) for _ in 1:n]
    for i in nz
        c = T[i, j] / Rrow[j]
        all(k -> _iszero_like(T[i, k] - c * Rrow[k]), 1:n) || return ngm.K
        Cvec[i] = c
    end
    return _small(reshape(Rrow, 1, n), ngm.Σ, reshape(Cvec, n, 1))
end
_iszero_like(x::Num) = symbolic_iszero(x)
_iszero_like(x::Number) = isapprox(x, 0; atol = 1e-12)
_small(Rm::AbstractMatrix{Num}, Σ, C) = tidy(-Rm * inv(Σ) * C)
_small(Rm::AbstractMatrix{<:Real}, Σ, C) = -(Rm / Σ) * C

"""
    reaction_is_transmission(k, infected_indices, netstoich)

Default classification of reaction `k` of a Catalyst `ReactionSystem` (available when
Catalyst is loaded): `true` if the net stoichiometry of the reaction (column `k` of
`netstoich`, rows indexed by species) increases the total number of individuals in the
infected species `infected_indices`. `S + I --> E + I` and `I --> I + J` are
transmissions, `E --> I` and `I --> R` are transitions.
"""
function reaction_is_transmission end
