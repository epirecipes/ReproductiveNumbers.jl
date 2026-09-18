# Building next-generation matrices.

"""
$(TYPEDSIGNATURES)

Given `T` and `Σ`, build the next-generation matrix with large domain `K_L = -T Σ⁻¹`, the
selection matrix `E` of states-at-infection (non-zero rows of `T`) and the next-generation
matrix `K = Eᵀ K_L E`. `F`, `G` and `method` record how the split was made (see
[`NextGenerationMatrix`](@ref)).
"""
function assemble(infected, uninfected, equilibrium, T::AbstractMatrix, Σ::AbstractMatrix;
        F = Num[], G = Num[], method::Symbol = :matrices,
        definitions::AbstractDict = Dict{Num, Num}())
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
        Dict{Num, Any}(equilibrium), T, Σ, K_L, E, K, infected[rows],
        Vector{Num}(F), Vector{Num}(G), method, Dict{Num, Num}(definitions))
end

_structural_zero(x::Num) = symbolic_iszero(x; numeric = false)
_structural_zero(x::Number) = iszero(x)

function _large_domain(T::AbstractMatrix{Num}, Σ::AbstractMatrix{Num})
    Σinv = try
        symbolic_inverse(Σ)
    catch err
        throw(ArgumentError("the transition matrix Σ could not be inverted symbolically " *
                            "($(sprint(showerror, err))); is every infected state eventually left?"))
    end
    return tidy(-T * Σinv)
end

function _large_domain(T::AbstractMatrix{<:Real}, Σ::AbstractMatrix{<:Real})
    return -T / Σ
end

function _restrict(K_L::AbstractMatrix{Num}, E)
    Matrix{Num}(tidy(Num.(E' * K_L * E); fractions = false))
end
_restrict(K_L::AbstractMatrix{<:Real}, E) = E' * K_L * E

const TERM_STRATEGIES = (:auto, :uninfected_dependence, :nonlinear_in_infected)

"""
    next_generation_matrix(sys, infected; equilibrium = nothing, transmission = :auto,
                           warn = true)

Construct the [`NextGenerationMatrix`](@ref) of the ModelingToolkit system `sys`
(a `System` built from ordinary differential equations, or a Catalyst `ReactionSystem`
when Catalyst is loaded).

# Arguments

  - `infected`: the infected state variables, as symbolic variables (`[E, I]`, `[sys.E, sys.I]`)
    or as `Symbol`s (`[:E, :I]`). Their order fixes the row and column order of `T`, `Σ`
    and `K_L`. [`suggest_infected`](@ref) proposes candidates.
  - `equilibrium`: the state at which to linearise, as a `Dict`/pairs mapping unknowns to
    values or expressions (`Dict(S => N)`). Infected states default to zero. When
    omitted, [`disease_free_equilibrium`](@ref) is used. Any other state may be given to
    obtain a reproduction number at that state (for instance an effective reproduction
    number with partial immunity); `Dict()` leaves every uninfected state symbolic.
  - `transmission`: how to decide which terms of the infected equations are new
    infections. The named strategies are `:auto` (default; a term is a transmission if it
    [`depends_on_uninfected`](@ref) states or is [`nonlinear_in_infected`](@ref) states),
    `:uninfected_dependence` and `:nonlinear_in_infected` (each rule alone). A `Function`
    is called on each additive term (a `Num`) and must return `Bool`. A `Vector` of
    expressions gives the new-infection rate `F_i` of each infected compartment explicitly
    (van den Driessche & Watmough's `𝓕`); the transitions are then the remainder of each
    equation. For Catalyst models the default is `:stoichiometry` (a reaction is a
    transmission if it increases the total number of infected individuals,
    [`reaction_is_transmission`](@ref)), a `Function` is called on each `Reaction`, and a
    `Vector{Bool}` or vector of reaction indices marks the transmission reactions; the
    term-based strategies above are also accepted and then applied to the network's ODEs.
  - `autonomous`: when `false`, explicit dependence of the equations on time is allowed (the
    independent variable then appears in `T` and `Σ`); this is what
    [`effective_reproduction_number`](@ref) uses for models with time-varying rates, and
    is not meaningful for `R₀` at an infection-free steady state.
  - `warn`: warn when an entry of `T` is manifestly negative, which means a term that
    *removes* individuals from an infected compartment was classified as a transmission
    (for example density-dependent death `-d (S + I) I`, which involves an uninfected
    state); such terms must be reassigned with an explicit `transmission`.

The strategy used and the resulting `F` and `G` are stored in the result (see
[`transmission_method`](@ref)) and shown when it is printed. Different choices of what
counts as a transmission give different next-generation matrices and different values of
`R₀`, but they all share the threshold property `R₀ > 1` if and only if the
infection-free steady state is unstable (Diekmann et al. 2010, section 3.1 and appendix A).
"""
function next_generation_matrix(sys::AbstractSystem, infected;
        equilibrium = nothing, transmission = :auto, warn::Bool = true,
        autonomous::Bool = true)
    sub = infected_subsystem(sys, infected; autonomous)
    x, y = sub.infected, sub.uninfected
    eq = _equilibrium(sys, infected, equilibrium, x, y)
    (F, G), method = _split(sub.f_infected, x, y, transmission)
    T, Σ = linearise(F, G, x, eq)
    warn && _warn_negative_transmissions(T, x, F)
    return assemble(x, y, eq, T, Σ; F, G, method)
end

function _warn_negative_transmissions(T::AbstractMatrix{Num}, x, F = Num[])
    # term level: a manifestly negative term classified as a transmission
    for (i, f) in enumerate(F), term in additive_terms(f)
        if manifestly_negative(term)
            @warn "the term $(term) removes individuals from $(x[i]) but was classified as a " *
                  "transmission, so T is not non-negative. Reassign it with the `transmission` " *
                  "keyword (see the documentation on choosing transmissions)."
        end
    end
    # entry level: catches negative linearised entries not visible term by term
    for j in axes(T, 2), i in axes(T, 1)
        if manifestly_negative(T[i, j])
            @warn "T[$(i), $(j)] = $(T[i, j]) is negative: a term that removes individuals from " *
                  "$(x[i]) was classified as a transmission. Reassign it with the `transmission` " *
                  "keyword (see the documentation on choosing transmissions)."
        end
    end
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
    if transmission === nothing || transmission === :auto
        return split_terms(f, term -> default_is_transmission(term, x, y)), :auto
    elseif transmission === :uninfected_dependence
        return split_terms(f, term -> depends_on_uninfected(term, y)), transmission
    elseif transmission === :nonlinear_in_infected
        return split_terms(f, term -> nonlinear_in_infected(term, x)), transmission
    elseif transmission isa Symbol
        throw(ArgumentError("unknown transmission strategy `$(repr(transmission))`; use one of $(TERM_STRATEGIES), a predicate, or a vector of new-infection rates"))
    elseif transmission isa Function
        return split_terms(f, transmission), :predicate
    elseif transmission isa AbstractVector
        length(transmission) == length(f) ||
            throw(ArgumentError("`transmission` must give one new-infection rate per infected state ($(length(f)))"))
        F = Num.(transmission)
        return (F, f .- F), :explicit
    else
        throw(ArgumentError("`transmission` must be a strategy name, a predicate on terms, or a vector of new-infection rates"))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a [`NextGenerationMatrix`](@ref) directly from a transmission matrix `T` and a
transition matrix `Σ` (or, equivalently, van den Driessche & Watmough's `F` and `V` with
`Σ = -V`). `infected` names the rows/columns, as symbolic variables or `Symbol`s. For
numeric matrices the sign conventions are checked with [`validate_decomposition`](@ref)
and a warning is issued if they fail (`check = false` disables this).
"""
function next_generation_matrix(T::AbstractMatrix, Σ::AbstractMatrix;
        infected = [Symbolics.variable(:x, i) for i in 1:size(T, 1)], check::Bool = true)
    x = Num[v isa Symbol ? Symbolics.variable(v) : Num(v) for v in infected]
    symbolic = !(_isnumeric(T) && _isnumeric(Σ))
    ngm = assemble(x, Num[], Dict{Num, Any}(v => 0 for v in x), _promote(T, symbolic),
        _promote(Σ, symbolic); method = :matrices)
    if check && !symbolic && !validate_decomposition(ngm)
        @warn "T and Σ violate the sign conventions of a next-generation matrix decomposition " *
              "(T ≥ 0, Σ with non-negative off-diagonal and non-positive diagonal entries, -Σ⁻¹ ≥ 0); " *
              "the spectral radius of K need not be R₀"
    end
    return ngm
end
# `Num <: Real`, so a symbolic matrix must be recognised before the numeric fallback.
function _isnumeric(A::AbstractMatrix)
    eltype(A) <: Real && !(eltype(A) <: Num) && !any(x -> x isa Num, A)
end
_promote(A, symbolic::Bool) = symbolic ? Matrix{Num}(Num.(A)) : Matrix{Float64}(A)

"""
    next_generation_matrix(F, V, x₀, p; infected = ..., check = true)

Numeric construction for models that are not written symbolically, following van den
Driessche & Watmough (2002): `F(x, p)` returns the vector of new-infection rates into the
infected compartments and `V(x, p)` the net outflow (`𝒱⁻ - 𝒱⁺`), so that the infected
subsystem is `ẋ = F(x, p) - V(x, p)`. Both are differentiated with ForwardDiff at the
infection-free state `x₀` (normally zeros), giving `T = ∂F/∂x` and `Σ = -∂V/∂x`. The
uninfected state enters through `p` or through closures.
"""
function next_generation_matrix(F::Function, V::Function, x₀::AbstractVector, p;
        infected = [Symbolics.variable(:x, i) for i in eachindex(x₀)], check::Bool = true)
    T = ForwardDiff.jacobian(x -> F(x, p), x₀)
    Σ = -ForwardDiff.jacobian(x -> V(x, p), x₀)
    ngm = next_generation_matrix(T, Σ; infected, check = false)
    ngm = assemble(ngm.infected, Num[], ngm.equilibrium, ngm.T, ngm.Σ; method = :functions)
    check && !validate_decomposition(ngm) &&
        @warn "F and V violate the sign conventions of a next-generation matrix decomposition"
    return ngm
end

"""
$(TYPEDSIGNATURES)

The next-generation matrix with small domain `K_S = -R Σ⁻¹ C` obtained from a rank
factorisation `T = C R` (Diekmann et al. 2010, section 3.3, equation 2.13). `R` holds a
maximal set of linearly independent rows of `T` and `C = T Rᵀ (R Rᵀ)⁻¹`, so `K_S` is
`r × r` with `r = rank T` and has the same non-zero eigenvalues as `K`. When the
states-at-infection are entered in fixed proportions (`r = 1`) it is the scalar `R₀`.
"""
function small_domain_matrix(ngm::NextGenerationMatrix)
    T = ngm.T
    basis = _independent_rows(T)
    (isempty(basis) || length(basis) == size(T, 1)) && return ngm.K
    Rm = T[basis, :]
    C = _factor_columns(T, Rm)
    return _small(Rm, ngm.Σ, C)
end

# Greedy selection of a maximal set of linearly independent rows.
function _independent_rows(T::AbstractMatrix)
    basis = Int[]
    for i in axes(T, 1)
        all(_structural_zero, T[i, :]) && continue
        candidate = vcat(basis, i)
        _full_row_rank(T[candidate, :]) && push!(basis, i)
    end
    return basis
end
function _full_row_rank(A::AbstractMatrix{Num})
    r, n = size(A)
    r > n && return false
    # some r×r minor must be non-zero
    for cols in _combinations(n, r)
        symbolic_iszero(det(A[:, cols])) || return true
    end
    return false
end
_full_row_rank(A::AbstractMatrix{<:Real}) = rank(A) == size(A, 1)
function _combinations(n, r)
    r == 0 && return [Int[]]
    out = Vector{Vector{Int}}()
    for c in _combinations(n, r - 1)
        start = isempty(c) ? 1 : c[end] + 1
        for j in start:n
            push!(out, vcat(c, j))
        end
    end
    return out
end
_factor_columns(T::AbstractMatrix{Num}, Rm) = tidy(T * Rm' * symbolic_inverse(Rm * Rm'))
_factor_columns(T::AbstractMatrix{<:Real}, Rm) = T * Rm' / (Rm * Rm')
_small(Rm::AbstractMatrix{Num}, Σ, C) = tidy(-Rm * symbolic_inverse(Σ) * C)
_small(Rm::AbstractMatrix{<:Real}, Σ, C) = -(Rm / Σ) * C

"""
$(TYPEDSIGNATURES)

The matrix `-Σ⁻¹` of expected sojourn times: entry `(i, j)` is the expected time an
individual now in infected state `j` will spend in state `i` over its remaining infected
life (Diekmann et al. 2010, section 3.1).
"""
function mean_sojourn_times(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}})
    tidy(-symbolic_inverse(ngm.Σ))
end
mean_sojourn_times(ngm::NextGenerationMatrix{<:AbstractMatrix{<:Real}}) = -inv(ngm.Σ)

"""
    reaction_is_transmission(k, infected_indices, netstoich)

Default classification of reaction `k` of a Catalyst `ReactionSystem` (available when
Catalyst is loaded): `true` if the net stoichiometry of the reaction (column `k` of
`netstoich`, rows indexed by species) increases the total number of individuals in the
infected species `infected_indices`. `S + I --> E + I` and `I --> I + J` are
transmissions, `E --> I` and `I --> R` are transitions.
"""
function reaction_is_transmission end
