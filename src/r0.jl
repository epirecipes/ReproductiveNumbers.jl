# Reproduction numbers.

"""
$(TYPEDSIGNATURES)

Substitute numerical values for the symbols of a symbolic [`NextGenerationMatrix`](@ref),
returning a numeric one. `p` may be a `Dict` or vector/tuple of pairs
`parameter => value`, or any object supporting the SymbolicIndexingInterface, such as an
`ODEProblem` built from the same system, from which parameter (and, if the equilibrium
left any state symbolic, initial state) values are read.
"""
function evaluate(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}}, p)
    syms = symbolic_variables([ngm.T; ngm.Σ])
    d = substitution_map(p, syms)
    T = to_number.(substitute.(ngm.T, Ref(d)))
    Σ = to_number.(substitute.(ngm.Σ, Ref(d)))
    eq = Dict{Num, Any}(k => _maybe_number(substitute(Num(v), d))
    for (k, v) in ngm.equilibrium)
    return assemble(ngm.infected, ngm.uninfected, eq, T, Σ; F = ngm.F, G = ngm.G,
        method = ngm.method)
end
evaluate(ngm::NextGenerationMatrix{<:AbstractMatrix{<:Real}}, p) = ngm
function _maybe_number(x)
    y = _fold(x)
    return y isa Number ? Float64(y) : Num(y)
end

"""
    basic_reproduction_number(ngm::NextGenerationMatrix; method = :auto, kwargs...)
    basic_reproduction_number(ngm::NextGenerationMatrix, p)
    basic_reproduction_number(sys, infected, [p]; kwargs...)

The basic reproduction number `R₀`, the spectral radius of the next-generation matrix `K`
(equivalently of `K_L`).

With a symbolic `ngm` and no parameter values, a closed-form expression is returned when
one exists. Two routes are available and `method` selects them: `:blocks` decomposes `K`
into irreducible blocks (see [`spectral_radius`](@ref)); `:small_domain` uses the
small-domain matrix [`small_domain_matrix`](@ref), which is smaller than `K` when `T` has
low rank. `:auto` (the default) tries both, the small domain first when `K` is larger than
`2 × 2` and `T` has low rank, and throws a [`NoClosedFormError`](@ref) if neither
succeeds. Set `ENV["JULIA_DEBUG"] = "ReproductiveNumbers"` to see the route taken. With parameter values `p` (see [`evaluate`](@ref)) a `Float64` is returned.

The three-argument form builds the next-generation matrix first; keyword arguments are
passed to [`next_generation_matrix`](@ref).
"""
function basic_reproduction_number(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}};
        numeric_rank_check::Bool = true, method::Symbol = :auto)
    blocks = () -> spectral_radius(ngm.K; numeric_rank_check)
    small = () -> begin
        K_S = small_domain_matrix(ngm)
        size(K_S, 1) < size(ngm.K, 1) ||
            throw(NoClosedFormError("T has full rank, so the small-domain matrix is K itself"))
        spectral_radius(K_S; numeric_rank_check)
    end
    if method === :blocks
        return blocks()
    elseif method === :small_domain
        return small()
    elseif method === :auto
        # try the smaller matrix first when K is large and T has low rank
        ladder = size(ngm.K, 1) > 2 && _independent_rows(ngm.T) |> length < size(ngm.K, 1) ?
                 [:small_domain => small, :blocks => blocks] :
                 [:blocks => blocks, :small_domain => small]
        return first_success(ladder; what = "closed-form R₀")
    end
    throw(ArgumentError("unknown method `$(repr(method))`; use :auto, :blocks or :small_domain"))
end
function basic_reproduction_number(ngm::NextGenerationMatrix{<:AbstractMatrix{<:Real}})
    spectral_radius(ngm.K)
end
function basic_reproduction_number(ngm::NextGenerationMatrix, p)
    basic_reproduction_number(evaluate(ngm, p))
end
function basic_reproduction_number(sys::AbstractSystem, infected; kwargs...)
    return basic_reproduction_number(next_generation_matrix(sys, infected; kwargs...))
end
function basic_reproduction_number(sys::AbstractSystem, infected, p; kwargs...)
    return basic_reproduction_number(next_generation_matrix(sys, infected; kwargs...), p)
end

"""
    type_reproduction_number(ngm::NextGenerationMatrix, types; kwargs...)
    type_reproduction_number(ngm::NextGenerationMatrix, types, p)

The type reproduction number `T_S` of Roberts & Heesterbeek (2003) for the set `S` of
states-at-infection `types` (symbolic variables, `Symbol`s or indices into
`states_at_infection(ngm)`): the expected number of new infections of the types in `S`
produced by one individual of a type in `S` over the whole course of an infection chain,
counting infections that pass through types outside `S`. With `P` the projection onto `S`,

    T_S = ρ( P K (I - (I - P) K)⁻¹ )

restricted to `S`; for a single type it is a scalar. Control that reduces transmission from
the types in `S` by a factor `1 - 1/T_S` eliminates the infection, and `T_S > 1` if and
only if `R₀ > 1`, provided the infection cannot persist among the remaining types alone
(`ρ((I - P) K) < 1`).
"""
function type_reproduction_number(ngm::NextGenerationMatrix, types; kwargs...)
    idx = _type_indices(ngm, types)
    M = _type_matrix(ngm.K, idx)
    return length(idx) == 1 ? M[1, 1] : _spectral(M; kwargs...)
end
function type_reproduction_number(
        ngm::NextGenerationMatrix{<:AbstractMatrix{Num}}, types, p)
    return type_reproduction_number(evaluate(ngm, p), types)
end
_spectral(M::AbstractMatrix{Num}; kwargs...) = spectral_radius(M; kwargs...)
_spectral(M::AbstractMatrix{<:Real}; kwargs...) = spectral_radius(M)

function _type_indices(ngm, types)
    sa = ngm.states_at_infection
    idx = Int[]
    for tname in (types isa Union{AbstractVector, Tuple} ? types : (types,))
        if tname isa Integer
            1 <= tname <= length(sa) ||
                throw(ArgumentError("type index $(tname) out of range 1:$(length(sa))"))
            push!(idx, tname)
        else
            i = tname isa Symbol ?
                findfirst(v -> Symbol(Symbolics.getname(_unwrap(v))) == tname, sa) :
                _findsym(tname, sa)
            i === nothing &&
                throw(ArgumentError("`$(tname)` is not a state-at-infection; these are $(sa)"))
            push!(idx, i)
        end
    end
    return unique(idx)
end

function _type_matrix(K::AbstractMatrix{Num}, idx)
    # no symbolic check that ρ((I - P) K) < 1 is possible; see the numeric method
    n = size(K, 1)
    P = zeros(Int, n, n)
    for i in idx
        P[i, i] = 1
    end
    Id = Matrix{Num}(I, n, n)
    M = tidy(Num.(P) * K * symbolic_inverse(Id - Num.(I - P) * K))
    return M[idx, idx]
end
function _type_matrix(K::AbstractMatrix{<:Real}, idx)
    n = size(K, 1)
    P = zeros(n, n)
    for i in idx
        P[i, i] = 1.0
    end
    Q = I - P
    spectral_radius(Q * K) < 1 ||
        throw(ArgumentError("the type reproduction number is undefined: the infection persists among the types outside the set (ρ((I - P)K) ≥ 1)"))
    M = P * K * inv(I - Q * K)
    return M[idx, idx]
end

"""
$(TYPEDSIGNATURES)

Check the sign conventions of a numeric decomposition (van den Driessche & Watmough 2002):
`T ≥ 0` entrywise, off-diagonal entries of `Σ` non-negative, diagonal entries of `Σ`
non-positive, and `-Σ⁻¹ ≥ 0`. A symbolic `ngm` is first evaluated at `p`. Returns
`true` if all checks pass; with `verbose = true` each failure is reported.
"""
function validate_decomposition(ngm::NextGenerationMatrix{<:AbstractMatrix{<:Real}};
        verbose::Bool = false, atol::Real = 1e-12)
    ok = true
    report(msg) = (ok = false; verbose && @warn msg)
    all(x -> x >= -atol, ngm.T) || report("T has negative entries")
    n = size(ngm.Σ, 1)
    for i in 1:n, j in 1:n
        if i == j
            ngm.Σ[i, i] <= atol || report("Σ[$i,$i] is positive")
        else
            ngm.Σ[i, j] >= -atol || report("Σ[$i,$j] is negative")
        end
    end
    all(x -> x >= -atol, -inv(ngm.Σ)) || report("-Σ⁻¹ has negative entries")
    return ok
end
function validate_decomposition(
        ngm::NextGenerationMatrix{<:AbstractMatrix{Num}}, p; kwargs...)
    return validate_decomposition(evaluate(ngm, p); kwargs...)
end
