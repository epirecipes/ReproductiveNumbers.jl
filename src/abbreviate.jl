# Readable forms of next-generation matrices through named abbreviations.

"""
    abbreviate(ngm; transmissions = false)

Rewrite a symbolic [`NextGenerationMatrix`](@ref) in terms of named epidemiological
quantities read off the transition matrix `Σ`, following the "epidemiological reasoning"
of Diekmann et al. (2010):

  - `τ_x = -1 / Σ[x, x]`, the mean sojourn time in infected state `x`;
  - `p_x_y = Σ[y, x] τ_x`, the probability that an individual leaving `x` moves to `y`
    (read "from `x` to `y`").

Since `Σ = (P - I) D⁻¹` with `D = diag(τ)`, every entry of `-Σ⁻¹`, `K_L` and `K` becomes a
polynomial in the `τ` and `p` (rational only when transitions form cycles), for instance
`K = [β τ_I p_E_I]` for the SEIR model. With `transmissions = true`, non-zero entries of
`T` that are not single symbols are also named `T_y_x` for `T[y, x]` (new infections in
`y` caused by `x`). Names are plain identifiers so that they print without quoting.

The result has the same fields as `ngm` with the matrices rewritten, and `definitions`
holding each new symbol and its expression in the original parameters. [`expand_definitions`](@ref)
reverses the abbreviation, [`evaluate`](@ref) accepts values for either the original
parameters or the abbreviations, and printing lists the definitions.
"""
function abbreviate(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}};
        transmissions::Bool = false)
    x = ngm.infected
    n = length(x)
    Σ = ngm.Σ
    defs = Dict{Num, Num}(ngm.definitions)
    newΣ = Matrix{Num}(zeros(Num, n, n))
    τ = Vector{Num}(undef, n)        # the symbols
    τexpr = [-1 / Σ[i, i] for i in 1:n]   # their definitions in the original parameters
    for i in 1:n
        τ[i] = _named(Symbol("τ_", _basename(x[i])), τexpr[i], defs)
        newΣ[i, i] = -1 / τ[i]
    end
    for i in 1:n, j in 1:n
        i == j && continue
        _structural_zero(Σ[j, i]) && continue
        pij = _named(
            Symbol("p_", _basename(x[i]), "_", _basename(x[j])), Σ[j, i] * τexpr[i], defs)
        newΣ[j, i] = pij / τ[i]
    end
    newT = Matrix{Num}(copy(ngm.T))
    if transmissions
        for i in 1:n, j in 1:n
            entry = ngm.T[i, j]
            (_structural_zero(entry) || Symbolics.issym(_unwrap(entry))) && continue
            newT[i, j] = _named(
                Symbol("T_", _basename(x[i]), "_", _basename(x[j])), entry, defs)
        end
    end
    return assemble(x, ngm.uninfected, ngm.equilibrium, newT, newΣ;
        F = ngm.F, G = ngm.G, method = ngm.method, definitions = defs)
end

# A fresh symbol `name` standing for `expr` (a single symbol is left as it is); the
# definition is recorded in `defs`. Names already used for a different expression get a
# numeric suffix.
function _named(name::Symbol, expr, defs::Dict{Num, Num})
    expr = tidy(Num(expr))
    Symbolics.issym(_unwrap(expr)) && return expr
    for (k, v) in defs
        symbolic_isequal(v, expr) && return k
    end
    base = name
    k = 1
    sym = Symbolics.variable(name)
    while haskey(defs, sym)
        k += 1
        sym = Symbolics.variable(Symbol(base, "_", k))
    end
    defs[sym] = expr
    return sym
end

"""
    abbreviate(ngm, definitions; solve_for = nothing)

Rewrite a symbolic [`NextGenerationMatrix`](@ref) in terms of user-chosen quantities.
`definitions` is a vector of pairs `symbol => expression` in the original parameters, for
example `[R_E => β / (γ + μ), p_E => σ / (σ + μ)]`. Each definition is solved for one of the
original parameters it contains (`β = R_E (γ + μ)`, or `σ = p_E μ / (1 - p_E)`), which is
then eliminated from the matrices; `solve_for` may give the parameter to eliminate for
each definition (a vector with one entry per definition, `nothing` entries meaning
"choose"), otherwise the first parameter for which the definition can be solved with
`Symbolics.symbolic_linear_solve` is used. Definitions that cannot be solved this way
(quadratic in every parameter, say) are rejected. The new symbols must not be existing parameters, and each
definition must eliminate a distinct parameter.

The result carries the definitions, see the one-argument [`abbreviate`](@ref); the two
forms can be combined by calling one after the other.
"""
function abbreviate(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}},
        definitions::Union{AbstractVector{<:Pair}, AbstractDict}; solve_for = nothing)
    defs = Dict{Num, Num}(ngm.definitions)
    params = symbolic_variables([ngm.T; ngm.Σ])
    pairs = collect(definitions)
    solve_for === nothing && (solve_for = fill(nothing, length(pairs)))
    length(solve_for) == length(pairs) ||
        throw(ArgumentError("`solve_for` must have one entry per definition"))
    subs = Dict{Num, Num}()
    eliminated = Num[]
    for ((sym, expr), target) in zip(pairs, solve_for)
        sym = Num(sym)
        expr = tidy(Num(expr))
        Symbolics.issym(_unwrap(sym)) ||
            throw(ArgumentError("`$(sym)` is not a symbol"))
        _findsym(sym, params) === nothing ||
            throw(ArgumentError("`$(sym)` is already a parameter of the model; choose a new name"))
        candidates = target === nothing ? symbolic_variables(expr) : [Num(target)]
        solution = nothing
        for θ in candidates
            depends_on(expr, [θ]) || continue
            _findsym(θ, eliminated) === nothing || continue
            sol = _solve_definition(sym, expr, θ)
            sol === nothing && continue
            # verify: substituting the solution back must give the new symbol
            symbolic_isequal(substitute(expr, Dict(θ => sol)), sym) || continue
            solution = (θ, sol)
            break
        end
        solution === nothing &&
            throw(ArgumentError("the definition `$(sym) = $(expr)` cannot be solved for a parameter to eliminate" *
                                (target === nothing ? "" : " (`$(target)`)")))
        θ, val = solution
        push!(eliminated, θ)
        # earlier eliminations may appear in this solution
        subs[θ] = tidy(substitute(val, subs))
        for (k, v) in subs
            subs[k] = tidy(substitute(v, Dict(θ => subs[θ])))
        end
        defs[sym] = expr
    end
    newT = tidy(substitute.(ngm.T, Ref(subs)))
    newΣ = tidy(substitute.(ngm.Σ, Ref(subs)))
    return assemble(ngm.infected, ngm.uninfected, ngm.equilibrium, newT, newΣ;
        F = ngm.F, G = ngm.G, method = ngm.method, definitions = defs)
end

# Solve `sym = expr` for the parameter `θ`: directly when `expr` is linear in `θ`, and after
# clearing the denominator when `expr` is a quotient (`σ / (σ + μ)` gives
# `σ = sym μ / (1 - sym)`). Returns `nothing` when neither works.
function _solve_definition(sym, expr, θ)
    attempts = Any[sym - expr]
    u = _unwrap(expr)
    if SymbolicUtils.isdiv(u)
        num, den = arguments(u)
        push!(attempts, sym * Num(den) - Num(num))
    end
    for lhs in attempts
        sol = try
            Symbolics.symbolic_linear_solve(lhs ~ 0, θ)
        catch
            nothing
        end
        sol === nothing && continue
        sol = tidy(Num(sol))
        depends_on(sol, [θ]) && continue
        return sol
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Substitute the [`abbreviate`](@ref)d definitions back into the matrices, returning a
[`NextGenerationMatrix`](@ref) in the original parameters with no definitions.
"""
function expand_definitions(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}})
    isempty(ngm.definitions) && return ngm
    d = Dict{Num, Any}(ngm.definitions)
    # definitions may refer to each other
    for _ in 1:length(d)
        for (k, v) in d
            d[k] = substitute(Num(v), d)
        end
    end
    T = tidy(substitute.(ngm.T, Ref(d)))
    Σ = tidy(substitute.(ngm.Σ, Ref(d)))
    return assemble(ngm.infected, ngm.uninfected, ngm.equilibrium, T, Σ;
        F = ngm.F, G = ngm.G, method = ngm.method)
end
expand_definitions(ngm::NextGenerationMatrix{<:AbstractMatrix{<:Real}}) = ngm

"""
$(TYPEDSIGNATURES)

Substitute the definitions of an abbreviated [`NextGenerationMatrix`](@ref) into an
expression (such as a reproduction number computed from it), giving it in the original
parameters.
"""
function expand_definitions(ngm::NextGenerationMatrix, x::Num)
    isempty(ngm.definitions) && return x
    d = Dict{Num, Any}(ngm.definitions)
    for _ in 1:length(d)
        for (k, v) in d
            d[k] = substitute(Num(v), d)
        end
    end
    return tidy(substitute(x, d); fractions = false)
end
