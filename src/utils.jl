# Symbolic helpers shared across the package.

const SymbolicLike = Union{Num, SymbolicUtils.BasicSymbolic}

"""
$(TYPEDSIGNATURES)

Strip the `Num` wrapper if present.
"""
_unwrap(x) = unwrap(x)

"""
$(TYPEDSIGNATURES)

Strip the `Num` wrapper and, if the result is a symbolic constant, return the underlying
Julia number. Symbolic expressions are returned unchanged.
"""
function _value(x)
    y = unwrap(x)
    y isa Number && return y
    SymbolicUtils.isconst(y) && return SymbolicUtils.unwrap_const(y)
    return y
end

"""
$(TYPEDSIGNATURES)

`true` if `x` and `y` are the same symbolic object, ignoring the `Num` wrapper.
"""
_symequal(x, y) = isequal(_unwrap(x), _unwrap(y))

"""
$(TYPEDSIGNATURES)

Index of `x` in the symbolic vector `xs`, or `nothing`.
"""
_findsym(x, xs) = findfirst(y -> _symequal(x, y), xs)

_negative_coefficient(x) = (v = _value(x); v isa Number && v < 0)

# A sum all of whose terms carry a negative coefficient, e.g. `-γ - μ`.
function _isnegsum(x)
    SymbolicUtils.isadd(x) || return false
    function isneg(a)
        _negative_coefficient(a) ||
            (SymbolicUtils.ismul(a) && _negative_coefficient(arguments(a)[1]))
    end
    return all(isneg, arguments(x))
end
_possum(x) = _unwrap(sum(-Num(a) for a in arguments(x)))

"""
$(TYPEDSIGNATURES)

Factor the sign out of an expression: returns `(s, y)` with `s ∈ {-1, 1}` and
`x = s * y`, where `y` has every all-negative sum replaced by its positive counterpart and
no negative numerical coefficient. Symbolics canonicalises `-(γ + μ)` to `-γ - μ`, so
the sign can only be cancelled against another sign in a product or quotient.
"""
function _pull_sign(x)
    if _isnegsum(x)
        return (-1, _possum(x))
    elseif SymbolicUtils.ismul(x)
        s = 1
        factors = Num[]
        for a in arguments(x)
            if _isnegsum(a)
                s = -s
                push!(factors, Num(_possum(a)))
            elseif _negative_coefficient(a)
                s = -s
                push!(factors, Num(-_value(a)))
            else
                push!(factors, Num(a))
            end
        end
        return (s, _unwrap(prod(factors)))
    elseif _negative_coefficient(x)
        return (-1, -_value(x))
    else
        return (1, x)
    end
end

"""
$(TYPEDSIGNATURES)

Cancel paired minus signs in a product or quotient. Symbolic matrix inversion produces
expressions such as `(-β) / (-γ - μ)`; this turns them into `β / (γ + μ)`.
"""
function _negsum(x)
    if SymbolicUtils.ismul(x)
        s, y = _pull_sign(x)
        return s == 1 ? y : _unwrap(-Num(y))
    elseif SymbolicUtils.isdiv(x)
        num, den = arguments(x)
        sn, fn = _pull_sign(num)
        sd, fd = _pull_sign(den)
        s = sn * sd
        return _unwrap((s == 1 ? Num(fn) : -Num(fn)) / Num(fd))
    end
    return x
end

"""
    tidy(x; fractions = true)

Cosmetic normalisation of a symbolic expression: cancel paired minus signs and, if
`fractions` is `true`, simplify products and quotients as rational functions (cancelling
common factors). Sums are handled term by term so that `a / b + c / d` keeps its shape. The
result is mathematically identical to the input.
"""
function tidy(x::Num; fractions::Bool = true)
    y = _value(x)
    y isa Number && return Num(y)
    if SymbolicUtils.isadd(y)
        return sum(tidy(Num(a); fractions) for a in arguments(y))
    end
    z = Postwalk(_negsum)(y)
    if fractions
        z = Postwalk(_negsum)(_unwrap(Symbolics.simplify_fractions(Num(z))))
    end
    return Num(z)
end
tidy(x::Number; kwargs...) = x
tidy(x::AbstractArray; kwargs...) = tidy.(x; kwargs...)

"""
$(TYPEDSIGNATURES)

Split an expression into its top-level additive terms. Products are *not* expanded, so
`β * (N - I) * I / N` stays a single term; a quotient whose numerator is a sum is split
term by term.
"""
function additive_terms(expr::Num)
    x = _unwrap(expr)
    if SymbolicUtils.isadd(x)
        return Num[Num(a) for a in arguments(x)]
    elseif SymbolicUtils.isdiv(x)
        num, den = arguments(x)
        if SymbolicUtils.isadd(num)
            return Num[Num(a) / Num(den) for a in arguments(num)]
        end
    end
    return Num[expr]
end
additive_terms(expr::Number) = Num[Num(expr)]

"""
$(TYPEDSIGNATURES)

The symbolic variables (unknowns and parameters) appearing in `x`, as a `Vector{Num}`.
"""
function symbolic_variables(x)
    vars = Num[]
    for v in get_variables(x)
        _findsym(v, vars) === nothing && push!(vars, Num(v))
    end
    return vars
end
function symbolic_variables(A::AbstractArray)
    vars = Num[]
    for x in A
        for v in symbolic_variables(x)
            _findsym(v, vars) === nothing && push!(vars, v)
        end
    end
    return vars
end

"""
$(TYPEDSIGNATURES)

`true` if `x` depends on any of the symbolic objects in `vars`.
"""
function depends_on(x, vars)
    isempty(vars) && return false
    return any(v -> _findsym(v, vars) !== nothing, get_variables(x))
end

# `true` if `x` is built only from `+`, `-`, `*`, `/` and integer powers, so that
# rational-function simplification can decide whether it vanishes.
function _isrational(x)
    y = _value(x)
    y isa Number && return true
    iscall(y) || return true
    op = operation(y)
    args = arguments(y)
    if op === (^)
        e = _value(args[2])
        return e isa Integer && _isrational(args[1])
    end
    op in (+, -, *, /) || return false
    return all(_isrational, args)
end

"""
$(TYPEDSIGNATURES)

Decide whether an expression is identically zero.

The expression is first evaluated at `nprobe` random points; a non-zero value proves that
it is not identically zero. If every probe vanishes and the expression is a rational
function of its symbols, it is simplified as a rational function and the numerator is
expanded; a literal zero then proves that it is identically zero. If that is inconclusive
(or the expression involves `sqrt`, `max`, ...) the answer is `true` only when `numeric` is
`true`: a rational function that is not identically zero vanishes on a set of measure
zero, so agreement at random points is overwhelming, but not formal, evidence.
"""
function symbolic_iszero(x; numeric::Bool = true, nprobe::Int = 4, atol::Real = 1e-9)
    y = _value(x)
    y isa Number && return iszero(y)
    vars = symbolic_variables(Num(y))
    for _ in 1:nprobe
        vals = Dict(v => 0.5 + rand() for v in vars)
        r = _fold(substitute(Num(y), vals))
        r isa Number || break
        abs(r) <= atol || return false
    end
    if _isrational(y)
        z = _value(Symbolics.simplify_fractions(Num(y)))
        z isa Number && return iszero(z)
        # A rational function vanishes iff its numerator does; expanding the numerator to
        # polynomial normal form is cheap, unlike a full `simplify`.
        num = SymbolicUtils.isdiv(z) ? arguments(z)[1] : z
        e = _value(Symbolics.expand(Num(num)))
        e isa Number && return iszero(e)
    end
    numeric || return false
    # all probes vanished (or could not be evaluated to a number)
    for _ in 1:nprobe
        vals = Dict(v => 0.5 + rand() for v in vars)
        r = _fold(substitute(Num(y), vals))
        r isa Number || return false
        abs(r) <= atol || return false
    end
    return true
end

"""
$(TYPEDSIGNATURES)

`true` if the two expressions agree as rational functions of their symbols, using
[`symbolic_iszero`](@ref) on their difference.
"""
symbolic_isequal(x, y; kwargs...) = symbolic_iszero(Num(x) - Num(y); kwargs...)

"""
$(TYPEDSIGNATURES)

Evaluate an expression whose leaves are all numbers. `substitute` leaves calls such as
`sqrt(6.0)` or `max(2.0, 3.0)` unevaluated; this folds them. Returns a `Number`, or the
(partially folded) symbolic expression if some symbol remains.
"""
function _fold(x)
    y = _value(x)
    y isa Number && return y
    if iscall(y)
        args = map(_fold, arguments(y))
        all(a -> a isa Number, args) && return operation(y)(args...)
    end
    return y
end

"""
$(TYPEDSIGNATURES)

Convert a symbolic expression in which every symbol has been substituted by a number into
a `Float64`. Throws an informative `ArgumentError` listing any unresolved symbols.
"""
function to_number(x)
    y = _fold(x)
    y isa Number && return Float64(y)
    vars = symbolic_variables(Num(y))
    throw(ArgumentError("expression `$(Num(y))` still contains the symbols $(vars); " *
                        "supply values for them (e.g. as `Dict(p => value)`)"))
end

"""
$(TYPEDSIGNATURES)

Build a `Dict{Num,Any}` substitution map from pairs, a `Dict`, or any object supporting the
SymbolicIndexingInterface (an `ODEProblem`, a solution, an integrator, ...), for the
symbols in `syms`.
"""
function substitution_map(p::AbstractDict, syms = nothing)
    return Dict{Num, Any}(Num(k) => v for (k, v) in p)
end
function substitution_map(p::AbstractVector{<:Pair}, syms = nothing)
    substitution_map(Dict(p), syms)
end
function substitution_map(p::Tuple{Vararg{Pair}}, syms = nothing)
    substitution_map(Dict(p...), syms)
end
function substitution_map(p, syms)
    syms === nothing &&
        throw(ArgumentError("cannot build a substitution map from `$(typeof(p))` without knowing which symbols are needed"))
    d = Dict{Num, Any}()
    for s in syms
        if is_parameter(p, s)
            d[s] = getp(p, s)(p)
        elseif is_variable(p, s)
            d[s] = getu(p, s)(p)
        end
    end
    return d
end
