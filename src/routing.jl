# Routing of the symbolic algorithms by problem size, with retries.
#
# Every symbolic step of the pipeline has more than one way of being done, with different
# costs: Laplace expansion gives compact inverses but is factorial in the size, LU is cubic;
# rational-function simplification is exact but can take minutes on large expressions;
# the closed form of a spectral radius may come from the blocks of K or from the smaller
# K_S. The helpers here choose an order of methods from the size of the problem and fall
# through to the next method when one fails. Set `ENV["JULIA_DEBUG"] = "ReproductiveNumbers"`
# to log the route taken.

"""
$(TYPEDSIGNATURES)

Try `methods`, a vector of `name => thunk` pairs, in order and return the first result;
each failure is logged at debug level and the next method is tried. If all fail, the last
exception is rethrown unless `onfail` is given, in which case `onfail(errors)` is called.
"""
function first_success(methods; what::AbstractString = "computation", onfail = nothing)
    errors = Pair{Symbol, Exception}[]
    for (name, thunk) in methods
        try
            result = thunk()
            @debug "$(what): method $(name) succeeded"
            return result
        catch err
            @debug "$(what): method $(name) failed" exception=err
            push!(errors, name => err)
        end
    end
    onfail === nothing && throw(last(errors).second)
    return onfail(errors)
end

"""
$(TYPEDSIGNATURES)

Number of nodes of a symbolic expression tree (a rough measure of its size).
"""
function expression_size(x)
    y = _value(x)
    y isa Number && return 1
    iscall(y) || return 1
    return 1 + sum(expression_size, arguments(y); init = 0)
end
expression_size(A::AbstractArray) = sum(expression_size, A; init = 0)

"""
Largest size (in expression nodes, see [`expression_size`](@ref)) of an entry for which
[`tidy`](@ref) still runs `Symbolics.simplify_fractions`, which is exact but can take minutes
on large rational functions. Larger entries only get their signs normalised.
"""
const TIDY_SIZE_LIMIT = Ref(2_000)

"""
Largest matrix size for which the symbolic inverse is attempted by Laplace expansion
(`inv(::Matrix{Num})`, compact results, factorial cost) before LU factorisation.
"""
const LAPLACE_INVERSE_LIMIT = Ref(4)

"""
    symbolic_inverse(A; method = :auto)

Inverse of a symbolic square matrix. `method` is `:laplace` (the adjugate formula used by
`inv(::Matrix{Num})`: compact expressions, factorial cost), `:lu` (`Symbolics.sym_lu`, cubic
cost) or `:auto`, which uses Laplace expansion for matrices up to
`LAPLACE_INVERSE_LIMIT[]` (4 by default) and LU above that, in each case retrying with
the other method if the first fails (a structurally singular pivot, say). Matrices larger
than eight are never expanded by Laplace.
"""
function symbolic_inverse(A::AbstractMatrix{Num}; method::Symbol = :auto)
    n = size(A, 1)
    n == 1 && return Matrix{Num}(reshape([1 / A[1, 1]], 1, 1))
    laplace = () -> Matrix{Num}(inv(A))
    lu_ = () -> Matrix{Num}(lu(A) \ Matrix{Num}(I, n, n))
    if method === :laplace
        return laplace()
    elseif method === :lu
        return lu_()
    elseif method === :auto
        ladder = if n <= LAPLACE_INVERSE_LIMIT[]
            [:laplace => laplace, :lu => lu_]
        elseif n <= 8
            [:lu => lu_, :laplace => laplace]
        else
            [:lu => lu_]
        end
        return first_success(ladder; what = "symbolic inverse of a $(n)×$(n) matrix")
    end
    throw(ArgumentError("unknown inverse method `$(repr(method))`; use :auto, :laplace or :lu"))
end
