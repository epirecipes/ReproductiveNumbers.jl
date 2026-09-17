# Spectral radii of non-negative matrices, symbolically where a closed form exists.

"""
$(TYPEDSIGNATURES)

Partition the indices of the square matrix `A` into the strongly connected components of
the directed graph with an edge `j → i` whenever `A[i, j]` is not structurally zero.
Components are returned in an order in which `A` is block upper triangular, so that the
spectrum of `A` is the union of the spectra of the diagonal blocks `A[c, c]`.
"""
function irreducible_blocks(A::AbstractMatrix)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("matrix must be square"))
    adj = [Int[] for _ in 1:n]
    for j in 1:n, i in 1:n
        i != j && !_structural_zero(A[i, j]) && push!(adj[j], i)
    end
    # Tarjan's algorithm; components are emitted in reverse topological order.
    index = zeros(Int, n)
    low = zeros(Int, n)
    onstack = falses(n)
    stack = Int[]
    comps = Vector{Vector{Int}}()
    counter = Ref(0)
    function strongconnect(v)
        counter[] += 1
        index[v] = low[v] = counter[]
        push!(stack, v)
        onstack[v] = true
        for w in adj[v]
            if index[w] == 0
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elseif onstack[w]
                low[v] = min(low[v], index[w])
            end
        end
        if low[v] == index[v]
            comp = Int[]
            while true
                w = pop!(stack)
                onstack[w] = false
                push!(comp, w)
                w == v && break
            end
            push!(comps, sort!(comp))
        end
    end
    for v in 1:n
        index[v] == 0 && strongconnect(v)
    end
    return comps
end

"""
$(TYPEDSIGNATURES)

`true` if `A` has rank at most one, tested through its `2 × 2` minors. For symbolic
matrices the minors are simplified as rational functions and, if `numeric` is `true`, also
checked at random parameter values (see [`symbolic_iszero`](@ref)).
"""
function isrankone(A::AbstractMatrix{Num}; numeric::Bool = true)
    n = size(A, 1)
    for i in 1:n, j in (i + 1):n, k in 1:n, l in (k + 1):n
        m = A[i, k] * A[j, l] - A[i, l] * A[j, k]
        symbolic_iszero(m; numeric) || return false
    end
    return true
end
isrankone(A::AbstractMatrix{<:Real}; kwargs...) = rank(A) <= 1

"""
$(TYPEDSIGNATURES)

The characteristic polynomial `det(λI - A)` of a symbolic square matrix, as an expression
in the symbolic variable `λ` (`Symbolics.variable(:λ)` by default).
"""
function characteristic_polynomial(A::AbstractMatrix{Num}; λ::Num = Symbolics.variable(:λ))
    n = size(A, 1)
    return tidy(Symbolics.expand(det(λ * Matrix{Num}(I, n, n) - A)); fractions = false)
end

"""
$(TYPEDSIGNATURES)

Dominant eigenvalue of an irreducible non-negative `2 × 2` matrix `[a b; c d]`,
`(a + d + √((a - d)² + 4 b c)) / 2` (equation 2.12 of Diekmann et al. 2010). The
discriminant is written as `(a - d)² + 4 b c` to make its non-negativity manifest.
"""
function spectral_radius_2x2(a, b, c, d)
    if _structural_zero(a) && _structural_zero(d)
        return sqrt(b * c)
    end
    return (a + d + sqrt((a - d)^2 + 4 * b * c)) / 2
end

function _irreducible_spectral_radius(B::AbstractMatrix{Num}; numeric::Bool)
    n = size(B, 1)
    n == 1 && return B[1, 1]
    isrankone(B; numeric) && return tr(B)
    n == 2 && return spectral_radius_2x2(B[1, 1], B[1, 2], B[2, 1], B[2, 2])
    throw(NoClosedFormError("no closed form is available for the dominant eigenvalue of the " *
                            "irreducible $(n)×$(n) block\n$(sprint(show, MIME("text/plain"), B))\n" *
                            "which is neither rank one nor 2×2; evaluate numerically with " *
                            "`basic_reproduction_number(ngm, parameter_values)`"))
end

"""
    spectral_radius(A::AbstractMatrix{Num}; numeric_rank_check = true)
    spectral_radius(A::AbstractMatrix{<:Real})

Spectral radius of a non-negative square matrix.

For a numeric matrix this is `maximum(abs, eigvals(A))`.

For a symbolic matrix a closed form is built when possible: `A` is permuted to block
triangular form by [`irreducible_blocks`](@ref), and for each irreducible diagonal block a
closed form is used when the block is `1 × 1`, has rank one (spectral radius equal to its
trace) or is `2 × 2` ([`spectral_radius_2x2`](@ref)). The spectral radius of `A` is then the
maximum over blocks, returned as a symbolic `max` when it cannot be resolved. If some
block admits no closed form a [`NoClosedFormError`](@ref) is thrown.

The rank-one test uses symbolic simplification and, if `numeric_rank_check` is `true`,
evaluation at random parameter values as a fallback.
"""
function spectral_radius(A::AbstractMatrix{Num}; numeric_rank_check::Bool = true)
    n = size(A, 1)
    n == 0 && return Num(0)
    vals = Num[]
    for comp in irreducible_blocks(A)
        ρ = tidy(_irreducible_spectral_radius(A[comp, comp]; numeric = numeric_rank_check);
            fractions = false)
        any(v -> symbolic_isequal(v, ρ; numeric = numeric_rank_check), vals) ||
            push!(vals, ρ)
    end
    # drop zeros when there are other candidates: ρ ≥ 0 for non-negative matrices
    nonzero = filter(v -> !_structural_zero(v), vals)
    isempty(nonzero) && return Num(0)
    length(nonzero) == 1 && return nonzero[1]
    return tidy(reduce(max, nonzero); fractions = false)
end
function spectral_radius(A::AbstractMatrix{<:Real})
    size(A, 1) == 0 && return 0.0
    return maximum(abs, eigvals(Matrix{Float64}(A)))
end
