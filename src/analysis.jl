# Further analysis: suggested infected sets, sensitivities, Perron vectors and the
# effective reproduction number along solutions.

"""
    suggest_infected(sys; all = false, maxstates = 14)

Propose the infected compartments of `sys`. A set `X` of states is a candidate if

  - setting every state in `X` to zero makes the equations of `X` vanish identically (the
    infection-free subspace is invariant),
  - some equation of `X` has a term that brings individuals into `X` (a term that depends
    on a state outside `X` or is non-linear in `X`, and is not manifestly negative),
  - no state of `X` leaves `X` only through contact with states outside `X` (which would
    make `X` a set of *susceptible* classes emptied by infection rather than a set of
    infected classes emptied by recovery and death), and
  - `X` is minimal with these properties.

Subsets are enumerated by increasing size, so the system may have at most `maxstates`
unknowns. By default the union of the minimal candidates is returned as one vector (so
that, for instance, the infected states of every strain of a multi-strain model are
included); with `all = true` the list of minimal candidates is returned instead. The
answer is a heuristic and should be checked against the model's meaning.
"""
function suggest_infected(sys::AbstractSystem; all::Bool = false, maxstates::Int = 14)
    states, rhs = ode_right_hand_sides(sys)
    n = length(states)
    n <= maxstates ||
        throw(ArgumentError("suggest_infected enumerates subsets of the $(n) unknowns; raise `maxstates` (currently $(maxstates)) if you really want this"))
    minimal = Vector{Vector{Int}}()
    for size in 1:(n == 1 ? 1 : n - 1), X in _subsets(n, size)
        Base.any(Y -> issubset(Y, X), minimal) && continue
        _invariant(states, rhs, X) || continue
        _has_inflow(states, rhs, X) || continue
        _outflow_within(states, rhs, X) || continue
        push!(minimal, X)
    end
    all && return [states[X] for X in minimal]
    return states[sort!(unique!(reduce(vcat, minimal; init = Int[])))]
end

_subsets(n, k) = _combinations(n, k)

# Every equation of X vanishes identically when the states of X are zero.
function _invariant(states, rhs, X)
    zero_X = Dict{Num, Any}(states[i] => 0 for i in X)
    return Base.all(i -> symbolic_iszero(substitute(rhs[i], zero_X); numeric = false), X)
end

# Some equation of X has a term that brings individuals into X from outside.
function _has_inflow(states, rhs, X)
    n = length(states)
    outside = states[setdiff(1:n, X)]
    inside = states[X]
    return Base.any(X) do i
        Base.any(additive_terms(rhs[i])) do term
            (depends_on(term, outside) || nonlinear_in_infected(term, inside)) &&
                !manifestly_negative(term)
        end
    end
end

# No state of X may have all of its outflow terms depend on states outside X: such a
# state is emptied by contact with outsiders, i.e. it is a susceptible class.
function _outflow_within(states, rhs, X)
    outside = states[setdiff(1:length(states), X)]
    for i in X
        negs = filter(manifestly_negative, additive_terms(rhs[i]))
        isempty(negs) && continue
        Base.all(term -> depends_on(term, outside), negs) && return false
    end
    return true
end

"""
    sensitivities(R0::Num, params = symbolic_variables(R0))
    sensitivities(ngm::NextGenerationMatrix, [p]; params = ...)

Partial derivatives `∂R₀/∂θ` of a closed-form reproduction number with respect to each
parameter `θ` in `params`, as a `Dict`. With a next-generation matrix the closed form is
taken from [`basic_reproduction_number`](@ref); with parameter values `p` the derivatives
are evaluated numerically. See also [`elasticities`](@ref).
"""
function sensitivities(R0::Num, params = symbolic_variables(R0))
    return Dict{Num, Num}(θ => tidy(Symbolics.derivative(R0, θ); fractions = false)
    for θ in params)
end
function sensitivities(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}};
        params = nothing, kwargs...)
    R0 = basic_reproduction_number(ngm; kwargs...)
    return sensitivities(R0, params === nothing ? symbolic_variables(R0) : params)
end
function sensitivities(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}}, p;
        params = nothing, kwargs...)
    s = sensitivities(ngm; params, kwargs...)
    d = substitution_map(p, collect(keys(s)))
    return Dict{Num, Float64}(θ => to_number(substitute(v, d)) for (θ, v) in s)
end

"""
    elasticities(R0::Num, params = symbolic_variables(R0))
    elasticities(ngm::NextGenerationMatrix, [p]; params = ...)

Elasticities `(θ / R₀) ∂R₀/∂θ`, the proportional change in `R₀` per proportional change in
each parameter, as a `Dict`. See [`sensitivities`](@ref).
"""
function elasticities(R0::Num, params = symbolic_variables(R0))
    # no rational-function simplification: R₀ may contain square roots, for which it
    # is extremely slow
    return Dict{Num, Num}(θ => tidy(θ / R0 * Symbolics.derivative(R0, θ); fractions = false)
    for θ in params)
end
function elasticities(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}};
        params = nothing, kwargs...)
    R0 = basic_reproduction_number(ngm; kwargs...)
    return elasticities(R0, params === nothing ? symbolic_variables(R0) : params)
end
function elasticities(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}}, p;
        params = nothing, kwargs...)
    e = elasticities(ngm; params, kwargs...)
    d = substitution_map(p, collect(keys(e)))
    return Dict{Num, Float64}(θ => to_number(substitute(v, d)) for (θ, v) in e)
end

"""
    perron_vectors(ngm::NextGenerationMatrix)
    perron_vectors(K::AbstractMatrix)

The right and left eigenvectors `(u, v)` of the next-generation matrix `K` for its
dominant eigenvalue `R₀`, indexed by the states-at-infection. `u`, normalised to sum to
one, is the distribution of states-at-infection among new infections once the epidemic
grows exponentially; `v`, normalised so that `vᵀu = 1`, holds the reproductive values of
the types. Numeric matrices are handled by `eigen`; symbolic ones in closed form when `K`
is `1 × 1`, `2 × 2` or rank one, otherwise a [`NoClosedFormError`](@ref) is thrown.
"""
perron_vectors(ngm::NextGenerationMatrix) = perron_vectors(ngm.K)
function perron_vectors(K::AbstractMatrix{<:Real})
    n = size(K, 1)
    n == 0 && return (Float64[], Float64[])
    ev = eigen(Matrix{Float64}(K))
    i = argmax(real.(ev.values))
    u = real.(ev.vectors[:, i])
    u ./= sum(u)
    evl = eigen(Matrix{Float64}(K'))
    j = argmax(real.(evl.values))
    v = real.(evl.vectors[:, j])
    v ./= dot(v, u)
    return (u, v)
end
function perron_vectors(K::AbstractMatrix{Num})
    n = size(K, 1)
    if n == 1
        return (Num[1], Num[1])
    elseif isrankone(K)
        # K = u vᵀ: pick a non-zero column as u and the matching row as v
        j = findfirst(j -> !all(_structural_zero, K[:, j]), 1:n)
        i = findfirst(i -> !_structural_zero(K[i, j]), 1:n)
        u = tidy(K[:, j] ./ sum(K[:, j]); fractions = false)
        v = tidy(K[i, :] ./ K[i, j]; fractions = false)
        return (u, tidy(v ./ dot(v, u)))
    elseif n == 2
        a, b, c, d = K[1, 1], K[1, 2], K[2, 1], K[2, 2]
        ρ = spectral_radius_2x2(a, b, c, d)
        u = Num[b, ρ - a]
        v = Num[c, ρ - a]
        if _structural_zero(b) && _structural_zero(c)
            throw(NoClosedFormError("K is diagonal; the Perron vector is the unit vector of the larger diagonal entry, which cannot be chosen symbolically"))
        end
        u = tidy(u ./ sum(u); fractions = false)
        v = tidy(v ./ dot(v, u); fractions = false)
        return (u, v)
    end
    throw(NoClosedFormError("no closed-form Perron vectors for a $(n)×$(n) matrix that is not rank one; evaluate numerically"))
end

"""
    effective_reproduction_number(sys, infected; equilibrium = Dict(), kwargs...)
    effective_reproduction_number(ngm)
    effective_reproduction_number(ngm, state)
    effective_reproduction_number(ngm, sol, t)
    effective_reproduction_number(ngm, sol, ts::AbstractVector)
    effective_reproduction_number(ngm, sol)

The effective reproduction number `R_t`: the expected number of new infections caused by an
infectious individual in a population in which some individuals are no longer susceptible.
It is computed exactly like `R₀`, as the spectral radius of the next-generation matrix, but
linearised at the *current* state instead of the infection-free steady state, so it is a
function of the uninfected states (for the SIR model `R_t = R₀ S(t)/N`).

  - `effective_reproduction_number(sys, infected)` builds the next-generation matrix with
    the uninfected states left symbolic (`equilibrium = Dict()`; give some of them to fix
    them) and returns `R_t` as a closed-form expression in those states when one exists
    (otherwise a [`NoClosedFormError`](@ref) is thrown). Rates that depend explicitly on
    time are allowed (`autonomous = false` by default here), so a seasonally forced `β(t)`
    gives `R_t = β(t) S(t)/(γ N)`. Other keyword arguments are those of
    [`next_generation_matrix`](@ref).
  - With a symbolic `ngm` built that way, `effective_reproduction_number(ngm)` is the same
    expression, and `effective_reproduction_number(ngm, state)` evaluates it at a state
    given as a `Dict` or pairs mapping (some of) the states and parameters to values or
    expressions; the result is a `Float64` when everything is numeric and an expression
    otherwise.
  - With an ODE solution `sol` (of the same system), `effective_reproduction_number(ngm, sol, t)`
    evaluates it at time `t`, reading the parameters from the solution and the state from
    `sol(t)`; a vector of times gives the trajectory, and `sol` alone uses `sol.t`. No closed
    form is needed here: the numeric spectral radius is used.

To obtain `R_t` directly from the solver, as a variable that can be indexed and plotted,
see [`add_effective_reproduction_number`](@ref).

!!! note "Which R_t?"

    This is the *instantaneous* reproduction number of the compartmental model, the
    threshold quantity for growth of the linearised infected subsystem at the current state
    (what `R₀ S(t)/N` denotes in textbooks). It is not the *case* reproduction number
    estimated from incidence data through a renewal equation, which averages over the
    infectious period of the cases infected at time `t`.
"""
function effective_reproduction_number(sys::AbstractSystem, infected;
        equilibrium = Dict{Num, Any}(), autonomous::Bool = false, kwargs...)
    ngm = next_generation_matrix(sys, infected; equilibrium, autonomous, kwargs...)
    return basic_reproduction_number(ngm)
end
function effective_reproduction_number(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}};
        kwargs...)
    return basic_reproduction_number(ngm; kwargs...)
end
function effective_reproduction_number(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}},
        state::Union{AbstractDict, AbstractVector{<:Pair}, Tuple{Vararg{Pair}}})
    d = substitution_map(state)
    syms = symbolic_variables([ngm.T; ngm.Σ])
    if Base.all(s -> haskey(d, s), syms)
        return basic_reproduction_number(ngm, d)
    end
    R = basic_reproduction_number(ngm)
    return tidy(substitute(R, d); fractions = false)
end
function effective_reproduction_number(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}},
        sol::AbstractTimeseriesSolution, t::Real)
    syms = symbolic_variables([ngm.T; ngm.Σ])
    d = Dict{Num, Any}()
    for s in syms
        if is_parameter(sol, s)
            d[s] = getp(sol, s)(sol)
        elseif is_variable(sol, s)
            d[s] = sol(t; idxs = s)
        elseif is_independent_variable(sol, s)
            d[s] = t   # time-varying rates
        end
    end
    return basic_reproduction_number(ngm, d)
end
function effective_reproduction_number(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}},
        sol::AbstractTimeseriesSolution, ts::AbstractVector)
    return [effective_reproduction_number(ngm, sol, t) for t in ts]
end
function effective_reproduction_number(ngm::NextGenerationMatrix{<:AbstractMatrix{Num}},
        sol::AbstractTimeseriesSolution)
    return effective_reproduction_number(ngm, sol, sol.t)
end

"""
    add_effective_reproduction_number(sys, infected; name = :Rt, kwargs...)

Return `(sys′, Rt)`: a compiled copy of the ModelingToolkit system `sys` with an observed
variable `Rt(t)` equal to the closed-form effective reproduction number
([`effective_reproduction_number`](@ref)) as a function of the state, and the variable
itself. Solutions of `sys′` then give the trajectory directly, `sol[Rt]`, and it can be
plotted with `idxs = Rt`. The system is rebuilt from `equations(sys)` with its initial
conditions and compiled with `mtkcompile`, so an uncompiled system may be given; rates
that depend explicitly on time are allowed. Keyword arguments are
those of [`next_generation_matrix`](@ref); the closed form must exist.
"""
function add_effective_reproduction_number(sys::AbstractSystem, infected;
        name::Symbol = :Rt, kwargs...)
    R = effective_reproduction_number(sys, infected; kwargs...)
    iv = get_iv(sys)
    Rt = only(@variables $(name)(iv))
    eqs = [collect(equations(sys)); Rt ~ R]
    newsys = System(eqs, iv; name = nameof(sys),
        initial_conditions = ModelingToolkit.initial_conditions(sys))
    return mtkcompile(newsys), Rt
end
