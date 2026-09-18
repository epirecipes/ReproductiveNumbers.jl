# Extracting and splitting the infected subsystem of a ModelingToolkit system.

# `full_equations` substitutes observed variables into the equations; it is only needed when
# there are observed variables, and avoiding it otherwise sidesteps an expression-cache bug
# in ModelingToolkitBase that surfaces when a Catalyst network is converted repeatedly.
function _equations(sys)
    isempty(ModelingToolkit.observed(sys)) ? equations(sys) : full_equations(sys)
end

"""
$(TYPEDSIGNATURES)

Return the ordinary differential equations of `sys` as a pair `(states, rhs)` with the
right-hand sides expressed purely in terms of the unknowns and parameters. Systems with
algebraic (observed) equations are compiled with `mtkcompile` first. Explicit dependence on
the independent variable is an error unless `autonomous = false`.
"""
function ode_right_hand_sides(sys::AbstractSystem; autonomous::Bool = true)
    eqs = _equations(sys)
    states = Num.(unknowns(sys))
    if !all(eq -> isdifferential(_unwrap(eq.lhs)), eqs)
        csys = try
            mtkcompile(sys)
        catch err
            throw(ArgumentError("`$(nameof(sys))` contains non-differential equations and " *
                                "could not be compiled with `mtkcompile`: $(sprint(showerror, err))"))
        end
        eqs = _equations(csys)
        states = Num.(unknowns(csys))
        all(eq -> isdifferential(_unwrap(eq.lhs)), eqs) ||
            throw(ArgumentError("`$(nameof(sys))` is not a system of ordinary differential equations"))
    end
    # Order the right-hand sides by the state that is differentiated on the left.
    rhs = Vector{Num}(undef, length(states))
    filled = falses(length(states))
    for eq in eqs
        var = Num(arguments(_unwrap(eq.lhs))[1])
        i = _findsym(var, states)
        i === nothing &&
            throw(ArgumentError("equation `$(eq)` differentiates a non-unknown"))
        rhs[i] = Num(eq.rhs)
        filled[i] = true
    end
    all(filled) ||
        throw(ArgumentError("some unknowns of `$(nameof(sys))` have no differential equation"))
    iv = Num(get_iv(sys))
    for r in rhs
        autonomous && depends_on(r, [iv]) &&
            throw(ArgumentError("the right-hand side `$(r)` depends explicitly on the independent variable; " *
                                "the basic reproduction number requires an autonomous system " *
                                "(the effective reproduction number does not)"))
    end
    return states, rhs
end

"""
$(TYPEDSIGNATURES)

Resolve `infected` (symbolic variables or `Symbol`s naming unknowns of `sys`) against the
unknowns `states`, returning the matching subset of `states` in the order given.
"""
function resolve_states(sys::AbstractSystem, states::AbstractVector, infected)
    out = Num[]
    for x in infected
        if x isa Symbol
            # match on the (un-namespaced) name so that both complete and incomplete
            # systems, and systems converted from Catalyst, are handled alike
            matches = findall(s -> _basename(s) == x, states)
            isempty(matches) &&
                throw(ArgumentError("no unknown of `$(nameof(sys))` is named `$(x)`; unknowns are $(states)"))
            length(matches) > 1 &&
                throw(ArgumentError("the name `$(x)` is ambiguous in `$(nameof(sys))`: it matches $(states[matches]); pass the symbolic variable instead"))
            i = matches[1]
        else
            i = _findsym(Num(x), states)
            i === nothing &&
                throw(ArgumentError("`$(x)` is not an unknown of `$(nameof(sys))`; unknowns are $(states)"))
        end
        push!(out, states[i])
    end
    allunique(string.(out)) ||
        throw(ArgumentError("infected states must be distinct, got $(out)"))
    return out
end

# Name of a state variable without namespacing: `sys₊E(t)` and `E(t)` both give `:E`.
function _basename(v)
    name = string(Symbolics.getname(_unwrap(v)))
    return Symbol(last(split(name, "₊")))
end

"""
$(TYPEDSIGNATURES)

Split the unknowns and right-hand sides of `sys` into the infected subsystem and the
uninfected subsystem. Returns a named tuple
`(infected, uninfected, f_infected, f_uninfected)`.
"""
function infected_subsystem(sys::AbstractSystem, infected; autonomous::Bool = true)
    states, rhs = ode_right_hand_sides(sys; autonomous)
    x = resolve_states(sys, states, infected)
    xi = [_findsym(v, states) for v in x]
    yi = setdiff(eachindex(states), xi)
    return (
        infected = x, uninfected = states[yi], f_infected = rhs[xi], f_uninfected = rhs[yi])
end

"""
$(TYPEDSIGNATURES)

Compute the infection-free steady state of `sys`: set the infected states to zero and solve
the uninfected subsystem `g(y) = 0` symbolically. Only steady states that are unique and
determined by a linear system are found automatically; otherwise (e.g. an SIR model without
demography, for which every `S` is a steady state) an `ArgumentError` asks you to supply the
equilibrium explicitly via the `equilibrium` keyword of [`next_generation_matrix`](@ref).

Returns a `Dict{Num,Any}` mapping every unknown to its equilibrium value.
"""
function disease_free_equilibrium(sys::AbstractSystem, infected)
    sub = infected_subsystem(sys, infected)
    zero_inf = Dict{Num, Any}(v => 0 for v in sub.infected)
    g = [Symbolics.simplify(substitute(f, zero_inf)) for f in sub.f_uninfected]
    y = sub.uninfected
    isempty(y) && return zero_inf
    if any(f -> symbolic_iszero(f; numeric = false), g)
        throw(ArgumentError("the infection-free steady state of `$(nameof(sys))` is not unique " *
                            "(an uninfected equation is identically zero when the infected states vanish); " *
                            "pass it explicitly with `equilibrium = Dict(...)`"))
    end
    sol = try
        Symbolics.symbolic_linear_solve(g .~ 0, y)
    catch err
        nothing
    end
    if sol === nothing
        return _polynomial_equilibrium(sys, g, y, zero_inf)
    end
    sol = sol isa AbstractVector ? sol : [sol]
    eq = copy(zero_inf)
    for (v, s) in zip(y, sol)
        s = tidy(Num(s))
        depends_on(s, y) &&
            throw(ArgumentError("the infection-free steady state of `$(nameof(sys))` is not unique; " *
                                "pass it explicitly with `equilibrium = Dict(...)`"))
        eq[v] = s
    end
    return eq
end

# Fallback for uninfected subsystems that are polynomial but not linear in the uninfected
# states (logistic host growth, say): solve with `Symbolics.symbolic_solve`, which needs
# Nemo to be loaded, and keep the candidate steady states with all uninfected states
# non-zero. Exactly one such candidate is accepted.
function _polynomial_equilibrium(sys, g, y, zero_inf)
    hint = "pass it explicitly with `equilibrium = Dict(...)`; polynomial steady states can be " *
           "found automatically when Nemo (and, for several coupled unknowns, Groebner) is loaded"
    sols = try
        length(y) == 1 ? [Dict(y[1] => r) for r in Symbolics.symbolic_solve(g[1], y[1])] :
        Symbolics.symbolic_solve(g, y)
    catch err
        throw(ArgumentError("could not solve for the infection-free steady state of `$(nameof(sys))` " *
                            "symbolically ($(sprint(showerror, err))); $(hint)"))
    end
    (sols === nothing || isempty(sols)) &&
        throw(ArgumentError("no infection-free steady state of `$(nameof(sys))` was found; $(hint)"))
    candidates = Dict{Num, Any}[]
    for sol in sols
        sol isa AbstractDict || continue
        vals = [tidy(Num(get(sol, v, v))) for v in y]
        any(v -> _structural_zero(v) || depends_on(v, y), vals) && continue
        push!(candidates, merge(zero_inf, Dict{Num, Any}(zip(y, vals))))
    end
    length(candidates) == 1 && return candidates[1]
    isempty(candidates) &&
        throw(ArgumentError("every infection-free steady state of `$(nameof(sys))` has a vanishing uninfected compartment; $(hint)"))
    throw(ArgumentError("`$(nameof(sys))` has several infection-free steady states: $(candidates); $(hint)"))
end

"""
$(TYPEDSIGNATURES)

`true` if the additive term depends on an uninfected state, for example `β S I / N`. This
is the `:uninfected_dependence` strategy of [`next_generation_matrix`](@ref).
"""
depends_on_uninfected(term, uninfected) = depends_on(term, uninfected)

"""
$(TYPEDSIGNATURES)

`true` if the additive term is non-linear in the infected states, for example
`β (N - I) I / N` after `S` has been eliminated. This is the `:nonlinear_in_infected`
strategy of [`next_generation_matrix`](@ref).
"""
function nonlinear_in_infected(term, infected)
    depends_on(term, infected) || return false
    H = jacobian(jacobian([term], infected)[1, :], infected)
    return !all(x -> symbolic_iszero(x; numeric = false), H)
end

"""
$(TYPEDSIGNATURES)

Default classification of an additive term of an infected equation as a transmission
(new-infection) term: [`depends_on_uninfected`](@ref) or [`nonlinear_in_infected`](@ref).
Linear terms in the infected states alone (`σ E`, `-γ I`, `p μ I`) are transitions, so
vertical transmission written as an ODE term is *not* recognised; write such models as
reaction networks or give the transmission terms explicitly.
"""
function default_is_transmission(term, infected, uninfected)
    return depends_on_uninfected(term, uninfected) || nonlinear_in_infected(term, infected)
end

"""
$(TYPEDSIGNATURES)

Split the infected right-hand sides `f` into transmission and transition parts,
`f = F + Vpart`, term by term, using the predicate `is_transmission(term)`.
"""
function split_terms(f::AbstractVector, is_transmission)
    F = Num[Num(0) for _ in f]
    G = Num[Num(0) for _ in f]
    for (i, fi) in enumerate(f)
        for term in additive_terms(fi)
            if is_transmission(term)
                F[i] += term
            else
                G[i] += term
            end
        end
    end
    return F, G
end

"""
$(TYPEDSIGNATURES)

Linearise the split infected subsystem at `equilibrium`, giving the transmission matrix
`T = ∂F/∂x` and the transition matrix `Σ = ∂G/∂x`.
"""
function linearise(F::AbstractVector, G::AbstractVector, x::AbstractVector, equilibrium)
    T = tidy(substitute.(jacobian(F, x), Ref(equilibrium)))
    Σ = tidy(substitute.(jacobian(G, x), Ref(equilibrium)))
    return T, Σ
end
