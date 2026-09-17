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
algebraic (observed) equations are compiled with `mtkcompile` first.
"""
function ode_right_hand_sides(sys::AbstractSystem)
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
        depends_on(r, [iv]) &&
            throw(ArgumentError("the right-hand side `$(r)` depends explicitly on the independent variable; " *
                                "next-generation matrices require an autonomous system"))
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
            i = findfirst(s -> _basename(s) == x, states)
            i === nothing &&
                throw(ArgumentError("no unknown of `$(nameof(sys))` is named `$(x)`; unknowns are $(states)"))
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
function infected_subsystem(sys::AbstractSystem, infected)
    states, rhs = ode_right_hand_sides(sys)
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
        throw(ArgumentError("could not solve for the infection-free steady state of `$(nameof(sys))` " *
                            "symbolically ($(sprint(showerror, err))); pass it explicitly with `equilibrium = Dict(...)`"))
    end
    sol === nothing &&
        throw(ArgumentError("could not solve for the infection-free steady state of `$(nameof(sys))`; " *
                            "pass it explicitly with `equilibrium = Dict(...)`"))
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

"""
$(TYPEDSIGNATURES)

Default classification of an additive term of an infected equation as a transmission
(new-infection) term: it depends on an uninfected state (for example `β S I / N`), or it
is non-linear in the infected states (for example `β (N - I) I / N` after `S` has been
eliminated). Linear terms in the infected states alone (`σ E`, `-γ I`) are transitions.
"""
function default_is_transmission(term, infected, uninfected)
    depends_on(term, uninfected) && return true
    depends_on(term, infected) || return false
    H = jacobian(jacobian([term], infected)[1, :], infected)
    return !all(x -> symbolic_iszero(x; numeric = false), H)
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
