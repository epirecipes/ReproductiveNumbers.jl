module ReproductiveNumbersCatalystExt

using ReproductiveNumbers
using ReproductiveNumbers: assemble, linearise, resolve_states, _equilibrium, _findsym,
                           tidy, ode_right_hand_sides, disease_free_equilibrium,
                           _warn_negative_transmissions, TERM_STRATEGIES
using Catalyst
using Catalyst: ReactionSystem, Reaction, reactions, species, oderatelaw, netstoichmat,
                nonreactions
using Symbolics: Num
import ReproductiveNumbers: next_generation_matrix, infected_subsystem,
                            disease_free_equilibrium, reaction_is_transmission

function ReproductiveNumbers.reaction_is_transmission(k::Int, infected_idx, ν)
    return sum(ν[i, k] for i in infected_idx; init = 0) > 0
end

_ode_system(rn::ReactionSystem) = Catalyst.ode_model(rn)

function infected_subsystem(rn::ReactionSystem, infected)
    return infected_subsystem(_ode_system(rn), infected)
end
function disease_free_equilibrium(rn::ReactionSystem, infected)
    return disease_free_equilibrium(_ode_system(rn), infected)
end

function next_generation_matrix(rn::ReactionSystem, infected;
        equilibrium = nothing, transmission = :stoichiometry, warn::Bool = true,
        autonomous::Bool = true,
        combinatoric_ratelaws = Catalyst.get_combinatoric_ratelaws(rn))
    sys = _ode_system(rn)
    # term-based strategies, and networks with coupled non-reaction equations, go
    # through the ODE route
    if transmission isa Symbol && transmission in TERM_STRATEGIES
        return next_generation_matrix(
            sys, infected; equilibrium, transmission, warn, autonomous)
    end
    if !isempty(nonreactions(rn))
        warn &&
            @warn "the reaction network has coupled non-reaction equations; transmissions are " *
                  "identified from the ODE terms (`:auto`) rather than from the stoichiometry"
        return next_generation_matrix(
            sys, infected; equilibrium, transmission = :auto, warn)
    end
    states, _ = ode_right_hand_sides(sys; autonomous)
    x = resolve_states(sys, states, infected)
    xi = [_findsym(v, states) for v in x]
    y = states[setdiff(eachindex(states), xi)]
    eq = _equilibrium(sys, infected, equilibrium, x, y)
    rxs = reactions(rn)
    spcs = Num.(species(rn))
    ν = netstoichmat(rn)
    sp_idx = [_findsym(v, spcs) for v in x]
    any(isnothing, sp_idx) &&
        throw(ArgumentError("infected states must be species of the reaction network"))
    istrans, method = _transmission_flags(rxs, transmission, sp_idx, ν)
    F = Num[Num(0) for _ in x]
    G = Num[Num(0) for _ in x]
    for (k, rx) in enumerate(rxs)
        rate = Num(oderatelaw(rx; combinatoric_ratelaw = combinatoric_ratelaws))
        for (i, si) in enumerate(sp_idx)
            c = ν[si, k]
            iszero(c) && continue
            if istrans[k]
                F[i] += c * rate
            else
                G[i] += c * rate
            end
        end
    end
    T, Σ = linearise(F, G, x, eq)
    warn && _warn_negative_transmissions(T, x, F)
    return assemble(x, y, eq, T, Σ; F, G, method)
end

function _transmission_flags(rxs, transmission, sp_idx, ν)
    n = length(rxs)
    if transmission === nothing || transmission === :stoichiometry
        return [reaction_is_transmission(k, sp_idx, ν) for k in 1:n], :stoichiometry
    elseif transmission isa Symbol
        throw(ArgumentError("unknown transmission strategy `$(repr(transmission))` for a reaction network; use `:stoichiometry`, one of $(TERM_STRATEGIES), a predicate on reactions, a `Vector{Bool}` or reaction indices"))
    elseif transmission isa Function
        return [Bool(transmission(rx)) for rx in rxs], :predicate
    elseif transmission isa AbstractVector{Bool}
        length(transmission) == n ||
            throw(ArgumentError("`transmission` must have one flag per reaction ($n)"))
        return collect(transmission), :explicit
    elseif transmission isa AbstractVector{<:Integer}
        all(k -> 1 <= k <= n, transmission) ||
            throw(ArgumentError("reaction indices must lie in 1:$n"))
        flags = falses(n)
        flags[transmission] .= true
        return flags, :explicit
    else
        throw(ArgumentError("`transmission` must be `:stoichiometry`, a term strategy, a predicate on reactions, a `Vector{Bool}` or reaction indices"))
    end
end
end
