module ReproductiveNumbersCatalystExt

using ReproductiveNumbers
using ReproductiveNumbers: assemble, linearise, resolve_states, _equilibrium, _findsym,
                           tidy,
                           ode_right_hand_sides, disease_free_equilibrium
using Catalyst
using Catalyst: ReactionSystem, Reaction, reactions, species, oderatelaw, netstoichmat
using Symbolics: Num
import ReproductiveNumbers: next_generation_matrix, infected_subsystem,
                            disease_free_equilibrium

"""
    default_is_transmission(rx::Reaction, infected_indices, netstoich)

A reaction is a transmission (creates new infected individuals) if its net stoichiometry
increases the total number of individuals in the infected states.
"""
function _reaction_is_transmission(k::Int, infected_idx, ν)
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
        equilibrium = nothing, transmission = nothing,
        combinatoric_ratelaws = Catalyst.get_combinatoric_ratelaws(rn))
    sys = _ode_system(rn)
    states, _ = ode_right_hand_sides(sys)
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
    istrans = _transmission_flags(rxs, transmission, sp_idx, ν)

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
    return assemble(x, y, eq, T, Σ)
end

function _transmission_flags(rxs, transmission, sp_idx, ν)
    n = length(rxs)
    if transmission === nothing
        return [_reaction_is_transmission(k, sp_idx, ν) for k in 1:n]
    elseif transmission isa Function
        return [Bool(transmission(rx)) for rx in rxs]
    elseif transmission isa AbstractVector{Bool}
        length(transmission) == n ||
            throw(ArgumentError("`transmission` must have one flag per reaction ($n)"))
        return collect(transmission)
    elseif transmission isa AbstractVector{<:Integer}
        all(k -> 1 <= k <= n, transmission) ||
            throw(ArgumentError("reaction indices must lie in 1:$n"))
        flags = falses(n)
        flags[transmission] .= true
        return flags
    else
        throw(ArgumentError("`transmission` must be `nothing`, a predicate on reactions, a `Vector{Bool}` or reaction indices"))
    end
end

end
