using ReproductiveNumbers, ModelingToolkit, Symbolics, Test
using ReproductiveNumbers: default_is_transmission, split_terms, linearise,
                           ode_right_hand_sides, resolve_states
using ModelingToolkit: t_nounits as t, D_nounits as D
include(joinpath(@__DIR__, "..", "models.jl"))

@testset "ode_right_hand_sides" begin
    m = seir_demography_model()
    states, rhs = ode_right_hand_sides(m.sys)
    @test length(states) == 4 && length(rhs) == 4
    i = findfirst(s -> isequal(s, m.E), states)
    @test symbolic_isequal(rhs[i], m.β * m.S * m.I / m.N - (m.σ + m.μ) * m.E)

    # non-autonomous systems are rejected
    @parameters β γ N
    @variables S(t) I(t)
    eqs = [
        D(S) ~ -β * (1 + sin(t)) * S * I / N, D(I) ~ β * (1 + sin(t)) * S * I / N - γ * I]
    sys = complete(System(eqs, t; name = :forced))
    @test_throws ArgumentError ode_right_hand_sides(sys)
end

@testset "algebraic equations are compiled away" begin
    @parameters β γ
    @variables S(t) I(t) R(t) N(t)
    eqs = [D(S) ~ -β * S * I / N, D(I) ~ β * S * I / N - γ * I, D(R) ~ γ * I, N ~ S + I + R]
    sys = complete(System(eqs, t; name = :sir_obs))
    ngm = next_generation_matrix(sys, [I]; equilibrium = Dict(S => 1000, R => 0))
    @test symbolic_isequal(basic_reproduction_number(ngm), β / γ)
end

@testset "resolve_states" begin
    m = seir_demography_model()
    states, _ = ode_right_hand_sides(m.sys)
    @test isequal(resolve_states(m.sys, states, [:E, :I]), [m.E, m.I])
    @test isequal(resolve_states(m.sys, states, [m.I, m.E]), [m.I, m.E])
    @test_throws ArgumentError resolve_states(m.sys, states, [:Z])
    @test_throws ArgumentError resolve_states(m.sys, states, [m.E, m.E])
end

@testset "default_is_transmission" begin
    m = seir_demography_model()
    inf = [m.E, m.I]
    uninf = [m.S, m.R]
    @test default_is_transmission(m.β * m.S * m.I / m.N, inf, uninf)
    @test !default_is_transmission(m.σ * m.E, inf, uninf)
    @test !default_is_transmission(-(m.γ + m.μ) * m.I, inf, uninf)
    @test !default_is_transmission(m.μ * m.N, inf, uninf)           # constant inflow
    # non-linear in the infected states alone (S eliminated)
    @test default_is_transmission(m.β * (m.N - m.I) * m.I / m.N, [m.I], Num[])
    @test !default_is_transmission(-m.γ * m.I, [m.I], Num[])
end

@testset "split_terms and linearise" begin
    m = seir_demography_model()
    x = [m.E, m.I]
    f = [m.β * m.S * m.I / m.N - (m.σ + m.μ) * m.E, m.σ * m.E - (m.γ + m.μ) * m.I]
    F, G = split_terms(f, term -> default_is_transmission(term, x, [m.S, m.R]))
    @test symbolic_isequal(F[1], m.β * m.S * m.I / m.N)
    @test symbolic_isequal(F[2], 0)
    @test symbolic_isequal(G[1], -(m.σ + m.μ) * m.E)
    T, Σ = linearise(F, G, x, Dict(m.S => m.N, m.R => 0, m.E => 0, m.I => 0))
    @test symbolic_isequal(T[1, 2], m.β)
    @test symbolic_isequal(T[1, 1], 0) && symbolic_isequal(T[2, 1], 0) &&
          symbolic_isequal(T[2, 2], 0)
    @test symbolic_isequal(Σ[1, 1], -(m.σ + m.μ))
    @test symbolic_isequal(Σ[2, 1], m.σ)
    @test symbolic_isequal(Σ[2, 2], -(m.γ + m.μ))
    @test symbolic_isequal(Σ[1, 2], 0)
end

@testset "disease_free_equilibrium" begin
    m = seir_demography_model()
    dfe = disease_free_equilibrium(m.sys, [m.E, m.I])
    @test symbolic_isequal(dfe[m.S], m.N)
    @test symbolic_isequal(dfe[m.R], 0)
    @test dfe[m.E] == 0 && dfe[m.I] == 0
    # SIR without demography has no unique infection-free steady state
    s = sir_model()
    @test_throws ArgumentError disease_free_equilibrium(s.sys, [s.I])
    @test_throws ArgumentError next_generation_matrix(s.sys, [s.I])
    # all states infected: empty uninfected subsystem
    e = sis_eliminated_model()
    dfe2 = disease_free_equilibrium(e.sys, [e.I])
    @test dfe2[e.I] == 0
end

@testset "infected_subsystem" begin
    m = seir_demography_model()
    sub = infected_subsystem(m.sys, [:E, :I])
    @test isequal(sub.infected, [m.E, m.I])
    @test length(sub.uninfected) == 2
    @test length(sub.f_infected) == 2 && length(sub.f_uninfected) == 2
end
