using ReproductiveNumbers, ModelingToolkit, Catalyst, Symbolics, LinearAlgebra, Random, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
include(joinpath(@__DIR__, "..", "models.jl"))
rng = MersenneTwister(12)

@testset "transmission strategies are recorded" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    @test ngm.method == :auto
    tm = transmission_method(ngm)
    @test tm.method == :auto
    @test symbolic_isequal(tm.F[1], m.β * m.S * m.I / m.N)
    @test symbolic_isequal(tm.G[2], m.σ * m.E - (m.γ + m.μ) * m.I)
    str = sprint(show, MIME("text/plain"), ngm)
    @test occursin("transmissions identified by: :auto", str)
    @test occursin("new infections in E(t)", str)
    for strategy in (:uninfected_dependence, :nonlinear_in_infected)
        n2 = next_generation_matrix(m.sys, [m.E, m.I]; transmission = strategy)
        @test n2.method == strategy
    end
    # only the uninfected-dependence rule finds βSI/N; only the non-linearity rule
    # finds β(N-I)I/N
    @test symbolic_isequal(
        basic_reproduction_number(m.sys, [m.E, m.I]; transmission = :uninfected_dependence),
        m.β * m.σ / ((m.σ + m.μ) * (m.γ + m.μ)))
    @test symbolic_isequal(
        basic_reproduction_number(m.sys, [m.E, m.I]; transmission = :nonlinear_in_infected),
        0)
    e = sis_eliminated_model()
    @test symbolic_isequal(
        basic_reproduction_number(e.sys, [e.I]; transmission = :nonlinear_in_infected),
        e.β / e.γ)
    @test symbolic_isequal(
        basic_reproduction_number(e.sys, [e.I]; transmission = :uninfected_dependence), 0)
    @test_throws ArgumentError next_generation_matrix(
        m.sys, [m.E, m.I]; transmission = :bogus)
    @test next_generation_matrix(m.sys, [m.E, m.I];
        transmission = term -> occursin("S(t)", string(term))).method == :predicate
    @test next_generation_matrix(
        m.sys, [m.E, m.I]; transmission = [m.β * m.S * m.I / m.N, 0]).method == :explicit
    # evaluated matrices keep the record
    num = evaluate(ngm, Dict(m.β => 0.5, m.σ => 0.25, m.γ => 0.2, m.μ => 0.01, m.N => 1.0))
    @test num.method == :auto && length(num.F) == 2
    @test next_generation_matrix([0 1.0; 0 0], [-1.0 0; 1 -1]).method == :matrices
end

@testset "negative transmission entries are reported" begin
    # density-dependent death of infecteds involves S, so the default rule misfiles it
    @parameters β γ d N
    @variables S(t) I(t)
    eqs = [D(S) ~ -β * S * I / N, D(I) ~ β * S * I / N - γ * I - d * (S + I) * I]
    sys = complete(System(eqs, t; name = :dd))
    @test_logs (:warn, r"negative") next_generation_matrix(
        sys, [I]; equilibrium = Dict(S => N))
    @test_logs next_generation_matrix(sys, [I]; equilibrium = Dict(S => N), warn = false)
    ngm = next_generation_matrix(sys, [I]; equilibrium = Dict(S => N),
        transmission = [β * S * I / N])
    @test symbolic_isequal(basic_reproduction_number(ngm), β / (γ + d * N))
end

@testset "ambiguous symbol names are rejected" begin
    @parameters β γ
    @variables I(t)
    @named a = System([D(I) ~ -γ * I], t)
    @named b = System([D(I) ~ -γ * I], t)
    @variables S(t)
    sys = complete(compose(System([D(S) ~ -β * S], t; name = :top), [a, b]))
    @test_throws ArgumentError next_generation_matrix(sys, [:I]; equilibrium = Dict(S => 1))
end

@testset "polynomial infection-free steady state with Nemo" begin
    import Nemo   # `using` would shadow `evaluate` in Main
    @parameters r K μ β γ
    @variables S(t) I(t)
    eqs = [D(S) ~ r * S * (1 - S / K) - μ * S - β * S * I,
        D(I) ~ β * S * I - (γ + μ) * I]
    sys = complete(System(eqs, t; name = :logistic))
    dfe = disease_free_equilibrium(sys, [I])
    @test symbolic_isequal(dfe[S], K * (1 - μ / r))
    @test symbolic_isequal(
        basic_reproduction_number(sys, [I]), β * K * (1 - μ / r) / (γ + μ))
    # a quadratic with two admissible roots is reported
    @parameters c
    eqs2 = [D(S) ~ (S - 1) * (S - c), D(I) ~ β * S * I - γ * I]
    sys2 = complete(System(eqs2, t; name = :two_roots))
    @test_throws ArgumentError disease_free_equilibrium(sys2, [I])
end

@testset "observed equations with parameters" begin
    @parameters β γ N₀
    @variables S(t) I(t) R(t) N(t) λ(t)
    eqs = [D(S) ~ -λ * S, D(I) ~ λ * S - γ * I, D(R) ~ γ * I, N ~ S + I + R, λ ~ β * I / N]
    sys = complete(System(eqs, t; name = :obs))
    ngm = next_generation_matrix(sys, [I]; equilibrium = Dict(S => N₀, R => 0))
    @test symbolic_isequal(basic_reproduction_number(ngm), β / γ)
end

@testset "Catalyst: combinatoric rate laws and coupled equations" begin
    rn = @reaction_network pair begin
        @parameters β γ N
        β / N, S + 2I --> 3I
        γ, I --> R
    end
    @unpack S, I, R, β, γ, N = rn
    ngm = next_generation_matrix(rn, [I]; equilibrium = Dict(S => N, R => 0))
    # the transmission rate is quadratic in I, so it vanishes on linearisation
    @test symbolic_isequal(ngm.T[1, 1], 0)
    @test symbolic_isequal(ngm.Σ[1, 1], -γ)
    ngm2 = next_generation_matrix(
        rn, [I]; equilibrium = Dict(S => N, R => 0), combinatoric_ratelaws = false)
    @test symbolic_isequal(ngm2.T[1, 1], 0)
    # term strategies on a network go through the ODE route
    seir = @reaction_network seir begin
        @parameters β σ γ μ N
        β / N, S + I --> E + I
        σ, E --> I
        γ, I --> R
        μ * N, 0 --> S
        μ, (S, E, I, R) --> 0
    end
    @unpack E = seir
    n_auto = next_generation_matrix(seir, [:E, :I]; transmission = :auto)
    @test n_auto.method == :auto
    @test symbolic_isequal(
        basic_reproduction_number(n_auto), basic_reproduction_number(seir, [:E, :I]))
    @test_throws ArgumentError next_generation_matrix(seir, [:E, :I]; transmission = :bogus)
    # coupled non-reaction equation: falls back to the ODE route with a warning
    @parameters k
    @variables T(t)
    coupled = @network_component coupled begin
        @parameters β γ N k
        @variables T(t)
        @equations D(T) ~ -k * T
        β / N, S + I --> 2I
        γ, I --> R
    end
    coupled = complete(coupled)
    @unpack S, I, R, N = coupled
    @test_logs (:warn, r"coupled") next_generation_matrix(
        coupled, [I]; equilibrium = Dict(S => N, R => 0, T => 0))
    n_c = next_generation_matrix(
        coupled, [I]; equilibrium = Dict(S => N, R => 0, T => 0), warn = false)
    @test n_c.method == :auto
end

@testset "type reproduction number of a proper subset of size two" begin
    Kn = [0.5 0.2 0.1; 0.3 0.4 0.2; 0.1 0.1 0.3]
    ngm = next_generation_matrix(Num.(Kn), Num.(-Matrix(I, 3, 3)))
    ngm_num = ReproductiveNumbers.evaluate(ngm, Dict())
    P = Diagonal([1.0, 1.0, 0.0])
    Q = I - P
    M = P * Kn * inv(I - Q * Kn)
    @test type_reproduction_number(ngm_num, [1, 2]) ≈ maximum(abs, eigvals(M[1:2, 1:2]))
    @test type_reproduction_number(ngm_num, [1, 2, 3]) ≈ basic_reproduction_number(ngm_num)
    @test (type_reproduction_number(ngm_num, [1, 2]) > 1) ==
          (basic_reproduction_number(ngm_num) > 1)
end

@testset "Latexify extension" begin
    using Latexify
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    s = string(latexify(ngm))
    @test occursin("\\Sigma", s) && occursin("R_0", s) && occursin("begin{aligned}", s)
    @test occursin("Sigma", sprint(show, MIME("text/latex"), ngm))
    # no closed form: R0 line is omitted
    @variables a b c d e f g h k
    n3 = next_generation_matrix(Num[a b c; d e f; g h k], Num.(-Matrix(I, 3, 3)))
    @test !occursin("R_0", string(latexify(n3)))
end

@testset "accessors and error printing" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    @test isequal(infected_states(ngm), [m.E, m.I])
    T, Σ = transmission_transition_matrices(ngm)
    @test T === ngm.T && Σ === ngm.Σ
    @test occursin("NoClosedFormError: boom", sprint(showerror, NoClosedFormError("boom")))
    # numeric sign guard
    @test_logs (:warn, r"sign conventions") next_generation_matrix(
        [0 -1.0; 0 0], [-1.0 0; 1 -1])
    # non-invertible symbolic Σ
    @variables β
    @test_throws ArgumentError next_generation_matrix(Num[0 β; 0 0], Num[0 0; 0 0])
end
