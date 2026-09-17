using ReproductiveNumbers, ModelingToolkit, Symbolics, LinearAlgebra, Random, Test
include(joinpath(@__DIR__, "..", "models.jl"))
rng = MersenneTwister(2)

@testset "SIR" begin
    m = sir_model()
    ngm = next_generation_matrix(m.sys, [m.I]; equilibrium = Dict(m.S => m.N, m.R => 0))
    @test ngm isa NextGenerationMatrix
    @test size(ngm.T) == (1, 1) && size(ngm.K) == (1, 1)
    @test symbolic_isequal(ngm.T[1, 1], m.β)
    @test symbolic_isequal(ngm.Σ[1, 1], -m.γ)
    @test symbolic_isequal(ngm.K_L[1, 1], m.β / m.γ)
    @test isequal(states_at_infection(ngm), [m.I])
    @test symbolic_isequal(basic_reproduction_number(ngm), m.β / m.γ)
    @test symbolic_isequal(
        basic_reproduction_number(m.sys, [:I]; equilibrium = [m.S => m.N]), m.β / m.γ)
    # partial substitution: effective reproduction number at general S
    ngm_s = next_generation_matrix(m.sys, [m.I]; equilibrium = Dict(m.R => 0))
    @test symbolic_isequal(basic_reproduction_number(ngm_s), m.β * m.S / (m.γ * m.N))
    # unknown symbol in equilibrium
    @test_throws ArgumentError next_generation_matrix(
        m.sys, [m.I]; equilibrium = Dict(m.β => 1))
    # show method
    str = sprint(show, MIME("text/plain"), ngm)
    @test occursin("NextGenerationMatrix (symbolic)", str)
end

@testset "SEIR with demography" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    @test isequal(ngm.infected, [m.E, m.I])
    @test isequal(states_at_infection(ngm), [m.E])
    @test ngm.E == reshape([1, 0], 2, 1)
    @test size(ngm.K_L) == (2, 2) && size(ngm.K) == (1, 1)
    R0 = m.β * m.σ / ((m.σ + m.μ) * (m.γ + m.μ))
    @test symbolic_isequal(ngm.K_L[1, 1], R0)
    @test symbolic_isequal(ngm.K_L[1, 2], m.β / (m.γ + m.μ))
    @test symbolic_isequal(ngm.K_L[2, 1], 0) && symbolic_isequal(ngm.K_L[2, 2], 0)
    @test symbolic_isequal(basic_reproduction_number(ngm), R0)
    # order of infected states does not change R0
    ngm2 = next_generation_matrix(m.sys, [:I, :E])
    @test symbolic_isequal(basic_reproduction_number(ngm2), R0)
    # the equilibrium was found automatically
    @test symbolic_isequal(ngm.equilibrium[m.S], m.N)
    # explicit transmission vector reproduces the default
    ngm3 = next_generation_matrix(
        m.sys, [m.E, m.I]; transmission = [m.β * m.S * m.I / m.N, 0])
    @test symbolic_isequal(basic_reproduction_number(ngm3), R0)
    # a predicate on terms
    ngm4 = next_generation_matrix(
        m.sys, [m.E, m.I]; transmission = term -> occursin("S(t)", string(term)))
    @test symbolic_isequal(basic_reproduction_number(ngm4), R0)
    # wrong-length transmission vector
    @test_throws ArgumentError next_generation_matrix(m.sys, [m.E, m.I]; transmission = [0])
    @test_throws ArgumentError next_generation_matrix(
        m.sys, [m.E, m.I]; transmission = :bogus)
end

@testset "Two latent categories (rank one)" begin
    m = two_latent_model()
    ngm = next_generation_matrix(m.sys, [m.E₁, m.E₂, m.I])
    @test isequal(states_at_infection(ngm), [m.E₁, m.E₂])
    @test size(ngm.K) == (2, 2)
    k(a, b, c) = a * b * m.β / (c * (m.γ + m.μ))
    @test symbolic_isequal(ngm.K[1, 1], k(m.p, m.ν₁, m.ν₁ + m.μ))
    @test symbolic_isequal(ngm.K[1, 2], k(m.p, m.ν₂, m.ν₂ + m.μ))
    @test symbolic_isequal(ngm.K[2, 1], k(1 - m.p, m.ν₁, m.ν₁ + m.μ))
    @test symbolic_isequal(ngm.K[2, 2], k(1 - m.p, m.ν₂, m.ν₂ + m.μ))
    @test symbolic_isequal(det(ngm.K), 0)
    R0 = k(m.p, m.ν₁, m.ν₁ + m.μ) + k(1 - m.p, m.ν₂, m.ν₂ + m.μ)
    @test symbolic_isequal(basic_reproduction_number(ngm), R0)
    K_S = small_domain_matrix(ngm)
    @test size(K_S) == (1, 1)
    @test symbolic_isequal(K_S[1, 1], R0)
end

@testset "Two host categories (2x2 quadratic formula)" begin
    m = two_host_model()
    ngm = next_generation_matrix(m.sys, [m.E₁, m.E₂, m.I₁, m.I₂])
    @test isequal(states_at_infection(ngm), [m.E₁, m.E₂])
    N = m.N₁ + m.N₂
    k(β, νj, νjμ, γi) = β * (m.N₁ + m.N₂ - N) # placeholder never used
    K11 = m.β₁₁ * m.N₁ / N * m.ν₁ / ((m.ν₁ + m.μ) * (m.γ₁ + m.μ))
    K12 = m.β₁₂ * m.N₁ / N * m.ν₂ / ((m.ν₂ + m.μ) * (m.γ₂ + m.μ))
    K21 = m.β₂₁ * m.N₂ / N * m.ν₁ / ((m.ν₁ + m.μ) * (m.γ₁ + m.μ))
    K22 = m.β₂₂ * m.N₂ / N * m.ν₂ / ((m.ν₂ + m.μ) * (m.γ₂ + m.μ))
    @test symbolic_isequal(ngm.K[1, 1], K11)
    @test symbolic_isequal(ngm.K[1, 2], K12)
    @test symbolic_isequal(ngm.K[2, 1], K21)
    @test symbolic_isequal(ngm.K[2, 2], K22)
    R0 = basic_reproduction_number(ngm)
    expected = (K11 + K22 + sqrt((K11 - K22)^2 + 4 * K12 * K21)) / 2
    for _ in 1:3
        vals = random_parameters(
            rng, [m.β₁₁, m.β₁₂, m.β₂₁, m.β₂₂, m.ν₁, m.ν₂, m.γ₁, m.γ₂, m.μ, m.N₁, m.N₂])
        @test ReproductiveNumbers.to_number(substitute(R0, vals)) ≈
              ReproductiveNumbers.to_number(substitute(expected, vals))
        @test basic_reproduction_number(ngm, vals) ≈
              ReproductiveNumbers.to_number(substitute(expected, vals))
    end
    # separable mixing gives rank one and R0 = trace
    sep = Dict(m.β₁₁ => 1.0, m.β₁₂ => 2.0, m.β₂₁ => 1.5, m.β₂₂ => 3.0)
    Ksep = substitute.(ngm.K, Ref(sep))
    @test ReproductiveNumbers.isrankone(Ksep)
    @test symbolic_isequal(spectral_radius(Ksep), tr(Ksep))
end

@testset "Ross-Macdonald (vector-host)" begin
    m = ross_macdonald_model()
    ngm = next_generation_matrix(
        m.sys, [m.I_H, m.I_V]; equilibrium = Dict(m.S_H => m.N_H, m.S_V => m.N_V))
    @test isequal(states_at_infection(ngm), [m.I_H, m.I_V])
    @test symbolic_isequal(ngm.K[1, 1], 0) && symbolic_isequal(ngm.K[2, 2], 0)
    @test symbolic_isequal(ngm.K[1, 2], m.a * m.b / m.μᵥ)
    @test symbolic_isequal(ngm.K[2, 1], m.a * m.c * m.N_V / (m.N_H * m.γ))
    R0 = basic_reproduction_number(ngm)
    @test symbolic_isequal(R0^2, m.a^2 * m.b * m.c * m.N_V / (m.N_H * m.γ * m.μᵥ))
    # host-to-host generations: count only vector->host as transmission
    ngm_h = next_generation_matrix(m.sys, [m.I_H, m.I_V];
        equilibrium = Dict(m.S_H => m.N_H, m.S_V => m.N_V),
        transmission = [m.a * m.b * m.S_H * m.I_V / m.N_H, 0])
    @test isequal(states_at_infection(ngm_h), [m.I_H])
    @test symbolic_isequal(basic_reproduction_number(ngm_h), R0^2)
end

@testset "Two strains (block structure -> max)" begin
    m = two_strain_model()
    ngm = next_generation_matrix(
        m.sys, [m.I₁, m.I₂]; equilibrium = Dict(m.S => m.N, m.R => 0))
    R0 = basic_reproduction_number(ngm)
    @test isequal(R0, max(m.β₁ / m.γ₁, m.β₂ / m.γ₂)) ||
          isequal(R0, max(m.β₂ / m.γ₂, m.β₁ / m.γ₁))
    @test basic_reproduction_number(
        ngm, Dict(m.β₁ => 2, m.β₂ => 3, m.γ₁ => 1, m.γ₂ => 2, m.N => 1)) ≈ 2.0
end

@testset "SIS with S eliminated" begin
    m = sis_eliminated_model()
    ngm = next_generation_matrix(m.sys, [m.I])
    @test symbolic_isequal(ngm.T[1, 1], m.β)
    @test symbolic_isequal(ngm.Σ[1, 1], -m.γ)
    @test symbolic_isequal(basic_reproduction_number(ngm), m.β / m.γ)
end

@testset "Vertical transmission needs an explicit transmission choice" begin
    m = vertical_transmission_model()
    # default: births of infected juveniles are classified as transitions
    ngm_default = next_generation_matrix(m.sys, [m.J, m.I])
    @test isequal(states_at_infection(ngm_default), [m.I])
    # explicit: vertical transmission counts as a new infection
    ngm = next_generation_matrix(m.sys, [m.J, m.I];
        transmission = [m.p * m.μ * m.I, m.β * m.S * m.I / m.N])
    @test isequal(states_at_infection(ngm), [m.J, m.I])
    @test symbolic_isequal(ngm.T[1, 2], m.p * m.μ)
    @test symbolic_isequal(ngm.T[2, 2], m.β)
    # both decompositions agree on the threshold
    for _ in 1:5
        vals = random_parameters(rng, [m.β, m.μ, m.ν, m.γ, m.N]; lo = 0.05, hi = 1.5)
        vals[m.p] = rand(rng)
        r_default = basic_reproduction_number(ngm_default, vals)
        r_explicit = basic_reproduction_number(ngm, vals)
        @test (r_default > 1) == (r_explicit > 1)
        @test validate_decomposition(ngm, vals)
        @test validate_decomposition(ngm_default, vals)
    end
end

@testset "Direct construction from T and Σ" begin
    @variables β γ σ
    T = Num[0 β; 0 0]
    Σ = Num[-σ 0; σ -γ]
    ngm = next_generation_matrix(T, Σ; infected = [:E, :I])
    @test symbolic_isequal(basic_reproduction_number(ngm), β / γ)
    ngm_num = next_generation_matrix([0 2.0; 0 0], [-1.0 0; 1.0 -0.5])
    @test basic_reproduction_number(ngm_num) ≈ 4.0
    @test_throws DimensionMismatch next_generation_matrix(Num[0 β], Σ)
    str = sprint(show, MIME("text/plain"), ngm_num)
    @test occursin("numeric", str)
end
