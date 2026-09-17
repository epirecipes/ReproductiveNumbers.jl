# Closed-form reproduction numbers from the literature, reproduced exactly.
using ReproductiveNumbers, ModelingToolkit, Symbolics, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
include(joinpath(@__DIR__, "..", "models.jl"))

@testset "Diekmann, Heesterbeek & Roberts (2010), Example 2.1" begin
    m = two_latent_model()
    R0 = basic_reproduction_number(m.sys, [m.E₁, m.E₂, m.I])
    paper = m.ν₁ * m.β * m.p / ((m.ν₁ + m.μ) * (m.γ + m.μ)) +
            (1 - m.p) * m.ν₂ * m.β / ((m.ν₂ + m.μ) * (m.γ + m.μ))
    @test symbolic_isequal(R0, paper)
end

@testset "Diekmann, Heesterbeek & Roberts (2010), Example 4.1 (two sexes)" begin
    @parameters β₁ β₂ ν₁ ν₂ γ₁ γ₂ μ N₁ N₂
    @variables S₁(t) E₁(t) I₁(t) S₂(t) E₂(t) I₂(t)
    # β₁: infections of females per infectious male per unit time when all females are
    # susceptible (and vice versa for β₂), which is the paper's normalisation.
    eqs = [D(S₁) ~ μ * N₁ - β₁ * S₁ * I₂ / N₁ - μ * S₁,
        D(E₁) ~ β₁ * S₁ * I₂ / N₁ - (ν₁ + μ) * E₁,
        D(I₁) ~ ν₁ * E₁ - (γ₁ + μ) * I₁,
        D(S₂) ~ μ * N₂ - β₂ * S₂ * I₁ / N₂ - μ * S₂,
        D(E₂) ~ β₂ * S₂ * I₁ / N₂ - (ν₂ + μ) * E₂,
        D(I₂) ~ ν₂ * E₂ - (γ₂ + μ) * I₂]
    sys = complete(System(eqs, t; name = :sti))
    ngm = next_generation_matrix(sys, [E₁, I₁, E₂, I₂])
    @test isequal(states_at_infection(ngm), [E₁, E₂])
    @test symbolic_isequal(ngm.K[1, 1], 0) && symbolic_isequal(ngm.K[2, 2], 0)
    @test symbolic_isequal(ngm.K[1, 2], β₁ * ν₂ / ((ν₂ + μ) * (γ₂ + μ)))
    @test symbolic_isequal(ngm.K[2, 1], β₂ * ν₁ / ((ν₁ + μ) * (γ₁ + μ)))
    R0 = basic_reproduction_number(ngm)
    paper = sqrt(β₁ * β₂ * ν₁ * ν₂ / ((ν₁ + μ) * (γ₁ + μ) * (ν₂ + μ) * (γ₂ + μ)))
    @test symbolic_isequal(R0^2, paper^2)
end

@testset "van den Driessche & Watmough (2002) treatment model" begin
    # Section 4.1 of vdD&W: SEIT tuberculosis model with treatment.
    @parameters β₁ β₂ d ν r₁ r₂ p q Λ N
    @variables S(t) E(t) I(t) T(t)
    eqs = [D(S) ~ Λ - d * S - β₁ * S * I / N,
        D(E) ~ β₁ * S * I / N + β₂ * T * I / N - (d + ν + r₁) * E + p * r₂ * I,
        D(I) ~ ν * E - (d + r₂) * I,
        D(T) ~ -d * T + r₁ * E + q * r₂ * I - β₂ * T * I / N]
    sys = complete(System(eqs, t; name = :tb))
    # DFE: S = Λ/d, T = 0, N = S at the DFE (N is treated as the equilibrium population)
    ngm = next_generation_matrix(sys, [E, I]; equilibrium = Dict(S => N, T => 0))
    @test symbolic_isequal(ngm.T[1, 2], β₁)
    @test symbolic_isequal(ngm.Σ[1, 1], -(d + ν + r₁))
    @test symbolic_isequal(ngm.Σ[1, 2], p * r₂)
    @test symbolic_isequal(ngm.Σ[2, 1], ν)
    @test symbolic_isequal(ngm.Σ[2, 2], -(d + r₂))
    paper = β₁ * ν / ((d + ν + r₁) * (d + r₂) - ν * p * r₂)
    @test symbolic_isequal(basic_reproduction_number(ngm), paper)
end

@testset "Ross-Macdonald" begin
    m = ross_macdonald_model()
    R0 = basic_reproduction_number(m.sys, [m.I_H, m.I_V];
        equilibrium = Dict(m.S_H => m.N_H, m.S_V => m.N_V))
    @test symbolic_isequal(R0^2, m.a^2 * m.b * m.c * (m.N_V / m.N_H) / (m.γ * m.μᵥ))
end

@testset "SEIR with demography (Anderson & May)" begin
    m = seir_demography_model()
    R0 = basic_reproduction_number(m.sys, [m.E, m.I])
    @test symbolic_isequal(R0, m.β * m.σ / ((m.σ + m.μ) * (m.γ + m.μ)))
end
