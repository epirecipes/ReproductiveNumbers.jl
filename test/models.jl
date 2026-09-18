# Model definitions shared by several test files.
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
SIR model without demography: `S => N` must be supplied as the equilibrium.
"""
function sir_model()
    @parameters β γ N
    @variables S(t) I(t) R(t)
    eqs = [D(S) ~ -β * S * I / N,
        D(I) ~ β * S * I / N - γ * I,
        D(R) ~ γ * I]
    sys = complete(System(eqs, t; name = :sir))
    return (; sys, S, I, R, β, γ, N)
end

"""
SEIR model with births and deaths; the infection-free steady state is unique.
"""
function seir_demography_model()
    @parameters β σ γ μ N
    @variables S(t) E(t) I(t) R(t)
    eqs = [D(S) ~ μ * N - β * S * I / N - μ * S,
        D(E) ~ β * S * I / N - (σ + μ) * E,
        D(I) ~ σ * E - (γ + μ) * I,
        D(R) ~ γ * I - μ * R]
    sys = complete(System(eqs, t; name = :seir))
    return (; sys, S, E, I, R, β, σ, γ, μ, N)
end

"""
section 2.1 of Diekmann et al. (2010): SEI with two latent categories.
"""
function two_latent_model()
    @parameters β p ν₁ ν₂ γ μ N
    @variables S(t) E₁(t) E₂(t) I(t) R(t)
    eqs = [D(S) ~ μ * N - β * S * I / N - μ * S,
        D(E₁) ~ p * β * S * I / N - (ν₁ + μ) * E₁,
        D(E₂) ~ (1 - p) * β * S * I / N - (ν₂ + μ) * E₂,
        D(I) ~ ν₁ * E₁ + ν₂ * E₂ - (γ + μ) * I,
        D(R) ~ γ * I - μ * R]
    sys = complete(System(eqs, t; name = :two_latent))
    return (; sys, S, E₁, E₂, I, R, β, p, ν₁, ν₂, γ, μ, N)
end

"""
section 2.2 of Diekmann et al. (2010): SEI with two host categories.
"""
function two_host_model()
    @parameters β₁₁ β₁₂ β₂₁ β₂₂ ν₁ ν₂ γ₁ γ₂ μ N₁ N₂
    @variables S₁(t) S₂(t) E₁(t) E₂(t) I₁(t) I₂(t)
    N = N₁ + N₂
    eqs = [D(S₁) ~ μ * N₁ - β₁₁ * S₁ * I₁ / N - β₁₂ * S₁ * I₂ / N - μ * S₁,
        D(S₂) ~ μ * N₂ - β₂₁ * S₂ * I₁ / N - β₂₂ * S₂ * I₂ / N - μ * S₂,
        D(E₁) ~ β₁₁ * S₁ * I₁ / N + β₁₂ * S₁ * I₂ / N - (ν₁ + μ) * E₁,
        D(E₂) ~ β₂₁ * S₂ * I₁ / N + β₂₂ * S₂ * I₂ / N - (ν₂ + μ) * E₂,
        D(I₁) ~ ν₁ * E₁ - (γ₁ + μ) * I₁,
        D(I₂) ~ ν₂ * E₂ - (γ₂ + μ) * I₂]
    sys = complete(System(eqs, t; name = :two_host))
    return (; sys, S₁, S₂, E₁, E₂, I₁, I₂, β₁₁, β₁₂, β₂₁, β₂₂, ν₁, ν₂, γ₁, γ₂, μ, N₁, N₂)
end

"""
Ross–Macdonald malaria model (host `H`, vector `V`).
"""
function ross_macdonald_model()
    @parameters a b c γ μᵥ N_H N_V
    @variables S_H(t) I_H(t) S_V(t) I_V(t)
    eqs = [D(S_H) ~ -a * b * S_H * I_V / N_H + γ * I_H,
        D(I_H) ~ a * b * S_H * I_V / N_H - γ * I_H,
        D(S_V) ~ μᵥ * N_V - a * c * S_V * I_H / N_H - μᵥ * S_V,
        D(I_V) ~ a * c * S_V * I_H / N_H - μᵥ * I_V]
    sys = complete(System(eqs, t; name = :ross_macdonald))
    return (; sys, S_H, I_H, S_V, I_V, a, b, c, γ, μᵥ, N_H, N_V)
end

"""
Two non-interacting strains: `R₀` is the maximum of the strain-specific values.
"""
function two_strain_model()
    @parameters β₁ β₂ γ₁ γ₂ N
    @variables S(t) I₁(t) I₂(t) R(t)
    eqs = [D(S) ~ -β₁ * S * I₁ / N - β₂ * S * I₂ / N,
        D(I₁) ~ β₁ * S * I₁ / N - γ₁ * I₁,
        D(I₂) ~ β₂ * S * I₂ / N - γ₂ * I₂,
        D(R) ~ γ₁ * I₁ + γ₂ * I₂]
    sys = complete(System(eqs, t; name = :two_strain))
    return (; sys, S, I₁, I₂, R, β₁, β₂, γ₁, γ₂, N)
end

"""
SIS model with `S` eliminated: the infection term is non-linear in `I` alone.
"""
function sis_eliminated_model()
    @parameters β γ N
    @variables I(t)
    eqs = [D(I) ~ β * (N - I) * I / N - γ * I]
    sys = complete(System(eqs, t; name = :sis))
    return (; sys, I, β, γ, N)
end

"""
Vertical transmission: infected mothers bear infected juveniles (`J`).
"""
function vertical_transmission_model()
    @parameters β μ ν γ p N
    @variables S(t) J(t) I(t)
    eqs = [D(S) ~ μ * N - p * μ * I - β * S * I / N - μ * S,
        D(J) ~ p * μ * I - (ν + μ) * J,
        D(I) ~ β * S * I / N + ν * J - (γ + μ) * I]
    sys = complete(System(eqs, t; name = :vertical))
    return (; sys, S, J, I, β, μ, ν, γ, p, N)
end

"""
Random parameter values in `(lo, hi)` for a list of symbolic parameters.
"""
function random_parameters(rng, params; lo = 0.1, hi = 2.0)
    return Dict(p => lo + (hi - lo) * rand(rng) for p in params)
end
