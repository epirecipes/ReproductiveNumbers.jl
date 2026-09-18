# Models of sections 4.2 and 4.3 of Diekmann, Heesterbeek & Roberts (2010) and the staged
# progression model of van den Driessche & Watmough (2002), checked for structure and
# internal consistency (the exact ODE systems are not reproduced in the paper, so the
# entries are not compared term by term).
using ReproductiveNumbers, ModelingToolkit, Symbolics, LinearAlgebra, Random, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
rng = MersenneTwister(21)
nz(A) = sort(abs.(filter(x -> abs(x) > 1e-9, eigvals(A))))

@testset "Section 4.2: sexually transmitted infection with vertical transmission" begin
    @parameters p μ ν₁ ν₂ γ₁ γ₂ β₁ β₂ N₁ N₂
    @variables S₁(t) J₁(t) I₁(t) S₂(t) J₂(t) I₂(t)
    eqs = [D(S₁) ~ μ * N₁ - p * μ * I₁ - β₁ * S₁ * I₂ / N₂ - μ * S₁,
        D(J₁) ~ p * μ * I₁ - (ν₁ + μ) * J₁,
        D(I₁) ~ ν₁ * J₁ + β₁ * S₁ * I₂ / N₂ - (γ₁ + μ) * I₁,
        D(S₂) ~ μ * N₂ - p * μ * I₂ - β₂ * S₂ * I₁ / N₁ - μ * S₂,
        D(J₂) ~ p * μ * I₂ - (ν₂ + μ) * J₂,
        D(I₂) ~ ν₂ * J₂ + β₂ * S₂ * I₁ / N₁ - (γ₂ + μ) * I₂]
    sys = complete(System(eqs, t; name = :sti_vertical))
    ngm = next_generation_matrix(sys, [J₁, I₁, J₂, I₂];
        transmission = [p * μ * I₁, β₁ * S₁ * I₂ / N₂, p * μ * I₂, β₂ * S₂ * I₁ / N₁])
    # T as in the paper: pμ into the juveniles, βᵢNᵢ/Nⱼ between the sexes, Σ diagonal
    @test symbolic_isequal(ngm.T[1, 2], p * μ)
    @test symbolic_isequal(ngm.T[3, 4], p * μ)
    @test symbolic_isequal(ngm.T[2, 4], β₁ * N₁ / N₂)
    @test symbolic_isequal(ngm.T[4, 2], β₂ * N₂ / N₁)
    # maturation of infected juveniles is a transition; the sexual-transmission terms
    # cancel exactly out of Σ (they are in T) and must be structurally zero
    @test symbolic_isequal(ngm.Σ[2, 1], ν₁)
    @test symbolic_isequal(ngm.Σ[4, 3], ν₂)
    @test ReproductiveNumbers._structural_zero(ngm.Σ[2, 4])
    @test ReproductiveNumbers._structural_zero(ngm.Σ[4, 2])
    @test length(states_at_infection(ngm)) == 4
    # det K = 0: only the two infectious states are sources, so K_S is 2×2
    @test symbolic_isequal(det(ngm.K), 0)
    K_S = small_domain_matrix(ngm)
    @test size(K_S) == (2, 2)
    R0 = basic_reproduction_number(ngm)
    for _ in 1:5
        vals = Dict(v => 0.1 + rand(rng) for v in [μ, ν₁, ν₂, γ₁, γ₂, β₁, β₂, N₁, N₂])
        vals[p] = rand(rng)
        num = evaluate(ngm, vals)
        @test nz(num.K) ≈ nz(ReproductiveNumbers.to_number.(substitute.(K_S, Ref(vals))))
        @test ReproductiveNumbers.to_number(substitute(R0, vals)) ≈
              basic_reproduction_number(num)
        @test validate_decomposition(num)
        @test sign(basic_reproduction_number(num) - 1) ==
              sign(maximum(real, eigvals(num.T + num.Σ)))
    end
end

@testset "Section 4.3: bovine viral diarrhoea with persistently infected animals" begin
    # E exposed, I transiently infectious, Z immune but pregnant with an infected foetus,
    # P persistently infected (born from Z, and from P at the reduced birth rate μ - a).
    @parameters β₁ β₂ ν γ μ p₁ p₂ α a b N
    @variables S(t) E(t) I(t) Z(t) P(t) R(t)
    horizontal = S * (β₁ * I + β₂ * P) / N
    eqs = [D(S) ~ μ * N - horizontal - μ * S,
        D(E) ~ horizontal - (ν + μ) * E,
        D(I) ~ ν * E - (γ + μ) * I,
        D(Z) ~ p₁ * γ * I - (α + μ) * Z,
        D(P) ~ p₂ * α * Z + (μ - a) * P - (μ + b) * P,
        D(R) ~ (1 - p₁) * γ * I + α * Z - μ * R]
    sys = complete(System(eqs, t; name = :bvd))
    ngm = next_generation_matrix(sys, [E, I, Z, P];
        transmission = [horizontal, 0, 0, p₂ * α * Z + (μ - a) * P])
    @test isequal(states_at_infection(ngm), [E, P])
    @test size(ngm.K) == (2, 2)
    @test symbolic_isequal(ngm.K[1, 1], β₁ * ν / ((ν + μ) * (γ + μ)))
    @test symbolic_isequal(ngm.K[1, 2], β₂ / (μ + b))
    @test symbolic_isequal(ngm.K[2, 1], ν * p₁ * γ * p₂ * α / ((ν + μ) * (γ + μ) * (α + μ)))
    @test symbolic_isequal(ngm.K[2, 2], (μ - a) / (μ + b))
    R0 = basic_reproduction_number(ngm)
    for _ in 1:5
        vals = Dict(v => 0.1 + rand(rng) for v in [β₁, β₂, ν, γ, μ, α, b, N])
        vals[p₁] = rand(rng)
        vals[p₂] = rand(rng)
        vals[a] = vals[μ] * rand(rng)
        num = evaluate(ngm, vals)
        @test ReproductiveNumbers.to_number(substitute(R0, vals)) ≈
              basic_reproduction_number(num)
        @test validate_decomposition(num)
        @test sign(basic_reproduction_number(num) - 1) ==
              sign(maximum(real, eigvals(num.T + num.Σ)))
    end
end

@testset "Staged progression (van den Driessche & Watmough 2002, section 4.3)" begin
    # n infectious stages in series with stage-specific infectivities; K stays 1×1.
    n = 4
    @parameters β[1:n] ν[1:n] d N
    @variables S(t) I(t)[1:n]
    βv, νv, Iv = collect(β), collect(ν), collect(I)
    λ = S * sum(βv[k] * Iv[k] for k in 1:n) / N
    eqs = [D(S) ~ -λ,
        D(Iv[1]) ~ λ - (ν[1] + d) * Iv[1]]
    for k in 2:n
        push!(eqs, D(Iv[k]) ~ νv[k - 1] * Iv[k - 1] - (νv[k] + d) * Iv[k])
    end
    sys = complete(System(eqs, t; name = :staged))
    ngm = next_generation_matrix(sys, Iv; equilibrium = Dict(S => N))
    @test isequal(states_at_infection(ngm), [Iv[1]])
    @test size(ngm.K) == (1, 1)
    @test size(ngm.K_L) == (n, n)
    # vdD&W equation (4.4): R0 = Σₖ βₖ Πⱼ<ₖ νⱼ / Πⱼ≤ₖ (νⱼ + d)
    paper = sum(βv[k] * prod(νv[j] for j in 1:(k - 1); init = 1) /
                prod(νv[j] + d for j in 1:k) for k in 1:n)
    @test symbolic_isequal(basic_reproduction_number(ngm), paper)
end
