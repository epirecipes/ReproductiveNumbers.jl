using ReproductiveNumbers, ModelingToolkit, Symbolics, LinearAlgebra, Random, Test
using OrdinaryDiffEqTsit5
include(joinpath(@__DIR__, "..", "models.jl"))
rng = MersenneTwister(4)

@testset "evaluate and numeric R0" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    p = Dict(m.β => 0.5, m.σ => 0.25, m.γ => 0.2, m.μ => 0.01, m.N => 1000.0)
    num = evaluate(ngm, p)
    @test num isa NextGenerationMatrix{Matrix{Float64}}
    @test evaluate(num, p) === num
    @test isequal(num.infected, ngm.infected)
    @test num.T ≈ [0 0.5; 0 0]
    @test num.Σ ≈ [-0.26 0; 0.25 -0.21]
    expected = 0.5 * 0.25 / (0.26 * 0.21)
    @test basic_reproduction_number(num) ≈ expected
    @test basic_reproduction_number(ngm, p) ≈ expected
    @test basic_reproduction_number(ngm, collect(p)) ≈ expected
    @test basic_reproduction_number(m.sys, [m.E, m.I], p) ≈ expected
    @test num.equilibrium[m.S] == 1000.0
    @test validate_decomposition(num)
    @test validate_decomposition(ngm, p; verbose = true)
    # missing parameter is reported
    @test_throws ArgumentError evaluate(ngm, Dict(m.β => 0.5))
end

@testset "parameter values from an ODEProblem" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    prob = ODEProblem(m.sys,
        [m.S => 999.0, m.E => 0.0, m.I => 1.0, m.R => 0.0,
            m.β => 0.5, m.σ => 0.25, m.γ => 0.2, m.μ => 0.01, m.N => 1000.0],
        (0.0, 100.0))
    @test basic_reproduction_number(ngm, prob) ≈ 0.5 * 0.25 / (0.26 * 0.21)
    # a state left symbolic in the equilibrium is read from the problem's initial state
    s = sir_model()
    ngm_s = next_generation_matrix(s.sys, [s.I]; equilibrium = Dict(s.R => 0))
    prob_s = ODEProblem(s.sys,
        [s.S => 500.0, s.I => 1.0, s.R => 0.0, s.β => 1.0, s.γ => 0.5, s.N => 1000.0],
        (0.0, 1.0))
    @test basic_reproduction_number(ngm_s, prob_s) ≈ 1.0 * 500 / (0.5 * 1000)
end

@testset "numeric agrees with symbolic for every model" begin
    for (m, inf, eq) in ((seir_demography_model(), (:E, :I), nothing),
        (two_latent_model(), (:E₁, :E₂, :I), nothing),
        (two_host_model(), (:E₁, :E₂, :I₁, :I₂), nothing))
        ngm = next_generation_matrix(m.sys, collect(inf))
        R0 = basic_reproduction_number(ngm)
        params = Num.(parameters(m.sys))
        for _ in 1:3
            vals = random_parameters(rng, params)
            haskey(vals, :p) && (vals[:p] = rand(rng))
            for k in keys(vals)
                string(k) == "p" && (vals[k] = rand(rng))
            end
            @test basic_reproduction_number(ngm, vals) ≈
                  ReproductiveNumbers.to_number(substitute(R0, vals))
        end
    end
end

@testset "type reproduction numbers" begin
    m = two_host_model()
    ngm = next_generation_matrix(m.sys, [m.E₁, m.E₂, m.I₁, m.I₂])
    vals = random_parameters(rng, Num.(parameters(m.sys)))
    K = evaluate(ngm, vals).K
    # single type formula: T_1 = K11 + K12 K21 / (1 - K22)
    if K[2, 2] < 1
        T1 = type_reproduction_number(ngm, 1, vals)
        @test T1 ≈ K[1, 1] + K[1, 2] * K[2, 1] / (1 - K[2, 2])
        @test type_reproduction_number(ngm, m.E₁, vals) ≈ T1
        @test type_reproduction_number(ngm, :E₁, vals) ≈ T1
        # symbolic version agrees
        T1s = type_reproduction_number(ngm, m.E₁)
        @test ReproductiveNumbers.to_number(substitute(T1s, vals)) ≈ T1
        # threshold agreement
        @test (T1 > 1) == (basic_reproduction_number(ngm, vals) > 1)
    end
    # the full set gives back R0
    @test type_reproduction_number(ngm, [1, 2], vals) ≈ basic_reproduction_number(ngm, vals)
    @test symbolic_isequal(
        type_reproduction_number(ngm, [m.E₁, m.E₂]), basic_reproduction_number(ngm))
    @test_throws ArgumentError type_reproduction_number(ngm, 3, vals)
    @test_throws ArgumentError type_reproduction_number(ngm, :Z, vals)
    @test_throws ArgumentError type_reproduction_number(ngm, m.I₁, vals)
    # undefined when the other types sustain transmission on their own
    big = copy(vals)
    big[m.β₂₂] = 1000.0
    @test_throws ArgumentError type_reproduction_number(ngm, 1, big)
end

@testset "NoClosedFormError falls back to numerics" begin
    @parameters β[1:3, 1:3] γ N[1:3]
    βm = collect(β)
    Nv = collect(N)
    @variables S(t)[1:3] I(t)[1:3]
    Sv = collect(S)
    Iv = collect(I)
    eqs = [[D(Sv[i]) ~ -Sv[i] * sum(βm[i, j] * Iv[j] / Nv[j] for j in 1:3) for i in 1:3];
           [D(Iv[i]) ~ Sv[i] * sum(βm[i, j] * Iv[j] / Nv[j] for j in 1:3) - γ * Iv[i]
            for i in 1:3]]
    sys = complete(System(eqs, t; name = :age))
    ngm = next_generation_matrix(sys, Iv; equilibrium = Dict(Sv[i] => Nv[i] for i in 1:3))
    @test size(ngm.K) == (3, 3)
    @test_throws NoClosedFormError basic_reproduction_number(ngm)
    B = [1.0 0.5 0.2; 0.5 1.5 0.4; 0.2 0.4 0.8]
    vals = Dict{Any, Any}(βm[i, j] => B[i, j] for i in 1:3, j in 1:3)
    for i in 1:3
        vals[Nv[i]] = 100.0
    end
    vals[γ] = 0.5
    @test basic_reproduction_number(ngm, vals) ≈ maximum(abs, eigvals(B / 0.5))
end
