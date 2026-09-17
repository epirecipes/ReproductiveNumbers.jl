# R0 > 1 iff the infection-free steady state is unstable (spectral bound of T + Σ > 0),
# and iff a simulated epidemic initially grows.
using ReproductiveNumbers, ModelingToolkit, Symbolics, LinearAlgebra, Random, Test
using OrdinaryDiffEqTsit5
include(joinpath(@__DIR__, "..", "models.jl"))
rng = MersenneTwister(5)

spectral_bound(J) = maximum(real, eigvals(J))

@testset "sign(R0 - 1) == sign(s(T + Σ))" begin
    models = [(seir_demography_model(), (:E, :I)),
        (two_latent_model(), (:E₁, :E₂, :I)),
        (two_host_model(), (:E₁, :E₂, :I₁, :I₂)),
        (vertical_transmission_model(), (:J, :I))]
    for (m, inf) in models
        ngm = next_generation_matrix(m.sys, collect(inf))
        params = Num.(parameters(m.sys))
        for _ in 1:20
            vals = random_parameters(rng, params; lo = 0.05, hi = 3.0)
            for k in keys(vals)
                string(k) == "p" && (vals[k] = rand(rng))
            end
            num = evaluate(ngm, vals)
            @test validate_decomposition(num)
            R0 = basic_reproduction_number(num)
            s = spectral_bound(num.T + num.Σ)
            @test sign(R0 - 1) == sign(s) || isapprox(R0, 1; atol = 1e-8)
        end
    end
end

@testset "R0 = 1 exactly when the Jacobian is singular" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    # choose β so that R0 = 1
    vals = Dict(m.σ => 0.3, m.γ => 0.2, m.μ => 0.02, m.N => 1.0)
    vals[m.β] = (0.3 + 0.02) * (0.2 + 0.02) / 0.3
    num = evaluate(ngm, vals)
    @test basic_reproduction_number(num) ≈ 1
    @test abs(det(num.T + num.Σ)) < 1e-12
    @test abs(det(num.K_L - I)) < 1e-12
end

@testset "simulated SIR epidemics grow iff R0 > 1" begin
    m = sir_model()
    ngm = next_generation_matrix(m.sys, [m.I]; equilibrium = Dict(m.S => m.N, m.R => 0))
    for (β, γ) in ((0.5, 0.25), (0.2, 0.25), (0.25, 0.25))
        p = [m.β => β, m.γ => γ, m.N => 1000.0]
        prob = ODEProblem(m.sys, [m.S => 999.0, m.I => 1.0, m.R => 0.0, p...], (0.0, 5.0))
        sol = solve(prob, Tsit5(); saveat = 1.0)
        R0 = basic_reproduction_number(ngm, prob)
        @test R0 ≈ β / γ
        I0, I5 = sol[m.I][1], sol[m.I][end]
        if R0 > 1.05
            @test I5 > I0
        elseif R0 < 0.95
            @test I5 < I0
        end
    end
end
