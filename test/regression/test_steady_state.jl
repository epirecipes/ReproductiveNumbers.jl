# R₀ is consistent with the long-run behaviour of the full non-linear model: sweeping a
# parameter, the steady state reached from a small introduction is the infection-free one
# exactly when R₀ ≤ 1, and an endemic equilibrium (found with NonlinearSolve) exists
# exactly when R₀ > 1. For the SEIR model with demography the endemic equilibrium has the
# classical property S* = N / R₀.
using ReproductiveNumbers, ModelingToolkit, Symbolics, Test
using SteadyStateDiffEq, OrdinaryDiffEqTsit5, NonlinearSolve
include(joinpath(@__DIR__, "..", "models.jl"))

@testset "SEIR with demography: steady state versus R₀ over a parameter sweep" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    base = Dict(m.σ => 0.2, m.γ => 0.1, m.μ => 0.02, m.N => 1000.0)
    for β in (0.03, 0.06, 0.12, 0.2, 0.4)
        p = merge(base, Dict(m.β => β))
        R0 = basic_reproduction_number(ngm, p)
        u0 = [m.S => 999.0, m.E => 0.0, m.I => 1.0, m.R => 0.0]
        # build from an ODEProblem: the direct System constructor of SteadyStateProblem
        # hits an internal error in ModelingToolkit 11.44
        ssprob = SteadyStateProblem(ODEProblem(m.sys, [u0; collect(p)], (0.0, 1.0)))
        ss = solve(ssprob, DynamicSS(Tsit5()); abstol = 1e-10, reltol = 1e-10)
        Istar = ss[m.I]
        Sstar = ss[m.S]
        # near the threshold the approach to the steady state is slow, hence the tolerances
        if R0 < 1
            @test Istar < 1e-3
            @test Sstar ≈ 1000.0 atol = 0.1
        else
            @test Istar > 1e-2
            @test Sstar ≈ 1000.0 / R0 rtol = 1e-3
        end
    end
end

@testset "SEIR with demography: endemic equilibrium from NonlinearSolve" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    base = Dict(m.σ => 0.2, m.γ => 0.1, m.μ => 0.02, m.N => 1000.0)
    for β in (0.05, 0.15, 0.3)
        p = merge(base, Dict(m.β => β))
        R0 = basic_reproduction_number(ngm, p)
        nlprob = NonlinearProblem(m.sys, [[m.S => 500.0, m.E => 50.0, m.I => 50.0, m.R => 400.0]; collect(p)])
        sol = solve(nlprob, NewtonRaphson(); abstol = 1e-10)
        if R0 > 1
            @test SciMLBase.successful_retcode(sol)
            @test sol[m.I] > 0
            @test sol[m.S] ≈ 1000.0 / R0 rtol = 1e-6
        else
            # the only non-negative equilibrium is the infection-free one
            @test abs(sol[m.I]) < 1e-6 || sol[m.I] < 0
        end
    end
end

@testset "Ross-Macdonald: endemic prevalence appears exactly when R₀ > 1" begin
    m = ross_macdonald_model()
    eq = Dict(m.S_H => m.N_H, m.S_V => m.N_V)
    ngm = next_generation_matrix(m.sys, [m.I_H, m.I_V]; equilibrium = eq)
    base = Dict(m.b => 0.5, m.c => 0.5, m.γ => 0.1, m.μᵥ => 0.1, m.N_H => 1000.0, m.N_V => 2000.0)
    for a in (0.05, 0.1, 0.2, 0.4)
        p = merge(base, Dict(m.a => a))
        R0 = basic_reproduction_number(ngm, p)
        u0 = [m.S_H => 999.0, m.I_H => 1.0, m.S_V => 2000.0, m.I_V => 0.0]
        ss = solve(SteadyStateProblem(ODEProblem(m.sys, [u0; collect(p)], (0.0, 1.0))),
                   DynamicSS(Tsit5()); abstol = 1e-10, reltol = 1e-10)
        @test (ss[m.I_H] > 1e-3) == (R0 > 1)
    end
end
