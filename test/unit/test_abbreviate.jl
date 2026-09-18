using ReproductiveNumbers, ModelingToolkit, Symbolics, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
include(joinpath(@__DIR__, "..", "models.jl"))

@testset "structural abbreviation of the SEIR model" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    ab = abbreviate(ngm)
    defs = ab.definitions
    names = Dict(string(k) => k for k in keys(defs))
    @test Set(keys(names)) == Set(["τ_E", "τ_I", "p_E_I"])
    @test symbolic_isequal(defs[names["τ_E"]], 1 / (m.σ + m.μ))
    @test symbolic_isequal(defs[names["τ_I"]], 1 / (m.γ + m.μ))
    @test symbolic_isequal(defs[names["p_E_I"]], m.σ / (m.σ + m.μ))
    τ_I, p = names["τ_I"], names["p_E_I"]
    @test isequal(ab.K[1, 1], m.β * τ_I * p) || isequal(ab.K[1, 1], m.β * p * τ_I)
    @test isequal(ab.Σ[1, 1], -1 / names["τ_E"])
    @test symbolic_isequal(basic_reproduction_number(ab), m.β * τ_I * p)
    # expanding recovers the original
    ex = expand_definitions(ab)
    @test isempty(ex.definitions)
    @test all(symbolic_isequal(ex.K[i, j], ngm.K[i, j]) for i in 1:1, j in 1:1)
    @test all(symbolic_isequal(ex.Σ[i, j], ngm.Σ[i, j]) for i in 1:2, j in 1:2)
    @test symbolic_isequal(expand_definitions(ab, basic_reproduction_number(ab)),
        basic_reproduction_number(ngm))
    # numbers: original parameters or abbreviations
    p_orig = Dict(m.β => 0.5, m.σ => 0.25, m.γ => 0.2, m.μ => 0.01, m.N => 1.0)
    @test basic_reproduction_number(ab, p_orig) ≈ basic_reproduction_number(ngm, p_orig)
    p_ab = Dict(m.β => 0.5, names["τ_E"] => 1 / 0.26,
        names["τ_I"] => 1 / 0.21, names["p_E_I"] => 0.25 / 0.26)
    @test basic_reproduction_number(ab, p_ab) ≈ basic_reproduction_number(ngm, p_orig)
    # printing lists the definitions
    str = sprint(show, MIME("text/plain"), ab)
    @test occursin("where", str) && occursin("τ_I = ", str)
    # idempotent: abbreviating again introduces nothing new
    @test length(abbreviate(ab).definitions) == 3
    # single-parameter rates are left alone
    s = sir_model()
    ab_sir = abbreviate(next_generation_matrix(
        s.sys, [s.I]; equilibrium = Dict(s.S => s.N, s.R => 0)))
    @test length(ab_sir.definitions) == 1        # τ_I = 1/γ
    @test symbolic_isequal(
        basic_reproduction_number(ab_sir), s.β * only(keys(ab_sir.definitions)))
end

@testset "structural abbreviation of a two-host model with named transmissions" begin
    m = two_host_model()
    ngm = next_generation_matrix(m.sys, [m.E₁, m.E₂, m.I₁, m.I₂])
    ab = abbreviate(ngm; transmissions = true)
    names = Dict(string(k) => k for k in keys(ab.definitions))
    @test haskey(names, "T_E₁_I₂") && haskey(names, "p_E₂_I₂") && haskey(names, "τ_I₁")
    @test isequal(ab.K[1, 2], names["T_E₁_I₂"] * names["p_E₂_I₂"] * names["τ_I₂"]) ||
          symbolic_isequal(expand_definitions(ab, ab.K[1, 2]), ngm.K[1, 2])
    @test all(symbolic_isequal(expand_definitions(ab, ab.K[i, j]), ngm.K[i, j])
    for i in 1:2, j in 1:2)
    vals = Dict(v => 0.3 + 0.1 * i for (i, v) in enumerate(Num.(parameters(m.sys))))
    @test basic_reproduction_number(ab, vals) ≈ basic_reproduction_number(ngm, vals)
    # the abbreviated closed form is much shorter
    @test ReproductiveNumbers.expression_size(basic_reproduction_number(ab)) <
          ReproductiveNumbers.expression_size(basic_reproduction_number(ngm))
end

@testset "user-defined abbreviations by elimination" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    @variables R_E p_E
    ab = abbreviate(ngm, [R_E => m.β / (m.γ + m.μ), p_E => m.σ / (m.σ + m.μ)])
    @test symbolic_isequal(basic_reproduction_number(ab), R_E * p_E)
    @test Set(string.(keys(ab.definitions))) == Set(["R_E", "p_E"])
    @test symbolic_isequal(expand_definitions(ab, basic_reproduction_number(ab)),
        basic_reproduction_number(ngm))
    p = Dict(m.β => 0.5, m.σ => 0.25, m.γ => 0.2, m.μ => 0.01, m.N => 1.0)
    @test basic_reproduction_number(ab, p) ≈ basic_reproduction_number(ngm, p)
    # choosing the parameter to eliminate
    ab2 = abbreviate(ngm, [R_E => m.β / (m.γ + m.μ)]; solve_for = [m.γ])
    @test symbolic_isequal(expand_definitions(ab2, basic_reproduction_number(ab2)),
        basic_reproduction_number(ngm))
    # combined with the structural form
    ab3 = abbreviate(abbreviate(ngm, [R_E => m.β / (m.γ + m.μ)]))
    @test symbolic_isequal(expand_definitions(ab3, basic_reproduction_number(ab3)),
        basic_reproduction_number(ngm))
    # errors
    @test_throws ArgumentError abbreviate(ngm, [m.β => m.σ / (m.σ + m.μ)])          # existing parameter
    @test_throws ArgumentError abbreviate(ngm, [R_E => m.β^2 / m.γ]; solve_for = [m.β])  # quadratic in β
    @test_throws ArgumentError abbreviate(ngm, [R_E => m.β / m.γ]; solve_for = [m.β, m.γ])
    @test_throws ArgumentError abbreviate(ngm, [R_E => 2])
    # numeric matrices pass through expand unchanged
    num = evaluate(ngm, p)
    @test expand_definitions(num) === num
end
