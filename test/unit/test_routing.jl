using ReproductiveNumbers, Symbolics, LinearAlgebra, Logging, Test
using ReproductiveNumbers: symbolic_inverse, expression_size, first_success,
                           TIDY_SIZE_LIMIT,
                           LAPLACE_INVERSE_LIMIT, tidy
include(joinpath(@__DIR__, "..", "models.jl"))

@testset "first_success" begin
    @test first_success([:bad => () -> error("x"), :good => () -> 42]; what = "demo") == 42
    @test_throws ErrorException first_success([:bad => () -> error("x")])
    @test first_success([:bad => () -> error("x")]; onfail = errs -> length(errs)) == 1
    # failures are logged at debug level
    @test_logs (:debug, r"failed") (:debug, r"succeeded") min_level=Logging.Debug first_success([
        :bad => () -> error("x"), :good => () -> 1])
end

@testset "expression_size" begin
    @variables a b c
    @test expression_size(a) == 1
    @test expression_size(2) == 1
    @test expression_size(a * b + c) > expression_size(a * b)
    @test expression_size(Num[a b; c a]) == 4
end

@testset "symbolic_inverse routes by size" begin
    @variables a b c d
    A = Num[a b; c d]
    for method in (:laplace, :lu, :auto)
        Ainv = symbolic_inverse(A; method)
        @test symbolic_isequal(Ainv[1, 1], d / (a * d - b * c))
        @test symbolic_isequal(Ainv[2, 1], -c / (a * d - b * c))
    end
    @test_throws ArgumentError symbolic_inverse(A; method = :bogus)
    # a 10x10 bidiagonal matrix: Laplace would be hopeless, LU is used
    n = 10
    @variables x[1:n]
    xv = collect(x)
    L = Num[i == j ? -xv[i] : (i == j + 1 ? xv[j] : 0) for i in 1:n, j in 1:n]
    stats = @timed symbolic_inverse(L)
    @test stats.time < 60
    Linv = stats.value
    @test symbolic_isequal(Linv[1, 1], -1 / xv[1])
    @test symbolic_isequal(Linv[2, 1], -1 / xv[2])
    # lowering the Laplace limit forces LU for the 2x2 as well
    old = LAPLACE_INVERSE_LIMIT[]
    LAPLACE_INVERSE_LIMIT[] = 0
    @test symbolic_isequal(symbolic_inverse(A)[1, 2], -b / (a * d - b * c))
    LAPLACE_INVERSE_LIMIT[] = old
end

@testset "tidy skips fraction simplification above the size limit" begin
    @variables a b c
    x = (a * b) / ((-b - c) * (-a))
    @test isequal(tidy(x), b / (b + c))
    old = TIDY_SIZE_LIMIT[]
    TIDY_SIZE_LIMIT[] = 1
    @test symbolic_isequal(tidy(x), b / (b + c))     # value unchanged
    @test_logs (:debug, r"skipping simplify_fractions") min_level=Logging.Debug match_mode=:any tidy(x)
    TIDY_SIZE_LIMIT[] = old
end

@testset "basic_reproduction_number routes" begin
    m = two_latent_model()
    ngm = next_generation_matrix(m.sys, [m.E₁, m.E₂, m.I])
    R0 = basic_reproduction_number(ngm)
    @test symbolic_isequal(basic_reproduction_number(ngm; method = :blocks), R0)
    @test symbolic_isequal(basic_reproduction_number(ngm; method = :small_domain), R0)
    @test_throws ArgumentError basic_reproduction_number(ngm; method = :bogus)
    h = two_host_model()
    n2 = next_generation_matrix(h.sys, [h.E₁, h.E₂, h.I₁, h.I₂])
    @test_throws NoClosedFormError basic_reproduction_number(n2; method = :small_domain)
    @test symbolic_isequal(
        basic_reproduction_number(n2; method = :blocks), basic_reproduction_number(n2))
    @test_logs (:debug, r"closed-form R₀") min_level=Logging.Debug match_mode=:any basic_reproduction_number(ngm)
end
