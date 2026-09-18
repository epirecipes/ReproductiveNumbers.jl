using ReproductiveNumbers, Symbolics, Test
using ReproductiveNumbers: additive_terms, tidy, symbolic_iszero, symbolic_variables,
                           depends_on, to_number, substitution_map

@variables a b c d x y

@testset "additive_terms" begin
    ts = additive_terms(a * x - b * y + c)
    @test length(ts) == 3
    @test length(additive_terms(a * (x + y))) == 1        # products are not expanded
    @test length(additive_terms((a * x + b * y) / c)) == 2 # quotients of sums are split
    @test length(additive_terms(a * x)) == 1
    @test additive_terms(3)[1] isa Num
end

@testset "tidy" begin
    @test isequal(tidy((-a) / (-b - c)), a / (b + c))
    @test isequal(tidy((a * b) / ((-c - d) * (-a - b))), (a * b) / ((c + d) * (a + b)))
    @test tidy(2) == 2
    @test isequal(tidy(a), a)
    # cancelling sums collapse, distinct fractions keep their shape
    @test isequal(tidy(a / b - a / b), 0)
    @test isequal(tidy((a * b) / c + (-a * b) / c), 0)
    @test isequal(tidy(a / b + c / d), a / b + c / d)
    @test isequal(tidy(a / b + c / d - c / d), a / b)
    @test size(tidy([a -b; -c d])) == (2, 2)
end

@testset "symbolic_iszero / symbolic_isequal" begin
    @test symbolic_iszero(0)
    @test !symbolic_iszero(1)
    @test symbolic_iszero(a - a)
    @test symbolic_iszero(a * b / (a * b) - 1)
    @test symbolic_iszero((a + b)^2 - a^2 - 2a * b - b^2)
    @test !symbolic_iszero(a - b)
    @test !symbolic_iszero(a - b; numeric = false)
    @test symbolic_isequal((a + b) / c, a / c + b / c)
    @test !symbolic_isequal(a / c, b / c)
end

@testset "symbolic_variables / depends_on" begin
    vs = symbolic_variables(a * x + b)
    @test length(vs) == 3
    @test depends_on(a * x, [x])
    @test !depends_on(a * x, [y])
    @test !depends_on(a * x, Num[])
    @test length(symbolic_variables([a * x, b * x])) == 3
end

@testset "to_number / substitution_map" begin
    @test to_number(substitute(a / b, Dict(a => 1.0, b => 4.0))) == 0.25
    @test to_number(Num(2)) === 2.0
    @test_throws ArgumentError to_number(a / b)
    m = substitution_map([a => 1, b => 2])
    @test m[a] == 1 && m[b] == 2
    @test substitution_map(Dict(a => 1))[a] == 1
    @test substitution_map((a => 1, b => 2))[b] == 2
    @test_throws ArgumentError substitution_map(1.0, nothing)
end
