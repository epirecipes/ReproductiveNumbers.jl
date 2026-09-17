using ReproductiveNumbers, Symbolics, LinearAlgebra, Random, Test
using ReproductiveNumbers: isrankone, spectral_radius_2x2

@variables a b c d e f
rng = MersenneTwister(1)

@testset "irreducible_blocks" begin
    @test irreducible_blocks([1.0 0; 0 2.0]) == [[1], [2]] ||
          irreducible_blocks([1.0 0; 0 2.0]) == [[2], [1]]
    @test irreducible_blocks([0.0 1.0; 1.0 0.0]) == [[1, 2]]
    @test irreducible_blocks(Num[a b; c d]) == [[1, 2]]
    blocks = irreducible_blocks(Num[a b 0; c d 0; 0 e f])
    @test sort(blocks) == [[1, 2], [3]]
    @test irreducible_blocks(zeros(3, 3)) |> length == 3
    @test_throws DimensionMismatch irreducible_blocks(zeros(2, 3))
end

@testset "isrankone" begin
    @test isrankone(Num[a*c a*d; b*c b*d])
    @test !isrankone(Num[a b; c d])
    @test isrankone([1.0 2.0; 2.0 4.0])
    @test !isrankone([1.0 0.0; 0.0 1.0])
    @test isrankone(Num[a;;])
end

@testset "spectral_radius_2x2" begin
    @test isequal(spectral_radius_2x2(Num(0), b, c, Num(0)), sqrt(b * c))
    ρ = spectral_radius_2x2(a, b, c, d)
    for _ in 1:5
        vals = Dict(a => rand(rng), b => rand(rng), c => rand(rng), d => rand(rng))
        M = [vals[a] vals[b]; vals[c] vals[d]]
        @test ReproductiveNumbers.to_number(substitute(ρ, vals)) ≈ maximum(abs, eigvals(M))
    end
end

@testset "spectral_radius: symbolic closed forms" begin
    @test isequal(spectral_radius(Num[a;;]), a)
    @test isequal(spectral_radius(Num[a*c a*d; b*c b*d]), a * c + b * d)  # rank one -> trace
    @test isequal(spectral_radius(Num[0 b; c 0]), sqrt(b * c))
    @test isequal(spectral_radius(Num[a 0; 0 d]), max(a, d))
    @test isequal(spectral_radius(Num[a 0; 0 a]), a)
    @test isequal(spectral_radius(Num[a b; 0 d]), max(a, d))     # block triangular
    @test isequal(spectral_radius(Num[0 0; 0 0]), 0)
    @test isequal(spectral_radius(Num[a 0; 0 0]), a)
    @test isequal(spectral_radius(Matrix{Num}(undef, 0, 0)), 0)
    # generic 2x2 matches the numeric spectral radius
    ρ = spectral_radius(Num[a b; c d])
    for _ in 1:5
        vals = Dict(a => rand(rng), b => rand(rng), c => rand(rng), d => rand(rng))
        M = [vals[a] vals[b]; vals[c] vals[d]]
        @test ReproductiveNumbers.to_number(substitute(ρ, vals)) ≈ maximum(abs, eigvals(M))
    end
    # 3x3 with a 2x2 irreducible block and a decoupled state
    ρ3 = spectral_radius(Num[a b 0; c d 0; 0 0 e])
    @test isequal(ρ3, max(spectral_radius_2x2(a, b, c, d), e)) ||
          isequal(ρ3, max(e, spectral_radius_2x2(a, b, c, d)))
    # full 3x3 has no closed form
    @test_throws NoClosedFormError spectral_radius(Num[a b c; d e f; a b e])
end

@testset "spectral_radius: numeric" begin
    @test spectral_radius([2.0 0; 0 1]) == 2.0
    @test spectral_radius([0 4.0; 1 0]) ≈ 2.0
    @test spectral_radius(zeros(0, 0)) == 0.0
    A = rand(rng, 4, 4)
    @test spectral_radius(A) ≈ maximum(abs, eigvals(A))
end

@testset "characteristic_polynomial" begin
    @variables λ
    p = characteristic_polynomial(Num[a b; c d]; λ)
    @test symbolic_isequal(p, λ^2 - (a + d) * λ + (a * d - b * c))
end
