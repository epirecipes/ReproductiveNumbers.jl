using ReproductiveNumbers, ModelingToolkit, Symbolics, LinearAlgebra, Random, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using ReproductiveNumbers: manifestly_negative, isrankone, _independent_rows, tidy
include(joinpath(@__DIR__, "..", "models.jl"))
rng = MersenneTwister(11)

@testset "manifestly_negative" begin
    @variables a b c
    @test manifestly_negative(-a)
    @test manifestly_negative(-a / (b + c))
    @test manifestly_negative((-a - b) / c)
    @test !manifestly_negative((-a) / (-b))
    @test !manifestly_negative(a * b)
    @test manifestly_negative(-2)
    @test !manifestly_negative(0)
end

@testset "symbolic_iszero is deterministic" begin
    @variables a b
    x = sqrt(a * b) - sqrt(a) * sqrt(b)
    r = [ReproductiveNumbers.symbolic_iszero(x) for _ in 1:5]
    @test all(r)
end

@testset "suggest_infected" begin
    m = seir_demography_model()
    @test Set(string.(suggest_infected(m.sys))) == Set(["E(t)", "I(t)"])
    @test length(suggest_infected(m.sys; all = true)) == 1
    s = sir_model()
    @test Set(string.(suggest_infected(s.sys))) == Set(["I(t)"])
    v = ross_macdonald_model()
    @test Set(string.(suggest_infected(v.sys))) == Set(["I_H(t)", "I_V(t)"])
    e = sis_eliminated_model()
    @test Set(string.(suggest_infected(e.sys))) == Set(["I(t)"])
    two = two_strain_model()
    # the two strains are separate minimal invariant sets; the union is suggested
    @test Set(string.(suggest_infected(two.sys))) == Set(["I₁(t)", "I₂(t)"])
    @test Set(Set.(string.(c) for c in suggest_infected(two.sys; all = true))) ==
          Set([Set(["I₁(t)"]), Set(["I₂(t)"])])
    l = two_latent_model()
    @test Set(string.(suggest_infected(l.sys))) == Set(["E₁(t)", "E₂(t)", "I(t)"])
    vt = vertical_transmission_model()
    @test Set(string.(suggest_infected(vt.sys))) == Set(["J(t)", "I(t)"])
end

@testset "mean_sojourn_times" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    S = mean_sojourn_times(ngm)
    @test symbolic_isequal(S[1, 1], 1 / (m.σ + m.μ))
    @test symbolic_isequal(S[2, 1], m.σ / ((m.σ + m.μ) * (m.γ + m.μ)))
    @test symbolic_isequal(S[2, 2], 1 / (m.γ + m.μ))
    @test symbolic_isequal(S[1, 2], 0)
    num = evaluate(ngm, Dict(m.β => 0.5, m.σ => 0.25, m.γ => 0.2, m.μ => 0.01, m.N => 1.0))
    @test mean_sojourn_times(num) ≈ -inv(num.Σ)
end

@testset "sensitivities and elasticities" begin
    m = seir_demography_model()
    ngm = next_generation_matrix(m.sys, [m.E, m.I])
    s = sensitivities(ngm)
    R0 = m.β * m.σ / ((m.σ + m.μ) * (m.γ + m.μ))
    @test symbolic_isequal(s[m.β], m.σ / ((m.σ + m.μ) * (m.γ + m.μ)))
    @test symbolic_isequal(s[m.γ], -R0 / (m.γ + m.μ))
    e = elasticities(ngm; params = [m.β, m.γ, m.σ])
    @test symbolic_isequal(e[m.β], 1)
    @test symbolic_isequal(e[m.γ], -m.γ / (m.γ + m.μ))
    @test symbolic_isequal(e[m.σ], m.μ / (m.σ + m.μ))
    p = Dict(m.β => 0.5, m.σ => 0.25, m.γ => 0.2, m.μ => 0.01, m.N => 1.0)
    sn = sensitivities(ngm, p)
    @test sn[m.β] ≈ 0.25 / (0.26 * 0.21)
    en = elasticities(ngm, p; params = [m.β])
    @test en[m.β] ≈ 1.0
    @test Set(keys(sensitivities(R0))) == Set([m.β, m.σ, m.γ, m.μ])
end

@testset "perron_vectors" begin
    # numeric
    K = [1.0 2.0; 3.0 4.0]
    u, v = perron_vectors(K)
    ρ = maximum(abs, eigvals(K))
    @test K * u ≈ ρ * u
    @test sum(u) ≈ 1
    @test K' * v ≈ ρ * v
    @test dot(v, u) ≈ 1
    @test perron_vectors(zeros(0, 0)) == (Float64[], Float64[])
    # symbolic 1x1
    @variables a b c d
    @test isequal(perron_vectors(Num[a;;]), (Num[1], Num[1]))
    # symbolic 2x2 agrees with numeric
    us, vs = perron_vectors(Num[a b; c d])
    vals = Dict(a => 1.0, b => 2.0, c => 3.0, d => 4.0)
    @test ReproductiveNumbers.to_number.(substitute.(us, Ref(vals))) ≈ u
    @test ReproductiveNumbers.to_number.(substitute.(vs, Ref(vals))) ≈ v
    # rank one: K = x yᵀ
    @variables x1 x2 y1 y2
    Kr = Num[x1*y1 x1*y2; x2*y1 x2*y2]
    ur, vr = perron_vectors(Kr)
    @test symbolic_isequal(ur[1], x1 / (x1 + x2))
    @test symbolic_isequal(dot(vr, ur), 1)
    @test symbolic_isequal(vr[2] / vr[1], y2 / y1)
    @test_throws NoClosedFormError perron_vectors(Num[a 0; 0 d])
    @variables e f g h k
    @test_throws NoClosedFormError perron_vectors(Num[a b c; d e f; g h k])
    # from a next-generation matrix
    m = two_host_model()
    ngm = next_generation_matrix(m.sys, [m.E₁, m.E₂, m.I₁, m.I₂])
    p = random_parameters(rng, Num.(parameters(m.sys)))
    un, vn = perron_vectors(evaluate(ngm, p))
    @test sum(un) ≈ 1
    @test evaluate(ngm, p).K * un ≈ basic_reproduction_number(ngm, p) * un
end

@testset "effective reproduction number R_t" begin
    s = sir_model()
    # symbolic R_t as a function of the state
    Rt = effective_reproduction_number(s.sys, [s.I])
    @test symbolic_isequal(Rt, s.β * s.S / (s.γ * s.N))
    ngm = next_generation_matrix(s.sys, [s.I]; equilibrium = Dict())
    @test symbolic_isequal(effective_reproduction_number(ngm), Rt)
    # at a given state: partial and complete
    @test symbolic_isequal(
        effective_reproduction_number(ngm, Dict(s.S => s.N / 2)), s.β / (2 * s.γ))
    @test effective_reproduction_number(
        ngm, [s.S => 500.0, s.β => 0.5, s.γ => 0.25, s.N => 1000.0]) ≈ 1.0
    # some states fixed at construction
    @test symbolic_isequal(
        effective_reproduction_number(s.sys, [s.I]; equilibrium = Dict(s.S => s.N)),
        s.β / s.γ)
    # along a solution
    prob = ODEProblem(s.sys,
        [s.S => 990.0, s.I => 10.0, s.R => 0.0, s.β => 0.5, s.γ => 0.25, s.N => 1000.0],
        (0.0, 30.0))
    sol = solve(prob, Tsit5(); saveat = 5.0)
    traj = effective_reproduction_number(ngm, sol)
    @test length(traj) == length(sol.t)
    @test traj[1] ≈ 2.0 * 0.99
    @test issorted(traj; rev = true)
    @test effective_reproduction_number(ngm, sol, 10.0) ≈ 2.0 * sol(10.0; idxs = s.S) / 1000
    @test effective_reproduction_number(ngm, sol, [0.0, 30.0]) ≈ [traj[1], traj[end]]
    # as an observed variable of the system
    sys2, Rtvar = add_effective_reproduction_number(s.sys, [s.I])
    @test string(Rtvar) == "Rt(t)"
    prob2 = ODEProblem(sys2,
        [s.S => 990.0, s.I => 10.0, s.R => 0.0, s.β => 0.5, s.γ => 0.25, s.N => 1000.0],
        (0.0, 30.0))
    sol2 = solve(prob2, Tsit5(); saveat = 5.0)
    @test sol2[Rtvar] ≈ traj
    @test sol2(10.0; idxs = Rtvar) ≈ effective_reproduction_number(ngm, sol, 10.0)
    sys3, Rvar3 = add_effective_reproduction_number(s.sys, [s.I]; name = :Reff)
    @test string(Rvar3) == "Reff(t)"
    # a model without closed form works numerically at a state but not symbolically
    @parameters β[1:3, 1:3] γ N[1:3]
    βm, Nv = collect(β), collect(N)
    @variables S(t)[1:3] I(t)[1:3]
    Sv, Iv = collect(S), collect(I)
    eqs = [[D(Sv[i]) ~ -Sv[i] * sum(βm[i, j] * Iv[j] / Nv[j] for j in 1:3) for i in 1:3];
           [D(Iv[i]) ~ Sv[i] * sum(βm[i, j] * Iv[j] / Nv[j] for j in 1:3) - γ * Iv[i]
            for i in 1:3]]
    age = complete(System(eqs, t; name = :age))
    @test_throws NoClosedFormError effective_reproduction_number(age, Iv)
    ngm_age = next_generation_matrix(age, Iv; equilibrium = Dict())
    B = [1.0 0.5 0.2; 0.5 1.5 0.4; 0.2 0.4 0.8]
    vals = Dict{Any, Any}(βm[i, j] => B[i, j] for i in 1:3, j in 1:3)
    for i in 1:3
        vals[Nv[i]] = 100.0
        vals[Sv[i]] = 50.0
    end
    vals[γ] = 0.5
    @test effective_reproduction_number(ngm_age, vals) ≈ maximum(abs, eigvals(B / 0.5)) / 2
    # SEIRS with reinfection: R_t counts reinfections of R
    @parameters βs σs κ γs μs ω Ns
    @variables Ss(t) Es(t) Is(t) Rs(t)
    λ = βs * Is / Ns
    seirs = complete(System(
        [D(Ss) ~ μs * Ns + ω * Rs - λ * Ss - μs * Ss,
            D(Es) ~ λ * Ss + σs * λ * Rs - (κ + μs) * Es,
            D(Is) ~ κ * Es - (γs + μs) * Is,
            D(Rs) ~ γs * Is - σs * λ * Rs - ω * Rs - μs * Rs],
        t;
        name = :seirs))
    Rt_s = effective_reproduction_number(seirs, [Es, Is])
    @test symbolic_isequal(Rt_s, βs * κ * (Ss + σs * Rs) / (Ns * (κ + μs) * (γs + μs)))
    # time-varying transmission: R₀ is refused, R_t is a function of t and S
    @parameters β₀ ε γt Nt
    @variables St(t) It(t) Rt_(t)
    βt = β₀ * (1 + ε * cos(2π * t / 365))
    seasonal = complete(System(
        [D(St) ~ -βt * St * It / Nt, D(It) ~ βt * St * It / Nt - γt * It, D(Rt_) ~ γt * It],
        t; name = :seasonal))
    @test_throws ArgumentError next_generation_matrix(
        seasonal, [It]; equilibrium = Dict(St => Nt))
    Rseason = effective_reproduction_number(seasonal, [It])
    @test symbolic_isequal(Rseason, βt * St / (γt * Nt))
    ngm_season = next_generation_matrix(
        seasonal, [It]; equilibrium = Dict(), autonomous = false)
    probs = ODEProblem(seasonal,
        [St => 990.0, It => 10.0, Rt_ => 0.0, β₀ => 0.3, ε => 0.5, γt => 0.2, Nt => 1000.0],
        (0.0, 400.0))
    sols = solve(probs, Tsit5(); saveat = 50.0)
    Rtraj = effective_reproduction_number(ngm_season, sols)
    expected = [0.3 * (1 + 0.5 * cos(2π * τ / 365)) * sols(τ; idxs = St) / (0.2 * 1000)
                for τ in sols.t]
    @test Rtraj ≈ expected
    sys_s, Rvs = add_effective_reproduction_number(seasonal, [It])
    probs2 = ODEProblem(sys_s,
        [St => 990.0, It => 10.0, Rt_ => 0.0, β₀ => 0.3, ε => 0.5, γt => 0.2, Nt => 1000.0],
        (0.0, 400.0))
    sols2 = solve(probs2, Tsit5(); saveat = 50.0)
    @test sols2[Rvs] ≈ expected
end

@testset "numeric construction from F and V (ForwardDiff)" begin
    # SEIR with demography, vdD&W style: x = (E, I), S fixed at N through p
    p = (β = 0.5, σ = 0.25, γ = 0.2, μ = 0.01, N = 1000.0)
    F(x, p) = [p.β * p.N * x[2] / p.N, 0.0]
    V(x, p) = [(p.σ + p.μ) * x[1], -p.σ * x[1] + (p.γ + p.μ) * x[2]]
    ngm = next_generation_matrix(F, V, zeros(2), p; infected = [:E, :I])
    @test ngm.method == :functions
    @test ngm.T ≈ [0 0.5; 0 0]
    @test ngm.Σ ≈ [-0.26 0; 0.25 -0.21]
    @test basic_reproduction_number(ngm) ≈ 0.5 * 0.25 / (0.26 * 0.21)
    @test isequal(states_at_infection(ngm), [Symbolics.variable(:E)])
    # violating the sign conventions warns
    Vbad(x, p) = [-(p.σ + p.μ) * x[1], (p.γ + p.μ) * x[2]]
    @test_logs (:warn, r"sign conventions") next_generation_matrix(F, Vbad, zeros(2), p)
    @test_logs next_generation_matrix(F, Vbad, zeros(2), p; check = false)
end

@testset "general small-domain matrix" begin
    @variables a b c d e
    # rank-two T with proportional rows
    T = Num[0 a 0 b; 0 0 0 0; 0 c 0 d; 0 0 0 0]
    @test ReproductiveNumbers._independent_rows(T) == [1, 3]
    Σ = Num[-1 0 0 0; 0 -2 0 0; 0 0 -3 0; 0 0 0 -4]
    ngm = next_generation_matrix(T, Σ)
    K_S = small_domain_matrix(ngm)
    @test size(K_S) == (2, 2)
    vals = Dict(a => 1.0, b => 2.0, c => 0.5, d => 1.5)
    num = evaluate(ngm, vals)
    @test sort(abs.(filter(x -> abs(x) > 1e-9, eigvals(num.K)))) ≈
          sort(abs.(filter(x -> abs(x) > 1e-9,
        eigvals(ReproductiveNumbers.to_number.(substitute.(K_S, Ref(vals)))))))
    # numeric small domain
    K_Sn = small_domain_matrix(num)
    @test size(K_Sn) == (2, 2)
    @test spectral_radius(K_Sn) ≈ basic_reproduction_number(num)
    # full rank: K itself is returned
    ngm2 = next_generation_matrix(Num[a b; c d], Num[-e 0; 0 -e])
    @test small_domain_matrix(ngm2) === ngm2.K
    # zero T
    ngm3 = next_generation_matrix(Num[0 0; 0 0], Num[-e 0; 0 -e])
    @test small_domain_matrix(ngm3) === ngm3.K
end

@testset "small domain unlocks a closed form: section 4.2-type model" begin
    # Sexually transmitted infection with vertical transmission: juveniles J₁, J₂ are born
    # infected from infected adults I₁, I₂, mature at ν, adults infect the other sex.
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
    @test length(states_at_infection(ngm)) == 4
    @test_throws NoClosedFormError spectral_radius(ngm.K)
    K_S = small_domain_matrix(ngm)
    @test size(K_S) == (2, 2)
    R0 = basic_reproduction_number(ngm)      # falls back to K_S
    vals = random_parameters(rng, [μ, ν₁, ν₂, γ₁, γ₂, β₁, β₂, N₁, N₂])
    vals[p] = 0.4
    @test ReproductiveNumbers.to_number(substitute(R0, vals)) ≈
          basic_reproduction_number(ngm, vals)
end

@testset "K and K_L share their non-zero spectrum" begin
    for (m, inf) in (
        (seir_demography_model(), (:E, :I)), (two_latent_model(), (:E₁, :E₂, :I)),
        (two_host_model(), (:E₁, :E₂, :I₁, :I₂)))
        ngm = next_generation_matrix(m.sys, collect(inf))
        num = evaluate(ngm, random_parameters(rng, Num.(parameters(m.sys))))
        nz(A) = sort(abs.(filter(x -> abs(x) > 1e-9, eigvals(A))))
        @test nz(num.K) ≈ nz(num.K_L)
    end
end

@testset "irreducible_blocks gives a block upper triangular permutation" begin
    for _ in 1:20
        n = 6
        A = (rand(rng, n, n) .< 0.25) .* rand(rng, n, n)
        blocks = irreducible_blocks(A)
        perm = vcat(blocks...)
        B = A[perm, perm]
        offsets = cumsum([0; length.(blocks)])
        # every non-zero entry lies in a diagonal block or above it
        for I in 1:n, J in 1:n
            B[I, J] == 0 && continue
            bi = findlast(o -> o < I, offsets)
            bj = findlast(o -> o < J, offsets)
            @test bi <= bj
        end
    end
end

@testset "2x2 threshold criterion (mirrors TwoByTwo.lean)" begin
    for _ in 1:50
        a, b, c, d = rand(rng, 4) .* 3
        ρ = spectral_radius([a b; c d])
        @test (ρ > 1) == (a + d > 2 || 1 - (a + d) + (a * d - b * c) < 0)
    end
end

@testset "symbolic inverse scales to a dozen compartments" begin
    # inv(::Matrix{Num}) uses Laplace expansion and would never finish here
    n = 6
    @parameters q σ γ C[1:n, 1:n] N[1:n]
    @variables S(t)[1:n] E(t)[1:n] I(t)[1:n]
    Cm, Nv, Sv, Ev, Iv = collect(C), collect(N), collect(S), collect(E), collect(I)
    λ = [q * sum(Cm[i, j] * Iv[j] / Nv[j] for j in 1:n) for i in 1:n]
    eqs = vcat([D(Sv[i]) ~ -λ[i] * Sv[i] for i in 1:n],
        [D(Ev[i]) ~ λ[i] * Sv[i] - σ * Ev[i] for i in 1:n],
        [D(Iv[i]) ~ σ * Ev[i] - γ * Iv[i] for i in 1:n])
    age = complete(System(eqs, t; name = :age12))
    stats = @timed next_generation_matrix(
        age, vcat(Ev, Iv); equilibrium = Dict(Sv[i] => Nv[i] for i in 1:n))
    ngm = stats.value
    @test stats.time < 120
    @test size(ngm.K_L) == (12, 12) && size(ngm.K) == (6, 6)
    @test symbolic_isequal(ngm.K[1, 2], q * Cm[1, 2] * Nv[1] / (Nv[2] * γ))
    vals = Dict{Any, Any}(Cm[i, j] => 1.0 + i + j for i in 1:n, j in 1:n)
    for i in 1:n
        vals[Nv[i]] = 100.0 + i
    end
    vals[q] = 0.05
    vals[σ] = 0.25
    vals[γ] = 0.2
    num = evaluate(ngm, vals)
    @test num.K_L ≈ -num.T / num.Σ
    @test mean_sojourn_times(num) ≈ -inv(num.Σ)
    Sm = mean_sojourn_times(ngm)
    @test symbolic_isequal(Sm[n + 1, 1], σ / (σ * γ))
    # explicit inverse helper
    @variables a b c d
    Ainv = ReproductiveNumbers.symbolic_inverse(Num[a b; c d])
    @test symbolic_isequal(Ainv[1, 1], d / (a * d - b * c))
    @test symbolic_isequal(ReproductiveNumbers.symbolic_inverse(Num[a;;])[1, 1], 1 / a)
end
