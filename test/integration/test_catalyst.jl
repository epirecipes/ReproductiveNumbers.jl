using ReproductiveNumbers, Catalyst, ModelingToolkit, Symbolics, Random, Test
rng = MersenneTwister(3)

@testset "SEIR reaction network" begin
    rn = @reaction_network seir begin
        @parameters β σ γ μ N
        β / N, S + I --> E + I
        σ, E --> I
        γ, I --> R
        μ * N, 0 --> S
        μ, (S, E, I, R) --> 0
    end
    @unpack S, E, I, R, β, σ, γ, μ, N = rn
    ngm = next_generation_matrix(rn, [E, I])
    @test isequal(states_at_infection(ngm), [E])
    @test symbolic_isequal(ngm.T[1, 2], β)
    @test symbolic_isequal(ngm.Σ[2, 1], σ)
    R0 = β * σ / ((σ + μ) * (γ + μ))
    @test symbolic_isequal(basic_reproduction_number(ngm), R0)
    @test symbolic_isequal(basic_reproduction_number(rn, [:E, :I]), R0)
    @test symbolic_isequal(disease_free_equilibrium(rn, [E, I])[S], N)
    sub = infected_subsystem(rn, [E, I])
    @test length(sub.f_infected) == 2
    # Catalyst and ModelingToolkit routes agree
    ngm_mtk = next_generation_matrix(Catalyst.ode_model(rn), [E, I])
    @test symbolic_isequal(basic_reproduction_number(ngm_mtk), R0)
end

@testset "SIR with S => N supplied" begin
    rn = @reaction_network sir begin
        @parameters β γ N
        β / N, S + I --> 2I
        γ, I --> R
    end
    @unpack S, I, R, β, γ, N = rn
    ngm = next_generation_matrix(rn, [I]; equilibrium = Dict(S => N, R => 0))
    @test symbolic_isequal(basic_reproduction_number(ngm), β / γ)
    @test_throws ArgumentError next_generation_matrix(rn, [I])
end

@testset "Vertical transmission is detected from the stoichiometry" begin
    rn = @reaction_network vertical begin
        @parameters β μ ν γ p N
        β / N, S + I --> 2I
        p * μ, I --> I + J          # infected birth: net +1 infected -> transmission
        (1 - p) * μ, I --> I + S
        μ, S --> 2S
        ν, J --> I
        γ + μ, I --> 0
        μ, (S, J) --> 0
    end
    @unpack S, J, I, β, μ, ν, γ, p, N = rn
    ngm = next_generation_matrix(rn, [J, I]; equilibrium = Dict(S => N))
    @test isequal(states_at_infection(ngm), [J, I])
    @test symbolic_isequal(ngm.T[1, 2], p * μ)
    @test symbolic_isequal(ngm.T[2, 2], β)
    @test symbolic_isequal(ngm.T[1, 1], 0) && symbolic_isequal(ngm.T[2, 1], 0)
    @test symbolic_isequal(ngm.Σ[2, 1], ν)
    @test symbolic_isequal(ngm.Σ[1, 1], -(ν + μ))
    # closed form: K is 2x2 with K = [pμν/((ν+μ)(γ+μ)), pμ/(γ+μ); βν/((ν+μ)(γ+μ)), β/(γ+μ)] which is rank one
    R0 = basic_reproduction_number(ngm)
    @test symbolic_isequal(R0, p * μ * ν / ((ν + μ) * (γ + μ)) + β / (γ + μ))
end

@testset "Overriding the transmission reactions" begin
    rn = @reaction_network malaria begin
        @parameters a b c γ μᵥ N_H N_V
        a * b / N_H, S_H + I_V --> I_H + I_V
        γ, I_H --> S_H
        a * c / N_H, S_V + I_H --> I_V + I_H
        μᵥ * N_V, 0 --> S_V
        μᵥ, (S_V, I_V) --> 0
    end
    @unpack S_H, I_H, S_V, I_V, a, b, c, γ, μᵥ, N_H, N_V = rn
    eq = Dict(S_H => N_H, S_V => N_V)
    ngm = next_generation_matrix(rn, [I_H, I_V]; equilibrium = eq)
    R0 = basic_reproduction_number(ngm)
    @test symbolic_isequal(R0^2, a^2 * b * c * N_V / (N_H * γ * μᵥ))
    # only vector -> host counts, by predicate
    ngm_pred = next_generation_matrix(rn, [I_H, I_V]; equilibrium = eq,
        transmission = rx -> any(isequal(I_H), rx.products) &&
            any(isequal(S_H), rx.substrates))
    @test symbolic_isequal(basic_reproduction_number(ngm_pred), R0^2)
    # ... by flags
    flags = [true, false, false, false, false, false]
    @test symbolic_isequal(
        basic_reproduction_number(next_generation_matrix(
            rn, [I_H, I_V]; equilibrium = eq, transmission = flags)),
        R0^2)
    # ... by index
    @test symbolic_isequal(
        basic_reproduction_number(next_generation_matrix(
            rn, [I_H, I_V]; equilibrium = eq, transmission = [1])),
        R0^2)
    @test_throws ArgumentError next_generation_matrix(
        rn, [I_H, I_V]; equilibrium = eq, transmission = [true])
    @test_throws ArgumentError next_generation_matrix(
        rn, [I_H, I_V]; equilibrium = eq, transmission = [99])
    @test_throws ArgumentError next_generation_matrix(
        rn, [I_H, I_V]; equilibrium = eq, transmission = "no")
end

@testset "Non-species infected states are rejected" begin
    rn = @reaction_network sir begin
        @parameters β γ N
        β / N, S + I --> 2I
        γ, I --> R
    end
    @unpack S, I, β = rn
    @test_throws ArgumentError next_generation_matrix(rn, [β]; equilibrium = Dict(S => 1))
end
