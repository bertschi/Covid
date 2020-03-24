## Try SEIR model using Turing.jl

using Turing
using DifferentialEquations
using Random
Random.seed!(123);

struct SEIRModel
    Λ
    μ
    β
    a
    γ
end

function seirdynamics!(du, u, θ::SEIRModel, t)
    S, E, I, R = u
    N = S + E + I + R
    du[1] = θ.Λ - θ.μ * S - θ.β * I / N * S
    du[2] = θ.β * I / N * S - (θ.μ + θ.a) * E
    du[3] = θ.a * E - (θ.γ + θ.μ) * I
    du[4] = θ.γ * I - θ.μ * R
end

function NegBinom2(μ, ϕ)
    NegativeBinomial(ϕ, ϕ / (μ + ϕ))
end
    
@model seir(infected, N, ::Type{T}=Vector{Float64}) where {T} = begin
    ## Model parameters
    beta ~ truncated(Normal(0, 5), 0, Inf)
    a ~ truncated(Normal(0, 5), 0, Inf)
    gamma ~ truncated(Normal(0, 5), 0, Inf)
    phi ~ truncated(Normal(0, 10), 0, Inf)
    
    x₀ = T(undef, 4)
    for i = 1:4
        x₀[i] ~ truncated(TDist(3), 0, Inf)
    end

    ## Integrate ODE and relate to observations
    u₀ = T(undef, 4)
    u₀[1] = N - x₀[1]
    for i = 2:4
        u₀[i] = x₀[i]
    end
    θ = SEIRModel(0.0, 0.0, beta, a, gamma)
    tspan = (0.0, length(infected))
    ode = ODEProblem(seirdynamics!, u₀, tspan, θ)
    sol = solve(ode, RK4())
    infected_hat = map(u -> u[3], sol(1:length(infected)).u)

    try
        for i = 1:length(infected)
            infected[i] ~ NegBinom2(infected_hat[i], phi)
        end
    catch
        Inf
    end
end
