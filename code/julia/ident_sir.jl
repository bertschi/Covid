## Model based on SIR dynamics

using DifferentialEquations
using StatsFuns
using Plots
plotly()

struct SIR
    ## Allow time varying infectivity
    β::Function
    γ::Float64
end

## SIR model
function sirdxdt!(du, u, θ::SIR, t)
    S, I, R = u
    N = S + I + R
    β = θ.β(t, u)
    du[1] = - β * S / N * I
    du[2] = β * S / N * I - θ.γ * I
    du[3] = θ.γ * I
end

function modelsim(I₀, θ::SIR, t, N)
    u₀ = [N - I₀, I₀, 0]
    ode = ODEProblem(sirdxdt!, big.(u₀), big.((minimum(t), maximum(t))), θ;
                     saveat = t, tstops = t,
                     abstol = big(10.0)^(-16), reltol = big(10.0)^(-16))
    solve(ode, RK4())
end

## Show that model is not identifiable
## First: Constant beta and fraction observed equals dynamic beta

cases(sol) = map(u -> u[2] + u[3], sol.u)
deaths(sol, cfr) = map(u -> cfr * u[3], sol.u)

I₀ = 1e-7
N = 1
γ = big(1 / 8)
## β₁(t, u) = γ * (3 + (2 - 3) * logistic(t - 70))
β₁(t, u) = γ * 3
t = 0:0.1:150
sol₁ = modelsim(I₀, SIR(β₁, γ), t, N)
plot(sol₁)

ρ = big(0.1)
β₂(t, u) = β₁(t, nothing) * ρ * sol₁(t)[1] / u[1] * sol₁(t)[2] / u[2]
sol₂ = modelsim(I₀, SIR(β₂, γ), t, N)
plot(sol₂)

cfr = big(0.01)
τ = 4.5
plot(t .+ τ,
     log10.(hcat(ρ .* cases(sol₁(t .+ τ)),
                 cases(sol₂(t)),
                 ρ .* deaths(sol₁(t .+ τ), cfr),
                 deaths(sol₂(t), cfr))))

plot(t, hcat(β₁.(t, nothing), β₂.(t, sol₂(t).u)))

## Erlang SIR model

struct ErlangSIR
    ## Allow time varying infectivity
    β::Function
    γ::Float64
    k::Int64
end

## SIR model
function erlangsirdxdt!(du, u, θ::ErlangSIR, t)
    @assert length(u) == θ.k + 2
    S = u[1]
    I = u[1 .+ (1:θ.k)]
    R = u[end]
    N = sum(u)
    β = θ.β(t, u)
    du[1] = - β * S / N * sum(I)
    du[2] = β * S / N * sum(I) - θ.k * θ.γ * I[1]
    for i = 2:θ.k
        du[i + 1] = θ.k * θ.γ * I[i - 1] - θ.k * θ.γ * I[i]
    end
    du[end] = θ.k * θ.γ * I[θ.k]
end

function modelsim(I₀, θ::ErlangSIR, t, N)
    u₀ = vcat(N - I₀, I₀, fill(0, θ.k - 1), 0)
    ode = ODEProblem(erlangsirdxdt!, u₀, (minimum(t), maximum(t)), θ;
                     saveat = t, tstops = t)
    solve(ode, RK4())
end

sole = modelsim(I₀, ErlangSIR(β₁, γ, 25), t, N)
plot(sole)

plot(t,
     map(u -> cfr * u[end], sole(t).u) ./
     map(u -> sum(u[2:end]), sole(t .- 1 / γ).u);
     ylim = (0, 3*cfr))
