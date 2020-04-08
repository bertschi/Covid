## Code showing identification issues

using Plots
using Distributions
using DifferentialEquations
using Optim
using DataFrames
using DataFramesMeta
using CSV

struct SIOUR
    βI::Function
    βO::Function
    βU::Function
    α::Float64
    γI::Float64
    γR::Float64
end

function siourdxdt!(du, u, θ::SIOUR, t)
    S, I, O, U, RO, RU = u
    N = S + I + O + U + RO + RU
    βI = θ.βI(t)
    βO = θ.βO(t)
    βU = θ.βU(t)
    du[1] = - (βI * I + βO * O + βU * U) * S / N 
    du[2] = (βI * I + βO * O + βU * U) * S / N - θ.γI * I
    du[3] = θ.α * θ.γI * I - θ.γR * O
    du[4] = (1 - θ.α) * θ.γI * I - θ.γR * U
    du[5] = θ.γR * O
    du[6] = θ.γR * U
end

function runmodel(dxdt, u₀, T, θ)
    tspan = (0.0, T)
    ode = ODEProblem(dxdt, u₀, tspan, θ)
    solve(ode, RK4())
end

R₀ = 3.0
τ = 10.0
θ = SIOUR(t -> R₀ / τ, t -> R₀ / τ, t -> R₀ / τ,
          0.1, 2 / τ, 2 / τ)
N = 1e8
u₀ = [N - 1, 1, 0, 0, 0, 0]
T = 200
sol = runmodel(siourdxdt!, u₀, T, θ)

plot(sol)
plot(sol.t, map(u -> u[3] + u[5], sol.u))

## Find best fitting sigmoid beta(t) for alpha = 1

function sigmoidal(t, x₁, x₂, β, τ)
    x₁ + (x₂ - x₁) / (1 + exp(- β * (t - τ)))
end

function NegBinom2(μ, ϕ)
    μ = max(zero(μ), μ)
    NegativeBinomial(ϕ, ϕ / (μ + ϕ))
end

function fitbeta(sol, θ)
    T = maximum(sol.t)
    u₀ = sol.u[1]
    Iobs = map(u -> u[3] + u[5], sol.u) ## cumulative observed cases

    params = log.([1, 0.1, 0.5, T / 2])
    function fit(params)
        β(t) = sigmoidal(t, map(exp, params)...)
        θfit = SIOUR(β, β, β, 1., θ.γI, θ.γR)
        runmodel(siourdxdt!, u₀, T, θfit)
    end
    function loss(params)
        model = fit(params)
        Ifit = map(u -> u[3] + u[5], model(sol.t).u)
        ## sum((log.(fit(params)[2:end]) .- log.(Iobs[2:end])).^2)
        sum((Ifit .- Iobs).^2)
    end
    opt = Optim.optimize(loss, params,
                         NelderMead(),
                         Optim.Options(iterations = 7500,
                                       show_trace = true,
                                       show_every = 100))
    opt, fit
end

opt, model = fitbeta(sol, θ)

df = vcat((@linq DataFrame(sol') |>
           transform(t = sol.t,
                     model = fill("True", length(sol)))),
          (@linq DataFrame(model(opt.minimizer)(sol.t)') |>
           transform(t = sol.t,
                     model = fill("Approx", length(sol)))),
          cols = :intersect)

df |> CSV.write("/tmp/foo.csv")

## library(tidyverse)
## library(ggthemes)

## df <- read_csv("/tmp/foo.csv")

## df %>%
##     ggplot(aes(t, x3 + x5,
##                color = model)) +
##     geom_line() +
##     theme_tufte() +
##     theme(legend.position = "top",
##           text = element_text(size = 12)) + 
##     scale_color_colorblind() +
##     labs(x = "Days",
##          y = TeX("Cumulated infections"),
##          color = "Model")

## ggsave("../../reports/figs/approx_infect.pdf")

## df %>%
##     ggplot(aes(t, x1,
##                color = model)) +
##     geom_line() +
##     theme_tufte() +
##     theme(legend.position = "top",
##           text = element_text(size = 12)) + 
##     scale_color_colorblind() +
##     labs(x = "Days",
##          y = TeX("Susceptible"),
##          color = "Model")

## ggsave("../../reports/figs/approx_suscept.pdf")
