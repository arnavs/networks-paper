#=
Implements the model described in Fogli and Veldkamp, "Germs, Social Networks, and Growth."
@author: Arnav Sood
@date: 2017-09-27

Provides an example of how to use the backend.
=#


# Import dependencies.
include("networksBackend.jl")
using networksBackend
using Distributions
using LightGraphs
using GraphPlot, Gadfly

# Endogenous parameters.
ξ = Bernoulli(0.2) # Guess based on paper. .
T = 500 # Paper page 21.
N = 400 # Paper page 22.
F1 = true
F2 = true
λ = 0.008
δ = 80
A₀ = 1 # Paper page 20.
π = 0.12
ϕ = 0.35
θ = Bernoulli(0.14)
nf = [4, 6]
m = [0, 0.2]
Δ = 7
Λ = 0.0008
Π = 0.12

# Bundle the three node distributions into one multivariate distribution.
α = NodeDist(θ, nf, m)

# Create model.
m = Model(ξ, T, N, F1, F2, λ, δ, A₀, π, ϕ, α, Δ, Λ, Π)

# Grab result.
@time result = runModel(m)
g = result.g

# Inspect result.
print(result.nodeList)
adjacency_matrix(g)

# Visualize.
gplot(g, layout=circular_layout)
