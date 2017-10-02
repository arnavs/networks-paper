#=
Implements the model described in Fogli and Veldkamp, "Germs, Social Networks, and Growth."
@author: Arnav Sood
@date: 2017-09-27

Runs the necessary simulations.
=#

# Import parameters.
include("networksBackend.jl")
include("networksReporting.jl")
using Distributions, networksBackend, LightGraphs, Gadfly, GraphPlot

# Set parameters and create holder objects.
ξ = Bernoulli(0.2)
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

# Run model 100 times, storing output in memory.
results = Array{ModelState}(100)
@time map!(entry -> runModel(m), results, results) # Can parallelize this eventually.

# Inspect an arbitrary element.
element = results[50]
g = element.g

print(element.nodeList)
adjacency_matrix(g)
gplot(g, layout=circular_layout)

# Run our specialized reporting functions.
rep = report(element)
print(rep)
