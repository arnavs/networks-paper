# Import dependencies.
include("networksFV.jl")
using networksFV
using Distributions
using LightGraphs
using GraphPlot, Gadfly

# In this test code, the α initial states distribution is a Benoulli throwaway. Until I get a handle on the multivariate dist. code, I edited the backend to pick with probability 0.5 between individualist 0 or 1,  and with uniform probability between mobilities 0.2, 0.7 and 0, and have initial degree = 4.

# Test, see how long it takes.
m = Model(Bernoulli(0.5), 20, 20, true, true, 0.4, 0.2, 3, 0.3, 0.3, Bernoulli(0.5), 3, 0.15, 200);

@time result = runModel(m)
g = result.g

# Inspect results.
print(result.nodeList)
adjacency_matrix(g)

# Visualize.
gplot(g, layout=circular_layout)
