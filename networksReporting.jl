#=
Implements the model described in Fogli and Veldkamp, "Germs, Social Networks, and Growth."
@author: Arnav Sood
@date: 2017-09-23
=#


include("networksBackend.jl")
using networksBackend, LightGraphs

function avgTechGrowth!(s::ModelState) # Log diffs last period's max tech to A₀.
    T = s.m.T
    A₀ = s.m.A₀
    maxTech = s.runtimeData.maxTech;
    avgTechGrowth = log((maxTech-A₀))/T
    return avgTechGrowth
end

function fracDisease!(s::ModelState)
    N = s.m.N
    numSick = count(isSick, s.nodeList)
    fracDisease = numSick/N
    return fracDisease
end

function avgPathLength!(s::ModelState)
    g = deepcopy(s.g) # Create a new copy of the state graph.
    g = Graph(g) # Convert to undirected to count paths.
    distances = floyd_warshall_shortest_paths(g).dists # Use the given algorithm to compute the distance matrix.
    avgPathLength = mean(distances) # Average over all pairs, including (i, i). Returns something like Matlab's Inf when any i, j have no path between them.
    return avgPathLength
end

function typeChangeProb!(s::ModelState)
    N = s.m.N
    T = s.m.T
    numChanges = s.runtimeData.typeChanges
    return numChanges / (N * T)
end

function avgLag!(s::ModelState)
    frontiers = s.runtimeData.frontierVec;
    lags = diff(frontiers)
    avgLag = mean(lags)
    return avgLag
end

function report(s::ModelState)
    avgTechGrowth = avgTechGrowth!(s)
    fracDisease = fracDisease!(s)
    avgPathLength = avgPathLength!(s)
    typeChangeProb = typeChangeProb!(s)
    avgLag = avgLag!(s)
    # halfDiffusion = halfDiffusion!(s)

    return avgTechGrowth, fracDisease, avgPathLength, typeChangeProb, avgLag
end
