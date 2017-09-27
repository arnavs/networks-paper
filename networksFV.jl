#=
Implements the model described in Fogli and Veldkamp, "Germs, Social Networks, and Growth."
@author: Arnav Sood
@date: 2017-09-27
=#

# Name module.
module networksFV

# Export types.
export Node, Model, DataHolder, ModelState, NodeDist

# Export functions.
export τ, isSick, Node, ourMod, runModel, update!, initialize, nodes

# Load dependencies. Alias some objects with our notation.
using Distributions
using LightGraphs
const Probability = Real
const nodes = vertices
import Base.rand

# Define the objects used in the model.
type Node
    θ::Bool # Collectivist flag.
    nf::Integer # Degree. Specifically, twice the invariant number of outgoing (left) bonds.
    m::Probability # Mobility probability.
    A::Number # Technology.
    mobilityFlag::Bool # Whether the agent has a mobility link or doesn't.
end

type Model
    ξ::Distribution # Shock distribution.
    T::Integer # Number of time periods.
    N::Integer # Number of agents.
    F1::Bool # Innovation feedback.
    F2::Bool # Infection feedback.
    λ::Real # Innovation Poisson parameter.
    δ::Real # Innovation gain.
    A₀::Real # Initial technology level.
    π::Probability # Disease transmission probability.
    ϕ::Probability # Technology transmission probability.
    α::Distribution # Distribution of initial node types.
    Δ::Integer # Spacing for individualist nodes.
    Λ::Real # Feedback 1 parameter.
    Π::Real # Feedback 2 parameter.
end

type DataHolder
    typeChanges::Integer
    maxTech::Real
    frontierVec::Array
end

# Create α.
type NodeDist <: Distribution{Multivariate, Discrete}
    θ::Distribution
    nf
    m
end

function rand(α::NodeDist, n::Int64=1)
    rtn = [rand(α.θ), rand(α.nf), rand(α.m)]
    return rtn
end

type ModelState
    m::Model
    g::DiGraph # Model state represented as a directed graph; left arrows are outgoing, and right arrows are incoming.
    nodeList::Array{Node} # Vector of nodes corresponding to each index.
    runtimeData::DataHolder # Captures data we need to track as the state evolves.
end

# Create accessors and constructors for these objects.
function τ(n::Node)
    return [n.θ, n.nf, n.m]
end

function isSick(n::Node)
    return iszero(n.A)
end

function Node(θ, nf, m)
    return Node(θ, nf, m, 0) # If unset, nodes initialize with A = 0.
end

function fillNode(m::Model)
    τ = rand(m.α) # Draw a node type.
    return Node(τ..., m.A₀, false) # ... unpacks f([x,y,z],γ) into f(x,y,z,γ).
end

function ModelState(m::Model, g::DiGraph, nodeList)
    return ModelState(m, g, nodeList, DataHolder(0, m.A₀, []))
end

# Modular arithmetic on [1...N] (instead of the default [0...N-1]).
function ourMod(n::Number, N::Integer)
    if n <= N && 0 < n
        return n
    elseif n >= N
        return ourMod(n - N, N)
    elseif n <= 0
        return ourMod(n + N, N)
    end
end

# Draw bonds based on type.
function setBonds(g::DiGraph, m::Model, indexSet::Array,
nodeList::Array{Node})
    N = m.N
    for index in indexSet
        nf = nodeList[index].nf
        if nodeList[index].θ == true
            # Collectivist.
            friends = collect(1:nf/2)
            for f in friends
                add_edge!(g, index, ourMod(index + f, N))
            end
        elseif nodeList[index].θ == false
            # Individualist.
            friends = [1]
            append!(friends, collect(m.Δ : m.Δ -2 + nf/2))
            for f in friends
                add_edge!(g, index, ourMod(index + f, N))
            end
        end
    end
end

# Initialize model.
function initialize(m::Model)
    # Create objects.
    nodeList = Array{Node}(m.N)
    g = DiGraph(m.N)
    # Fill the objects.
    map!(n -> fillNode(m), nodeList, nodeList)
    setBonds(g, m, collect(1:m.N), nodeList)
    # Create a ModelState around them.
    return ModelState(m, g, nodeList)
end

# Do technology transmission.
function updateTech!(s::ModelState, ϕ::Distribution)
    # Create a changes list, so that we run each node's changes independently.
    # deepcopy() ensures that we're actually creating an independent object.
    newList = deepcopy(s.nodeList)
    for nodeIndex in nodes(s.g)
        friendsList = all_neighbors(s.g, nodeIndex)
        techList = [s.nodeList[f].A for f in friendsList]
        # Do the shock for each.
        map!(n -> rand(ϕ) == 1 ? n : 0, techList, techList)
        bestNeighbor = maximum(techList)
        # Make sure the new tech is actually better than the current one.
        if bestNeighbor >= s.nodeList[nodeIndex].A
            newList[nodeIndex].A = bestNeighbor
        end
    end
    # Replace.
    s.nodeList = newList
end

# Do disease transmission.
function updateDisease!(s::ModelState, π::Distribution)
    newList = deepcopy(s.nodeList)
    for nodeIndex in nodes(s.g)
        numSick = count(isSick, [s.nodeList[n] for n in all_neighbors(s.g, nodeIndex)])
        draw = rand(π, numSick)
        # If at least one shock is positive, we infect the current node.
        if 1 ∈ draw
            newList[nodeIndex].A = 0
        end
    end
    s.nodeList = newList
end

# Tech innovation.
function updateInnovation!(s::ModelState, λ::Real)
    shock = Poisson(λ)
    # No need to store changes in a copy, since these events are individual.
    for node in s.nodeList
        draw = rand(shock)
        newTech = (1 + draw * s.m.δ) * node.A
        node.A = newTech
    end
end

# Check technology adoption lag.
function lagCheck!(s::ModelState, t::Integer)
    oldMax = s.runtimeData.maxTech;
    newMax = maximum(n.A for n in s.nodeList);
    if newMax > oldMax
        print_with_color(32, "New tech max attained at time $(t)\n")
        s.runtimeData.maxTech = newMax
        push!(s.runtimeData.frontierVec, t)
    end
end

# Perform tech frontier-related reporting.
function checkFrontiers!(s::ModelState, t::Integer)
    lagCheck!(s::ModelState, t::Integer)
end

# Update types.
function updateTypes!(s::ModelState, ξ::Array)
    newList = deepcopy(s.nodeList)
    changedIndices = find(ξ)
    for nodeIndex in changedIndices
        friendsList = all_neighbors(s.g, nodeIndex)
        techList = [s.nodeList[f].A for f in friendsList]
        A, loc = findmax(techList)
        # Increment type-change flag as necessary.
        if τ(s.nodeList[loc]) != τ(s.nodeList[nodeIndex])
            s.runtimeData.typeChanges += 1
        end
        # Replace type.
        newList[nodeIndex] = Node(τ(s.nodeList[loc])..., A, false) # false indicates that nodes start without mobility links.
    end
    s.nodeList = newList
end

# Type-changed agents update bonds.
function updateBonds!(s::ModelState, ξ::Array)
    changedIndices = find(ξ)
    for nodeIndex in changedIndices
        # Break all outgoing bonds.
        neighbors = deepcopy(out_neighbors(s.g, nodeIndex))
        for neighborIndex in neighbors
            rem_edge!(s.g, nodeIndex, neighborIndex)
        end
    end
    # Redraw new outgoing bonds, based on the node's type.
    setBonds(s.g, s.m, changedIndices, s.nodeList)
end

# Redraw bonds for mobile agents.
function redrawBonds!(s::ModelState, mobList::Array)
    newGraph = deepcopy(s.g)
    for nodeIndex in mobList
        shock = Bernoulli(s.nodeList[nodeIndex].m)
        draw = rand(shock)
        if s.nodeList[nodeIndex].mobilityFlag == false && draw == 1 # No existing mobility links, and a valid shock.
            neighbors = deepcopy(out_neighbors(s.g, nodeIndex)) # Neighbors from old graph.
            neighborIndex = rand(neighbors) # Pick one.
            rem_edge!(newGraph, nodeIndex, neighborIndex) # Break the corresponding bond on the new graph.
            unconnected = setdiff(collect(1:s.m.N), all_neighbors(newGraph, nodeIndex))
            # Remove the actual node.
            unconnected = setdiff(unconnected, [nodeIndex])
            newNeighbor = rand(unconnected)
            add_edge!(newGraph, nodeIndex, newNeighbor) # Add the new edge on the new graph.
            println("Replaced edge $(neighborIndex) with edge $(newNeighbor) for $(nodeIndex)")
            s.nodeList[nodeIndex].mobilityFlag = true; # Update mobility flag.
        end
    end
    # Replace the graph.
    s.g = newGraph
end

# Mobile agents redraw bonds.
function updateMobile!(s::ModelState)
    mobilities = [n.m for n in s.nodeList]
    mobileAgents = find(mobilities)
    println("Mobile agents:")
    println(mobileAgents)
    redrawBonds!(s, mobileAgents)
end

function update!(s::ModelState, t::Integer)
    # Unpack parameters, draw ξ shock vector, etc.
    m = s.m
    g = s.g
    nodeList = s.nodeList
    π = Bernoulli(m.π)
    Π = m.Π
    ϕ = Bernoulli(m.ϕ)
    ξ = rand(m.ξ, m.N) # This is an Array, instead of a Distribution!
    λ = m.λ
    Λ = m.Λ
    N = m.N
    # Update parameters for any feedbacks.
    if m.F1 == true
        numSick = count(isSick, nodeList)
        avgSick = numSick/N
        λ = 1 - cdf(Normal(0, 1), avgSick/Λ)
    end
    if m.F2 == true
        techList = [n.A for n in nodeList]
        meanTech = mean(techList)
        ϕ = Bernoulli(1-cdf(Normal(0, 1), meanTech/Π))
    end
    updateTech!(s,ϕ)
    updateDisease!(s,π)
    updateInnovation!(s, λ)
    checkFrontiers!(s, t::Integer)
    updateTypes!(s, ξ)
    updateBonds!(s, ξ)
    updateMobile!(s)
end

# Run model.
function runModel(m::Model)
    state = initialize(m)
    for t in 1:m.T # So initialization is at time 0.
        update!(state, t)
    end
    return state
end

# End module.
end
