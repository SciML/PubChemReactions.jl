using PubChemReactions, Catalyst
using Test

@testset "PubChemReactions.jl" begin
    @testset "graph" begin include("graph.jl") end
end