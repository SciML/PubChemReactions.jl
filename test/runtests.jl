using PubChemReactions, Catalyst
using Test

@testset "PubChemReactions.jl" begin
    @testset "graph" begin include("graph.jl") end
    @testset "species" begin include("species.jl") end
end