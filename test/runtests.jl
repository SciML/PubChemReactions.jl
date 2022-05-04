using PubChemReactions, Catalyst
using Test

@parameters t

@testset "PubChemReactions.jl" begin
    @testset "graph" begin include("graph.jl") end
    @testset "species" begin include("species.jl") end
    @testset "pathway" begin include("pathway.jl") end
end