using PubChemReactions, Catalyst, OrdinaryDiffEq
using Test

@parameters t

@testset "PubChemReactions.jl" begin
    @testset "graph" begin
        include("graph.jl")
    end
    @testset "species" begin
        include("species.jl")
    end
    # @testset "balance" begin include("balance.jl") end
    @testset "pathway" begin
        include("pathway.jl")
    end
    @testset "glycolysis" begin
        include("glycolysis.jl")
    end
end
