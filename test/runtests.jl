using PubChemReactions, Catalyst, OrdinaryDiffEq
using Test
using SafeTestsets

@testset "PubChemReactions.jl" begin
    @safetestset "graph" begin
        include("graph.jl")
    end
    @safetestset "species" begin
        include("species.jl")
    end
    # @safetestset "balance" begin include("balance.jl") end
    @safetestset "pathway" begin
        include("pathway.jl")
    end
    @safetestset "glycolysis" begin
        include("glycolysis.jl")
    end
end
