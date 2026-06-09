using Pkg
using Test

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    using PubChemReactions, Catalyst, OrdinaryDiffEq

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
end

if GROUP == "QA"
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.instantiate()
    include(joinpath(@__DIR__, "qa", "qa.jl"))
end
