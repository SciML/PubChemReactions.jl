using PubChemReactions, Catalyst
using Test

@testset "PubChemReactions.jl" begin
    @testset "rhea" begin include("rhea.jl") end
end