using PubChemReactions, Aqua, JET
using Test

@testset "Aqua" begin
    Aqua.test_all(PubChemReactions)
end

@testset "JET" begin
    JET.test_package(PubChemReactions; target_defined_modules = true)
end
