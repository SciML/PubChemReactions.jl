using Groebner
using Symbolics
using PubChemReactions
using Catalyst

@parameters t
@variables x[1:6]
rxns = []

@species CH4(t) [cid = 297, save = true, load = true]
@species O2(t) [cid = 977, save = true, load = true]
@species CO2(t) [cid = 280, save = true, load = true]
@species H2O(t) [cid = 962, save = true, load = true]

@species Dichromate(t) [cid = 24503, save = true, load = true]
@species FerrousIon(t) [cid = 27284, save = true, load = true]
@species Hydron(t) [cid = 1038, save = true, load = true]

@species ChromiumIII(t) [cid = 27668, save = true, load = true]
@species FerricIon(t) [cid = 29936, save = true, load = true]

@species P2I4(t) [cid = 83484, save = true, load = true]
@species P4(t) [cid = 123286, save = true, load = true]
@species PH4I(t) [cid = 166618, save = true, load = true]
@species H3PO4(t) [cid = 1004, save = true, load = true]
@species Hydroxide(t) [cid = 961, save = true, load = true]

rn = @reaction_network begin
    c1, X --> 2X
    c2, X --> 0
    c3, 0 --> X
end c1 c2 c3;

@test_throws Exception isbalanced(rn) # not species

rxn = Reaction(1, [CH4, O2], [CO2, H2O])
comb_rxn = rxn
brxn = balance(rxn)
push!(rxns, rxn)
@test !isbalanced(rxn)
@test isbalanced(brxn)
@test brxn.substoich == [1, 2]
@test brxn.prodstoich == [1, 2]

# ion test
rxn = Reaction(1, [Dichromate, FerrousIon, Hydron], [ChromiumIII, FerricIon, H2O])
push!(rxns, rxn)
brxn = balance(rxn)
@test isbalanced(brxn)

# this example tests the fact that just dividing by the minimum doesn't guarantee integer solutions
rxn = Reaction(1, [P2I4, P4, H2O], [PH4I, H3PO4])
push!(rxns, rxn)
brxn = balance(rxn)
@test isbalanced(brxn)

# this one was honestly a bitch
rxn = Reaction(1, [Hydroxide, Hydron], [H2O])
h2o_rxn = rxn
push!(rxns, rxn)
brxn = balance(rxn)
@test isbalanced(brxn)

# this reaction actually contains two reactions, more work will be required to separate out.
# TODO try this one in Mathematica
# 2 MnO4– + 5 H2O2 + 6 H+ → 2 Mn2+ + 5 O2 + 8 H2O
# 2 H2O2 → O2 + 2 H2O
s = @species var"Permanganate"(t) [cid = 24401, save = true, load = true] var"Hydrogen peroxide"(t) [
    cid = 784, save = true, load = true,
] var"Hydron"(t) [
    cid = 1038, save = true, load = true,
] var"Manganese(2+)"(t) [
    cid = 27845, save = true, load = true,
] var"Oxygen"(t) [
    cid = 977, save = true, load = true,
] var"Water"(t) [cid = 962, save = true, load = true]
rxn = Reaction(1, s[1:3], s[4:6])
push!(rxns, rxn)
@test_throws InexactError balance(rxn)

# rxn = PubChemReactions.get_reaction("H+ + L-Arginine + NADPH + O2 = NADP+ + L-Citrulline + Nitric oxide")
# @test !isbalanced(rxn)
# brxn = balance(rxn)
# @test isbalanced(brxn)

# @species var"Ozone"(t) [cid = 24823, save = true, load = true]
# @species var"Hydroxide"(t) [cid = 961, save = true, load = true]
# @species var"Hydrogenperoxide(1-)"(t) [cid = 18500, save = true, load = true]
# @species var"Oxygen"(t) [cid = 977, save = true, load = true]

# sps = reduce(vcat, [[O3, OH], [HO2, O2]])

# rxn = Reaction(1, [O3, OH], [HO2, O2])
# eqs = balance_eqs(rxn)
# weqs = eqs_to_mathematica(eqs)
