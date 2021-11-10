using PubChemReactions, Catalyst, Test
using OrdinaryDiffEq
using Graphs

C6H12O6 = PubChemReactions.search_compound("glucose")
H2O = PubChemReactions.search_compound("water")

df = PubChemReactions.get_biochem_rxns(C6H12O6, H2O)
eqs = df[!,:Equation]
eq = first(eqs) # I happen to know this one has no prefixed stoich data


h2o = PubChemReactions.search_compound("water")
@test PubChemReactions.isspecies(h2o)
@test PubChemReactions.get_graph(h2o) isa PubChemReactions.AtomBondGraph

# redox (testing balancing with ionic compounds works)
cr2o7_2minus = PubChemReactions.search_compound("dichromate")
fe_2plus = PubChemReactions.search_compound("ferrous ion")
hydron = PubChemReactions.search_compound("H+") # test bondless species work okay
cr_3plus = PubChemReactions.search_compound("cr3+")
fe_3plus = PubChemReactions.search_compound("ferric ion")
h2o = PubChemReactions.search_compound("water")

sps = [cr2o7_2minus, fe_2plus, hydron, cr_3plus, fe_3plus, h2o]
@test PubChemReactions.get_charge.(sps) == [-2, 2, 1, 3, 3, 0]

subs = [cr2o7_2minus, fe_2plus, hydron]
prods = [cr_3plus, fe_3plus, h2o]

rxn = PubChemReactions.get_balanced_reaction(subs, prods)
@test PubChemReactions.isbalanced(rxn)

@parameters t k
# some reactions that work
rxn_strs = [
    # eq,
    # "alpha,alpha-trehalose + H2O = 2 D-glucose",
    # "ATP + D-Glucose = ADP + D-Glucose 6-phosphate",
    "ch4 + o2 = co2 + h2o",
    "H+ + L-Arginine + NADPH + O2 = NADP+ + L-Citrulline + Nitric oxide",
    # "K+ + H2O + ATP + Na+ = K+ + Phosphate + ADP + Na+"
]

for rxn_str in rxn_strs 
    @info rxn_str
    rxn = PubChemReactions.get_reaction(rxn_str)
    @named rs = ReactionSystem([rxn], t, PubChemReactions.species(rxn), [k]; checks=false)
    @test PubChemReactions.isbalanced(rs)
end

rn = @reaction_network begin
  c1, X --> 2X
  c2, X --> 0
  c3, 0 --> X
end c1 c2 c3;

@test_throws Exception PubChemReactions.isbalanced(rn) # missing atom counts 
