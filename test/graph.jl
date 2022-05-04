using PubChemReactions, Catalyst, Test
using OrdinaryDiffEq
using Graphs

C6H12O6 = PubChemReactions.search_compound("glucose")
H2O = PubChemReactions.search_compound("water")

df = PubChemReactions.get_biochem_rxns(C6H12O6, H2O)
eqs = df[!, :Equation]
eq = first(eqs) # I happen to know this one has no prefixed stoich data

h2o = PubChemReactions.search_compound("water")
@test PubChemReactions.isspecies(h2o)
@test PubChemReactions.get_graph(h2o) isa PubChemReactions.AtomBondGraph

@species αGlc(t) [cid = 79025, save = true, load = true]
@species Glc(t) [cid = 5793, save = true, load = true]
g1 = get_graph(αGlc).g
g2 = get_graph(Glc).g
@test Graphs.Experimental.has_isomorph(g1, g2)
