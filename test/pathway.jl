using PubChemReactions, Catalyst, OrdinaryDiffEq
# sadly it seems only reactions that are balanced with all 1 stoich values are working with `get_pathway`

@variables t
# PathBank Citric Acid Cycle
pid = "PathBank:SMP0000057"
rxns, failed = PubChemReactions.get_pathway(pid)
@named rs = ReactionSystem(rxns, Catalyst.DEFAULT_IV)
@test_broken length(failed) == 0
@test length(reactions(rs)) == 7
@test isbalanced(rs)

sys = convert(ODESystem, rs)
sts = states(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())


# # Reactome Glycolysis 
pid = "Reactome:R-HSA-70171"
rxns, failed = PubChemReactions.get_pathway(pid)
@named rs2 = ReactionSystem(rxns, Catalyst.DEFAULT_IV)
@test_broken length(failed) == 0
@test length(reactions(rs2)) == 8
@test isbalanced(rs2)

sys = convert(ODESystem, rs2)
sts = states(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())

f = failed[4]
s, p = PubChemReactions.subs_prods_from_pathway_rxn(f)
rxn = Reaction(1, s, p)
@info get_name.(species(rxn))
A = atom_matrix(rxn)
