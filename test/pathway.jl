using PubChemReactions, Catalyst, OrdinaryDiffEq
# sadly it seems only reactions that are balanced with all 1 stoich values are working with `get_pathway`

@variables t
# PathBank Citric Acid Cycle
pid = "PathBank:SMP0000057"
rxns = PubChemReactions.get_pathway(pid)
@named rs = ReactionSystem(rxns, Catalyst.DEFAULT_IV)
# Catalyst v15+: ReactionSystems must be `complete`d before conversion to other systems,
# and the resulting ODESystem must also be completed before constructing an ODEProblem.
rs = complete(rs)

sys = complete(convert(ODESystem, rs))
sts = unknowns(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())

# # Reactome Glycolysis
pid = "Reactome:R-HSA-70171"
rxns2 = PubChemReactions.get_pathway(pid)
@named rs2 = ReactionSystem(rxns2, Catalyst.DEFAULT_IV)
rs2 = complete(rs2)

sys = complete(convert(ODESystem, rs2))
sts = unknowns(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())
