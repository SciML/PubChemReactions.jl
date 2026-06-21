using PubChemReactions, Catalyst, OrdinaryDiffEq
# sadly it seems only reactions that are balanced with all 1 stoich values are working with `get_pathway`

@variables t
# PathBank Citric Acid Cycle
pid = "PathBank:SMP0000057"
rxns = PubChemReactions.get_pathway(pid)
@named rs = ReactionSystem(rxns, Catalyst.DEFAULT_IV)
# Catalyst v16: systems must be `complete`d before building a problem, and an
# `ODEProblem` is constructed directly from the ReactionSystem (`convert(ODESystem, ...)`
# is deprecated). `states` was renamed to `unknowns`.
rs = complete(rs)
sts = unknowns(rs)
prob = ODEProblem(rs, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())

# # Reactome Glycolysis
pid = "Reactome:R-HSA-70171"
rxns2 = PubChemReactions.get_pathway(pid)
@named rs2 = ReactionSystem(rxns2, Catalyst.DEFAULT_IV)
rs2 = complete(rs2)
sts = unknowns(rs2)
prob = ODEProblem(rs2, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())
