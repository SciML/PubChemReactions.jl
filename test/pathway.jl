using PubChemReactions, Catalyst, OrdinaryDiffEq
# sadly it seems only reactions that are balanced with all 1 stoich values are working with `get_pathway`

@variables t
# PathBank Citric Acid Cycle
pid = "PathBank:SMP0000057"
rxns = PubChemReactions.get_pathway(pid)
@named rs = ReactionSystem(rxns, Catalyst.DEFAULT_IV)

sys = convert(ODESystem, rs)
sts = states(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())

# # Reactome Glycolysis 
pid = "Reactome:R-HSA-70171"
rxns2 = PubChemReactions.get_pathway(pid)
@named rs2 = ReactionSystem(rxns2, Catalyst.DEFAULT_IV)

sys = convert(ODESystem, rs2)
sts = states(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())
