using PubChemReactions, Catalyst, OrdinaryDiffEq
# sadly it seems only reactions that are balanced with all 1 stoich values are working with `get_pathway`

@variables t
# PathBank Citric Acid Cycle
pid = "PathBank:SMP0000057"
rxns, failed = PubChemReactions.get_pathway(pid)
@named rs = ReactionSystem(rxns, Catalyst.DEFAULT_IV)
# @test_broken length(failed) == 0
# @test length(reactions(rs)) == 7
# @test isbalanced(rs)

sys = convert(ODESystem, rs)
sts = states(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())


# # Reactome Glycolysis 
pid = "Reactome:R-HSA-70171"
rxns2, failed2 = PubChemReactions.get_pathway(pid)
r2, f2 = rxns2, failed2
@named rs2 = ReactionSystem(rxns, Catalyst.DEFAULT_IV)
# @test_broken length(failed) == 0
# @test length(reactions(rs2)) == 8
# @test isbalanced(rs2)

sys = convert(ODESystem, rs2)
sts = states(sys)
prob = ODEProblem(sys, sts .=> 1.0, (0, 100.0))
sol = OrdinaryDiffEq.solve(prob, Tsit5())

f = failed[4]
s, p = PubChemReactions.parse_pathway_reaction(f)
rxn = Reaction(1, s, p)
@info get_name.(species(rxn))
# eqs = balance_eqs(rxn)
# {12, 3, 9, 4, 9, 11}
# isbalanced(Reaction(1, s, p, [12, 3, 9], [4, 9, 11]))

A = atom_matrix(rxn)
xs = PubChemReactions.parse_pathway_reaction.(failed)
frxns = map(x -> Reaction(1, x...), xs)
# eqs_to_mathematica.(balance_eqs.(frxns))

# xs = PubChemReactions.parse_pathway_reaction.(f2)