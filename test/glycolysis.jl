using PubChemReactions, Catalyst, OrdinaryDiffEq, Test
using InvertedIndices
using Unitful, Plots

# below was attempting to do this all manually
# ss = @species glucose(t) g6p(t) f6p(t) var"f1,6bp"(t) gadp(t) [cid = 729] dhap(t)
# ss2 = @species var"1,3bpg"(t) [cid = 683] var"3pg"(t) [cid = 439183] var"2pg"(t) [cid = 59] pep(t) [cid = 1005] pyr(t) [cid = 1060]
# ss3 = @species atp(t) hydron(t) adp(t) var"nad+"(t) [cid = 5893]
# S = vcat(ss, ss2, ss3)
# get_name.(S) .=> get_cid.(S)

# @parameters k[1:10]
# r1 = Reaction(k[1], [glucose, atp], [hydron, adp, g6p])
# r2 = Reaction((k[2], k[3]), [g6p], [f6p])
# r3 = Reaction(k[4], [f6p, atp], [var"f1,6bp", hydron, adp])
# r4 = Reaction(k[5], [var"f1,6bp"], [gadp, dhap])
# r5 = Reaction(k[6], [dhap], [gadp])

pid = "Reactome:R-HSA-70171"
jrxns = PubChemReactions.pathway_json(pid)
aligned = PubChemReactions.is_reacts_prods_cids_aligned.(jrxns)
all_rxns = PubChemReactions.get_pathway(pid)

# Positional indices into the PubChem pathway listing for the core glycolysis reactions
# shown in Wikipedia's table. These have shifted over time as upstream data is updated;
# they match the current 21-reaction listing for "Reactome:R-HSA-70171".
wiki_rxn_idxs = [4, 6, 12, 13, 14, 15, 18, 19, 20, 21]

# this checks that all the species in the reactions have cids
@test all(aligned[wiki_rxn_idxs])

wiki_rxns = PubChemReactions.pathway_reaction_to_reaction.(jrxns[wiki_rxn_idxs])
@parameters t
@named rs = ReactionSystem(wiki_rxns, t)
rs = complete(rs)
rxns = wiki_rxns

# Historically reaction 3 was unbalanced and worked around by injecting a Hydron; with
# the current upstream data the full system already balances, so no workaround is needed.
@test isbalanced(rs)
good = rxns[Not(3)]
@test all(isbalanced.(good))

# Catalyst v16: build the ODEProblem directly from the completed ReactionSystem
# (`convert(ODESystem, ...)` is deprecated) and use `unknowns` (was `states`).
sts = unknowns(rs)
prob = ODEProblem(rs, sts .=> 1.0, (0, 100.0))
sol = solve(prob, Rosenbrock23())
# Smoke test: a successful integration produces a non-trivial trajectory.
@test successful_retcode(sol)
@test length(sol.t) > 1

# Exercise the plot recipe (regression check); the result is not asserted.
plot(sol)

# next steps are adding bidirectionality to applicable reactions, reaction rates, and enzymes/hill rates

# cids = unique(reduce(vcat, map(x -> x.cids, jrxns)))
# proteins = filter(!isnothing, unique(reduce(vcat, map(x -> get(x, :protacxns, nothing), jrxns))))
