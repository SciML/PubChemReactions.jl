
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

# cids = unique(reduce(vcat, map(x -> x.cids, jrxns)))
# proteins = filter(!isnothing, unique(reduce(vcat, map(x -> get(x, :protacxns, nothing), jrxns))))

# this is the indices from the reactions from wikipedia's table of the main reactions in glycolysis
wiki_rxn_idxs = [4, 6, 12, 13, 14, 15, 18, 20, 21, 22]

# this checks that all the species in the reactions have cids 
@test all(aligned[wiki_rxn_idxs])

wiki_rxns = PubChemReactions.pathway_reaction_to_reaction.(jrxns[wiki_rxn_idxs])
@parameters t
@named rs = ReactionSystem(wiki_rxns, t)
sys = convert(ODESystem, rs)
prob = ODEProblem(sys, ones(length(states(sys))), (0, 100))
sol = solve(prob, Rosenbrock23())

# using Plots
# plot(sol)
