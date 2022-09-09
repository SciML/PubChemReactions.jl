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

# this is the indices from the reactions from wikipedia's table of the main reactions in glycolysis
wiki_rxn_idxs = [4, 6, 12, 13, 14, 15, 18, 20, 21, 22]

# this checks that all the species in the reactions have cids 
@test all(aligned[wiki_rxn_idxs])

wiki_rxns = PubChemReactions.pathway_reaction_to_reaction.(jrxns[wiki_rxn_idxs])
@parameters t
@named rs = ReactionSystem(wiki_rxns, t)
rxns = wiki_rxns

# not sure why reaction 3 is not balanced 
@test_broken isbalanced(rs)
good = rxns[Not(3)]
@test all(isbalanced.(good))

r = rxns[3]
Hydron = states(rs)[3]
new_r = Reaction(r.rate, r.substrates, [r.products..., Hydron])
rxns[1:2, 4:end]
new_rxns = [rxns[1:2]..., new_r, rxns[4:end]...]

@named new_rs = ReactionSystem(new_rxns, t)
@test isbalanced(new_rs)

sys = convert(ODESystem, new_rs)
prob = ODEProblem(sys, ones(length(states(sys))), (0, 100))
sol = solve(prob, Rosenbrock23())

sts = states(new_rs)

# for st in sts
#     ex = :((; $st) = $new_rs)
#     @info ex
# end

(; Adenosinetriphosphate, var"alpha-D-Glucopyranose", Hydron, var"alpha-D-glucose 6-phosphate(2-)", var"Adenosine-diphosphate", var"beta-D-fructofuranose 6-phosphate(2-)", var"Fructose 1,6-bisphosphate", var"D-glyceraldehyde 3-phosphate(2-)", var"Glycerone phosphate(2-)", var"Diphosphopyridine nucleotide", var"Hydrogen phosphate", var"NADH dianion", var"3-phosphonato-D-glyceroyl phosphate(4-)", var"3-phosphonato-D-glycerate(3-)", var"2-phosphonato-D-glycerate(3-)", Phosphonatoenolpyruvate, Water, Pyruvate) = new_rs


@unpack Adenosinetriphosphate, var"alpha-D-Glucopyranose", Hydron, var"alpha-D-glucose 6-phosphate(2-)", var"Adenosine-diphosphate", var"beta-D-fructofuranose 6-phosphate(2-)", var"Fructose 1,6-bisphosphate", var"D-glyceraldehyde 3-phosphate(2-)", var"Glycerone phosphate(2-)", var"Diphosphopyridine nucleotide", var"Hydrogen phosphate", var"NADH dianion", var"3-phosphonato-D-glyceroyl phosphate(4-)", var"3-phosphonato-D-glycerate(3-)", var"2-phosphonato-D-glycerate(3-)", Phosphonatoenolpyruvate, Water, Pyruvate = new_rs

defaults = [
    # these are from wikipedia https://en.wikipedia.org/wiki/Glycolysis#Free_energy_changes
    Adenosinetriphosphate => 1.85u"mM",
    var"Adenosine-diphosphate" => 0.14u"mM",
    var"alpha-D-Glucopyranose" => 5.0u"mM",
    var"alpha-D-glucose 6-phosphate(2-)" => 0.083u"mM",
    var"beta-D-fructofuranose 6-phosphate(2-)" => 0.014u"mM",
    var"Fructose 1,6-bisphosphate" => 0.031u"mM",
    var"D-glyceraldehyde 3-phosphate(2-)" => 0.019u"mM",
    var"Glycerone phosphate(2-)" => 0.14u"mM",
    var"Hydrogen phosphate" => 1.0u"mM",
    var"3-phosphonato-D-glyceroyl phosphate(4-)" => 0.001u"mM",
    var"3-phosphonato-D-glycerate(3-)" => 0.12u"mM",
    var"2-phosphonato-D-glycerate(3-)" => 0.03u"mM",
    Phosphonatoenolpyruvate => 0.023u"mM",
    Pyruvate => 0.051u"mM",

    # im guessing these, definitely wrong
    Hydron => 1.0u"mM",
    var"Diphosphopyridine nucleotide" => 1.0u"mM",
    var"NADH dianion" => 1.0u"mM",
    Water => 1.0u"mM"
]
@named def_rs = ReactionSystem(new_rxns, t; defaults)
@test ModelingToolkit.validate(def_rs)

new_sys = convert(ODESystem, def_rs)
new_prob = ODEProblem(new_sys, defaults, (0, 100))
@test_throws Unitful.DimensionError solve(new_prob, Rosenbrock23())
@test_throws Unitful.DimensionError solve(new_prob, Tsit5())

defaults = first.(defaults) .=> ustrip.(last.(defaults))
@named def_rs = ReactionSystem(new_rxns, t; defaults)

new_sys = convert(ODESystem, def_rs)
new_prob = ODEProblem(new_sys, defaults, (0, 100))
sol = solve(new_prob, Rosenbrock23())
plot(sol)
atp_ = sol[Adenosinetriphosphate]
@test_broken atp_[end]

# next steps are adding bidirectionality to applicable reactions, reaction rates, and enzymes/hill rates

# cids = unique(reduce(vcat, map(x -> x.cids, jrxns)))
# proteins = filter(!isnothing, unique(reduce(vcat, map(x -> get(x, :protacxns, nothing), jrxns))))
