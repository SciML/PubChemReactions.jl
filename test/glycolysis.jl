# no enzymes

@variables t

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
rxns2, failed2 = PubChemReactions.get_pathway(pid)
cids = unique(reduce(vcat, map(x -> x.cids, jrxns)))
proteins = filter(!isnothing, unique(reduce(vcat, map(x -> get(x, :protacxns, nothing), jrxns))))

rxn_str = jr.reaction

function pc_pathway_rxn_to_rp_cids(jr)
    rcids = get(jr, :cidsreactant, [])
    pcids = get(jr, :cidsproduct, [])
    to_arr.((rcids, pcids))
end

to_arr(xs) = isa(xs, AbstractArray) ? xs : [xs]

"check that the reaction str `jr.reaction`, when parsed, has the same length as the rcids and pcids from the json"
function is_reacts_prods_cids_aligned(jr)
    rxn_str = jr.reaction
    rp = pc_pathway_rxn_to_rp_cids(jr)
    sp = rhea_to_reacts_prods(rxn_str)
    length.(rp) == length.(sp)
end
aligned = is_reacts_prods_cids_aligned.(jrxns)

# I am ignoring the hexokinase reactions here
"GCK + GCKR ⟶ GCK1:GKRP complex"
"GCK1:GKRP complex ⟶ GCK + GCKR"
"GCK1:GKRP complex ⟶ GCK + GCKR"


 "1 , GCK + GCKR ⟶ GCK1:GKRP complex"
 "2 , GCK1:GKRP complex ⟶ GCK + GCKR"
 "3 , GCK1:GKRP complex ⟶ GCK + GCKR"
 "4 , ATP + Glc ⟶ ADP + G6P + H+"
#  "5 , ADP + Glc ⟶ AMP + G6P"
 "6 , G6P ⟶ Fru(6)P"
 "7 , GlcN6P + H2O ⟶ Fru(6)P + NH4+"
 "8 , H2O + pPF2K-Pase complex ⟶ PFKFB1 dimer + Pi"
"9 , ATP + PFKFB1 dimer ⟶ ADP + pPF2K-Pase complex"
 "10 , ATP + Fru(6)P ⟶ ADP + D-Fructose 2,6-bisphosphate"
 "11 , D-Fructose 2,6-bisphosphate + H2O ⟶ Fru(6)P + Pi"
 "12 , ATP + Fru(6)P ⟶ ADP + F1,6PP"
 "13 , F1,6PP ⟶ DHAP + GA3P"
 "14 , DHAP ⟶ GA3P"
 "15 , GA3P + NAD + Pi ⟶ 1,3BPG + H+ + NADH"
 "16 , 1,3BPG ⟶ 2,3BPG"
 "17 , 1,3BPG + G6P ⟶ 3PG + G1,6BP + H+"
 "18 , 1,3BPG + ADP ⟶ 3PG + ATP"
 "19 , 3PG + H2O ⟶ Glycerol + Pi"
 "20 , 3PG ⟶ 2PG"
 "21 , 2PG ⟶ H2O + PEP"
 "22 , ADP + H+ + PEP ⟶ ATP + PYR"

# 1	Glucose + ATP4− → Glucose-6-phosphate2− + ADP3− + H+	−16.7	−34 "4" "ATP + Glc ⟶ ADP + G6P + H+"

# 2	Glucose-6-phosphate2− → Fructose-6-phosphate2−	1.67	−2.9  "6" 

# 3	Fructose-6-phosphate2− + ATP4− → Fructose-1,6-bisphosphate4− + ADP3− + H+	−14.2	−19 "12" 


# 4	Fructose-1,6-bisphosphate4− → Dihydroxyacetone phosphate2− + Glyceraldehyde-3-phosphate2−	23.9	−0.23 "13"

# 5	Dihydroxyacetone phosphate2− → Glyceraldehyde-3-phosphate2−	7.56	2.4  "14"

# 6	Glyceraldehyde-3-phosphate2− + Pi2− + NAD+ → 1,3-Bisphosphoglycerate4− + NADH + H+	6.30	−1.29 "15"
# 7	1,3-Bisphosphoglycerate4− + ADP3− → 3-Phosphoglycerate3− + ATP4−	−18.9	0.09 "18"

# 8	3-Phosphoglycerate3− → 2-Phosphoglycerate3−	4.4	0.83 "20"
# 9	2-Phosphoglycerate3− → Phosphoenolpyruvate3− + H2O	1.8	1.1 "21"
# 10	Phosphoenolpyruvate3− + ADP3− + H+ → Pyruvate− + ATP4−	−31.7	−23.0 "22"

wiki_rxn_idxs = [4, 6, 12, 13, 14, 15, 18, 20, 21, 22]
jrxns[wiki_rxn_idxs]
aligned[wiki_rxn_idxs]

function jr_to_reaction(jr)
    rcids, pcids = pc_pathway_rxn_to_rp_cids(jr)

    rs = species_from_cid.(rcids)
    ps = species_from_cid.(pcids)
    @info rs, ps
    Reaction(1, rs, ps)
end

r = jr_to_reaction(jr)
rxns = jr_to_reaction.(jrxns[wiki_rxn_idxs])
@parameters t
@named rs = ReactionSystem(rxns, t)
