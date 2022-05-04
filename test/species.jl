using PubChemReactions, Catalyst, Test
using Graphs, Symbolics
using Unitful # maybe testdep

@variables t
@species CO(t)
co_cpd = getmetadata(CO, Compound);
@test co_cpd.name == "Cobalt"

# "CO" didn't resolve right
@species CO(t) [cid = 281]
co_cpd = getmetadata(CO, Compound);
@test co_cpd.name == "Carbon monoxide"

# test manually passing in PubChem data
@variables CO(t)
CO = PubChemReactions.tospecies(CO; jsons=PubChemReactions.get_json_and_view_from_cid(281))
@test getmetadata(CO, Compound).name == "Carbon monoxide"

@species N2_gas(t) [name = "N2"]
n2 = getmetadata(N2_gas, Compound);
@test n2.name == "Nitrogen"
ag = getmetadata(N2_gas, PubChemReactions.AtomBondGraph)
@test nv(ag.g) == 2 # test it's N2 and not elemental Nitrogen
@test PubChemReactions.get_molecular_formula(N2_gas) == "N2"

@species CO2(t) = 5 [unit = u"mol/L"]
@test Symbolics.getdefaultval(CO2) == 5

# reaction with charged species
@species Hplus(t) [name="H+"] OHminus(t) [name="OH-"] H2O(t)
@test PubChemReactions.get_charge.([Hplus, OHminus, H2O]) == [1, -1, 0]
@test isbalanced(balance(Reaction(1, [Hplus, OHminus], [H2O])))

# save/load 
@species Water(t)
cid = get_cid(Water)
p = joinpath(PubChemReactions.COMPOUNDS_DIR, string(cid))
isdir(p) && rm(p;recursive=true)
@test !isdir(p)
@species Water(t) [save=true]
@test isdir(p)
@species Water(t) [load=true]
