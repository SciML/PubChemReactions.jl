using PubChemReactions, Catalyst, Test
using OrdinaryDiffEq
using MetaGraphs, LightGraphs, NetworkLayout

atp = PubChemReactions.symbolic_species_from_name("atp");

foo = PubChemReactions.parse_rhea_equation("H2O + lactose = D-galactose + D-glucose")
rxn = Reaction(1, foo...; only_use_rate=true)
@test PubChemReactions.isbalanced(rxn)
@named rn = ReactionSystem([rxn], Catalyst.DEFAULT_IV, PubChemReactions.species(rxn), [])
@test PubChemReactions.isbalanced(rn)

foo = PubChemReactions.parse_rhea_equation("alpha,alpha-trehalose + H2O = 2 D-glucose")
rxn = Reaction(1, foo...; only_use_rate=true)
@named rn = ReactionSystem([rxn], Catalyst.DEFAULT_IV, PubChemReactions.species(rxn), [])
@test PubChemReactions.isbalanced(rn)

repressilator = @reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end α K n δ γ β μ;

@test_throws Exception PubChemReactions.isbalanced(repressilator) # missing counts 

h2o = species"water"
@test PubChemReactions.isspecies(h2o)
@test PubChemReactions.get_graph(h2o) isa AtomBondGraph


ch4 = species"ch4"
o2 = species"o2"
co2 = species"co2"
h2o = species"h2o"

subs = [ch4, o2]
prods = [co2, h2o]
gs, gs2 = map(x->get_graph(x).g, subs), map(x->get_graph(x).g, prods)
rxn = PubChemReactions.get_balanced_reaction(subs, prods)
@test PubChemReactions.isbalanced(rxn)

# redox 
cr2o7_2minus = species"dichromate"
fe_2plus = species"ferrous ion"
hydron = species"H+" # test bondless species work okay
cr_3plus = species"cr3+"
fe_3plus = species"ferric ion"
h2o = species"water"

sps = [cr2o7_2minus, fe_2plus, hydron, cr_3plus, fe_3plus, h2o]
@test PubChemReactions.get_charge.(sps) == [-2, 2, 1, 3, 3, 0]

subs = [cr2o7_2minus, fe_2plus, hydron]
prods = [cr_3plus, fe_3plus, h2o]

rxn = PubChemReactions.get_balanced_reaction(subs, prods)
@test PubChemReactions.isbalanced(rxn)

eq = "ATP + D-Glucose = ADP + D-Glucose 6-phosphate"
parsed_eq = PubChemReactions.parse_rhea_equation(eq)
rxn = reaction"ch4 + o2 = co2 + h2o"
@parameters t k=1
@named rs = ReactionSystem([rxn], t, PubChemReactions.species(rxn), [k]; checks=false)
@test PubChemReactions.isbalanced(rs)

sys = convert(ODESystem, rs)
# prob = ODEProblem(sys, [1,1,1,1], (0, 10.))
# sol = solve(prob, Tsit5())
# plot(sol)

# random reactions that I think work
rxn = reaction"H+ + L-Arginine + NADPH + O2 = NADP+ + L-Citrulline + Nitric oxide"
# rxn2 = reaction"Phosphate + Maltoheptaose <-> alpha-D-Glucose 1-phosphate + Maltohexaose"
rxn3 = reaction"K+ + H2O + ATP + Na+ = K+ + Phosphate + ADP + Na+ "


# test viz 
s = species"atp"
atomplot(s)
atomplot2(s)
atomplot3(s)

# this demonstrates the inconsistency for namings in PubChem
# you cant just use the @species_str for C3H8 and have it know its propane, since there are multiple possible configs
# alkanes = String[]
# for i in 2:5
#     push!(alkanes, "C$(i)H$(2i+2)")
# end

# sps = PubChemReactions.get_species.(alkanes)


# scratch of net energy and thermo stuff
# using Unitful
# rxn = reaction"ch4 + o2 = co2 + h2o"
# PubChemReactions.is_mass_conserved(rxn)
# c = 299_792_458u"m / s"
# m = u"g" * -(PubChemReactions.is_mass_conserved(rxn)...)
# E = m*c^2
# uconvert(Unitful.eV, E)

"below is calc for E=mc^2 for the energy released by D-T reaction"
# rxn = reaction"Deuterium + Tritium = Helium"
# D = species"Deuterium"
# atomplot2(D)
# T = species"Tritium"
# He = species"Helium"
# dm = PubChemReactions.get_mass(D)
# dt = PubChemReactions.get_mass(T)
# dhe = PubChemReactions.get_mass(He)


# Dmass = 2.01410177811u"u"
# Tmass = 3.01604928u"u"
# Hemass = 4.002603254u"u"
# nmass = 1.00866491588u"u"
# a = Dmass + Tmass
# b = Hemass + nmass
# @test uconvert(Unitful.MeV, (a - b)*c^2) ≈ 17.589299u"MeV"


# scratch pad for handling ionic bonds, pubchem does provide it if the compound has ionic bonds
# s = species"nacl"
# j = get_j(s)
# jv = get_jv(s)
# chg = j.atoms.charge


# s2 = species"O2"
# j2 = get_j(s2)
# jv = get_jv(s2)
# chg = j2.atoms.charge

# s2 = species"Iron(II) Oxide"
# j2 = get_j(s2)
# jv = get_jv(s2)
# chg = j2.atoms.charge
