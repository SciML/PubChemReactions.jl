using PubChemReactions, Catalyst, Test

atp = PubChemReactions.symbolic_species_from_name("atp");

foo = PubChemReactions.parse_rhea_equation2("H2O + lactose = D-galactose + D-glucose")
rxn = Reaction(1, foo...; only_use_rate=true)
@named rn = ReactionSystem([rxn], Catalyst.DEFAULT_IV, PubChemReactions.rxnspecies(rxn), [])
@test PubChemReactions.isbalanced(rn)

foo = PubChemReactions.parse_rhea_equation2("alpha,alpha-trehalose + H2O = 2 D-glucose")
rxn = Reaction(1, foo...; only_use_rate=true)
@named rn = ReactionSystem([rxn], Catalyst.DEFAULT_IV, PubChemReactions.rxnspecies(rxn), [])
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
