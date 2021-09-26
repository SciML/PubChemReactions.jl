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
