C6H12O6 = PubChemReactions.gen_sym("glucose")
H2O = PubChemReactions.gen_sym("water")

res = PubChemReactions.species_info(C6H12O6)
res[:InformationList][:Information]

df = PubChemReactions.get_biochem_rxns(C6H12O6, H2O)
eqs = df[!,:Equation]
eq = first(eqs) # I happen to know this one has no prefixed stoich data
reactants, products, rstoich, pstoich = PubChemReactions.parse_rhea_equation(eq)
rxn = Reaction(1, reactants, products, rstoich, pstoich; only_use_rate=true) # arbitrarily assigning constant rate ()

@test rxn isa Catalyst.Reaction

eq = eqs[19] # includes stoich 
rxn = Reaction(1, PubChemReactions.parse_rhea_equation(eq)...; only_use_rate=true) # arbitrarily assigning constant rate ()
@test rxn isa Catalyst.Reaction
