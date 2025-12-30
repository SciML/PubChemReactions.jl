# PubChemReactions.jl

Generation of reaction networks from PubChem data


### Example

Let us suppose that you have two components, **Glucose** and **Water**, and you want to find potential BioChemical reactions & details associated with it.

We'll just need the name of the species we are interested and it will generate a symbolic variable with appropiate metadata:

```julia
using PubChemReactions, Catalyst
C6H12O6 = PubChemReactions.search_compound("glucose")
H2O = PubChemReactions.search_compound("water")
```

We can have some details about our components as well (returns a JSON Object currently):

```julia
julia> C6H12O6.metadata
Base.ImmutableDict{DataType, Any} with 4 entries:
  CompoundCharge => CompoundCharge(0)
  AtomBondGraph  => AtomBondGraph({24, 24} …
  Compound       => Compound("D-Glucose", 5793, {…
  VariableSource => (:variables, :glucose)
```

Now, let us find potential BioChemical reactions in which these species occur:

```julia
df = PubChemReactions.get_biochem_rxns(C6H12O6, H2O)
eqs = df[!, :Equation]
```

### Generating a Reaction Network using PubChem biochem equations

```julia
eq = first(eqs)
reactants, products, rstoich, pstoich = PubChemReactions.parse_rhea_equation(eq)

# arbitrarily assigning constant rate of 1. 
# this is where we'd like to do some lookup on the reactants and products to make an educated guess about the rate law
rxn = Reaction(1, reactants, products, rstoich, pstoich; only_use_rate = true)
```
