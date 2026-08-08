# Usage

PubChemReactions attaches PubChem records and atom-bond graphs to Symbolics species.
Use its `@species` macro when building Catalyst reaction systems so the resulting symbolic
variables retain the metadata required by the balancing and plotting helpers.

## Create a species

```julia
using Catalyst: Reaction
using Symbolics: @variables
using PubChemReactions
using PubChemReactions: @species

@variables t
@species water(t) [cid = 962]
```

The macro fetches the PubChem record for the supplied name or `cid`, attaches the metadata,
and returns an ordinary Symbolics variable. The metadata can then be inspected through the
stable accessors.

```julia
get_cid(water)
get_name(water)
get_charge(water)
atom_counts(water)
```

## Check a reaction

Create a Catalyst reaction from metadata-bearing species, then check its elemental balance
or derive its symbolic balance equations.

```julia
rxn = Reaction(1, [water], [water])
isbalanced(rxn)
balance_eqs(rxn)
```

Network access is required when records are not saved locally. Call `save_species` to cache
a species under `COMPOUNDS_DIR`; `load_species` constructs a species that reads that cache
when it is next resolved.
