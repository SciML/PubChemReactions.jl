using OrdinaryDiffEq
using Catalyst
using PubChemReactions
using PubChemReactions: atom_counts
using PeriodicTable, Unitful
# https://en.wikipedia.org/wiki/Ammonia_production
# to_elements(aids) = map(i -> elements[i], aids)

t = only(@parameters t)

# steam reformation of H2
# function Reform()

@species CH4(t) [cid = 297, save = true, load = true]
@species H2O(t) [cid = 962, save = true, load = true]
@species CO(t) [cid = 281, save = true, load = true]
@species H2(t) [cid = 783, save = true, load = true]
@species O2(t) [cid = 977, save = true, load = true]
@species N2(t) [cid = 947, save = true, load = true]
@species NH3(t) [cid = 222, save = true, load = true]
@species CO2(t) [cid = 280, save = true, load = true]
@species N(t) [cid = 57370662, save = true, load = true]

subs, prods = ([N2, H2], [NH3])
rxn = Reaction(1, subs, prods)
ss = PubChemReactions.get_smiles.(subs)
ps = PubChemReactions.get_smiles.(prods)
ammonia_rxn = balance(rxn)
