"""
    Compound(name, cid, json, json_view)

Metadata container for a PubChem compound record attached to a symbolic species.

Stores the compound name, PubChem compound identifier, PUG record JSON, and
PUG-View JSON for later property queries and graph construction.

# Fields

- `name::String`: Compound name from the PubChem PUG-View record.
- `cid::Int`: PubChem compound identifier.
- `json::JSON3.Object`: The compound record returned by the PubChem PUG API.
- `json_view::JSON3.Object`: The detailed record returned by the PubChem PUG-View API.

# Constructor

`Compound(name, cid, json, json_view)` constructs the metadata container from
the name, identifier, and the two PubChem JSON objects. The constructor stores
these values without modifying the JSON objects.

# Returns

A `Compound` value suitable for attaching to a symbolic species with
`SymbolicUtils.setmetadata`.

# Examples

```julia
using JSON3
using PubChemReactions

record = JSON3.read("{}")
compound = Compound("Water", 962, record, record)
compound.cid
```
"""
struct Compound
    name::String
    cid::Int
    json::JSON3.Object
    json_view::JSON3.Object
end

struct AtomBondGraph
    g::SimpleGraph
    atoms::Vector{Pair{Int, Int}}
end

struct CompoundCharge
    charge::Int
end

function build_atom_graph(n_vertices, bond_pairs)
    g = SimpleGraph(n_vertices)
    for bp in bond_pairs
        add_edge!(g, bp)
    end
    return g
end

reaction_species(rxn::Reaction) = unique(reduce(vcat, (rxn.substrates, rxn.products)))
reaction_species(rxns::Vector{<:Reaction}) = unique(reduce(vcat, map(reaction_species, rxns)))
