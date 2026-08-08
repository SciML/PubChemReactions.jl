"""
    Compound(name, cid, json, json_view)

Metadata container for a PubChem compound record attached to a symbolic species.

`Compound` is attached internally while constructing a PubChem species. Use
accessors such as [`get_cid`](@ref), [`get_name`](@ref), and [`get_graph`](@ref)
instead of depending on its storage representation.

# Fields
- `name::String`: PubChem record title.
- `cid::Int`: PubChem compound identifier.
- `json::AbstractDict{Symbol, Any}`: PUG compound record used to build structural metadata.
- `json_view::AbstractDict{Symbol, Any}`: PUG-View record used to retrieve displayed
  properties.

# Examples
```julia
using SymbolicUtils: getmetadata

species = PubChemReactions.search_compound("water")
compound = getmetadata(species, PubChemReactions.Compound)
compound.cid
```
"""
struct Compound
    name::String
    cid::Int
    json::AbstractDict{Symbol, Any}
    json_view::AbstractDict{Symbol, Any}
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
