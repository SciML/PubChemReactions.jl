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
