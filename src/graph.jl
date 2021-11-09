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
    g
end

Catalyst.species(rxn::Reaction) = unique(reduce(vcat, (rxn.substrates, rxn.products)))
Catalyst.species(rxns::Vector{<:Reaction}) = unique(reduce(vcat, map(species, rxns)))

"makes a Num with the Compound metadata"
function search_compound(cname)
    j, jview = PubChemReactions.get_json_and_view_from_cname(cname) #; verbose=true);
    g, atom_pairs = compound_json_to_simplegraph(j)
    csym = Symbol(cname)
    csym = Symbolics.unwrap(first(@variables $csym(Catalyst.DEFAULT_IV)))
    csym = setmetadata(csym, PubChemReactions.Compound, PubChemReactions.Compound(jview.Record.RecordTitle, jview.Record.RecordNumber, j, jview))
    csym = setmetadata(csym, PubChemReactions.AtomBondGraph, PubChemReactions.AtomBondGraph(g, atom_pairs))
    csym = setmetadata(csym, PubChemReactions.CompoundCharge, PubChemReactions.CompoundCharge(j.PC_Compounds[1].charge))
    csym
end

"generate a chemical species and plot it"
macro species_str(cname)
    s = search_compound(cname)
    atomplot(s)
    s
end

function isspecies(s)
    hasmetadata(s, AtomBondGraph)
end
