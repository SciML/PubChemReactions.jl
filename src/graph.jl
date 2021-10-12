# struct Compound2
#     name::String
#     cid::Int
#     g::SimpleGraph
#     atom_pairs::Vector{Pair{Int, Int}}
# end

struct Compound2
    name::String
    cid::Int
    json::JSON3.Object
    json_view::JSON3.Object
    # g::SimpleGraph
    # atom_pairs::Vector{Pair{Int, Int}}
end

struct AtomBondGraph
    g::SimpleGraph
    atoms::Vector{Pair{Int, Int}}
end

struct CompoundCharge
    charge::Int
end

function get_json_from_cname(cname::AbstractString; verbose=false)
    cname = HTTP.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/record/JSON/"#?record_type=3d"
    verbose && @info input_url
    res = HTTP.get(input_url)
    if res.status == 200
        return JSON3.read(String(res.body))
    else
        error("Cannot Find CID of the species $cname.")
    end
end

macro pc_str(x)
    get_json_from_cname(x)
end

function get_json_and_view_from_cname(cname::AbstractString; verbose=false)
    cname = HTTP.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/record/JSON/"#?record_type=3d" # FIX
    
    verbose && @info input_url
    res = HTTP.get(input_url)
    if res.status == 200
        j = JSON3.read(String(res.body))
        cid = j.PC_Compounds[1].id.id.cid
        input_url2 = "$(PUG_VIEW_URL)/data/compound/$(cid)/JSON"
        j2 = JSON3.read(String(HTTP.get(input_url2).body))
        return j, j2
    else
        error("Cannot Find CID of the species $cname.")
    end
end

macro pc2_str(x)
    get_json_and_view_from_cname(x)
end

macro reaction_str(x)
    get_reaction(x)
end

function get_json_from_cid(cid; verbose=false)
    input_url = "$(PUG_URL)/compound/cid/$(cid)/record/JSON"
    verbose && @info input_url
    res = HTTP.get(input_url)
    if res.status == 200
        return JSON3.read(String(res.body))
    else
        error("Cannot Find CID with id $cid.")
    end
end

function compound_json_to_simplegraph(j)
    compound = j.PC_Compounds[1]
    atom_pairs = compound.atoms.aid .=>compound.atoms.element
    if haskey(compound, :bonds) 
        bonds = compound.bonds
        bond_pairs = bonds.aid1 .=> bonds.aid2
    else  
        bond_pairs = []
    end
    g = build_atom_graph(length(atom_pairs), bond_pairs)
    g, atom_pairs
end

function build_atom_graph(n_vertices, bond_pairs)
    g = SimpleGraph(n_vertices)
    for bp in bond_pairs
        add_edge!(g, bp)
    end
    g
end

function symbolic_species_from_name(cname)
    j, jview = PubChemReactions.get_json_and_view_from_cname(cname) #; verbose=true);
    g, atom_pairs = compound_json_to_simplegraph(j)
    csym = Symbol(cname)
    csym = Symbolics.unwrap(first(@variables $csym(Catalyst.DEFAULT_IV)))
    csym = setmetadata(csym, PubChemReactions.Compound2, PubChemReactions.Compound2(jview.Record.RecordTitle, jview.Record.RecordNumber, j, jview))
    csym = setmetadata(csym, PubChemReactions.AtomBondGraph, PubChemReactions.AtomBondGraph(g, atom_pairs))
    csym = setmetadata(csym, PubChemReactions.CompoundCharge, PubChemReactions.CompoundCharge(j.PC_Compounds[1].charge))
    csym
end

get_species(s) = symbolic_species_from_name(s)

Catalyst.species(rxn::Reaction) = unique(reduce(vcat, (rxn.substrates, rxn.products)))
Catalyst.species(rxns::Vector{<:Reaction}) = unique(reduce(vcat, map(species, rxns)))

function parse_rhea_equation(eq::AbstractString)
    reactants, products = PubChemReactions.rhea_to_reacts_prods(eq)
    rs = map(PubChemReactions.make_stoich_from_rhea, reactants)
    ps = map(PubChemReactions.make_stoich_from_rhea, products)
    rstoich, reactants = first.(rs), last.(rs)
    pstoich, products = first.(ps), last.(ps)

    symbolic_species_from_name.(reactants), symbolic_species_from_name.(products), rstoich, pstoich
end

"generate a chemical species"
macro species_str(cname)
    s = symbolic_species_from_name(cname)
    atomplot3(s)
    s
end

function isspecies(s)
    hasmetadata(s, AtomBondGraph)
end

get_graph(s) = isspecies(s) ? getmetadata(s, AtomBondGraph) : error("no graph for var $s")
get_charge(s) = isspecies(s) ? getmetadata(s, CompoundCharge).charge : error("no charge for var $s")

elements(s) = unique(last.(get_graph(s).atoms))
elements(s::Vector) = Set(reduce(vcat, elements.(s)))


function get_reaction(eq)
    x = parse_rhea_equation(eq)
    get_balanced_reaction(x[1], x[2])
end

function get_cid(s)
    getmetadata(s, Compound2).cid
end

function get_j(s)
    getmetadata(s, Compound2).json.PC_Compounds[1]
end

function get_jv(s)
    getmetadata(s, Compound2).json_view.Record
end

function get_name(s)
    getmetadata(s, Compound2).name
end

function get_mass(s)
    j = get_j(s)
    for p in j.props
        if p.urn["label"] == "Mass" && p.urn["name"] == "Exact"
            return parse(Float64, p.value["sval"])
        end
    end 
    error("not found")
end

function get_mf(s)
    jv = get_jv(s)
    for sec in jv.Section
        if sec.TOCHeading == "Names and Identifiers"
            for sec2 in sec.Section
                if sec2.TOCHeading == "Molecular Formula"
                    return sec2.Information[1].Value.StringWithMarkup[1].String
                end
            end
        end
    end 
    error("not found")
end

" fix this, its misleading because it doesn't return "
function is_mass_conserved(rxn)
    subs = sum(rxn.substoich .* get_mass.(rxn.substrates))
    prods = sum(rxn.prodstoich .* get_mass.(rxn.products))
    subs, prods
end


function netmass(rxn)
    subs = sum(rxn.substoich .* get_mass.(rxn.substrates))
    prods = sum(rxn.prodstoich .* get_mass.(rxn.products))
    prods - subs
end