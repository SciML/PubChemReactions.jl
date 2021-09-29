# struct Compound2
#     name::String
#     cid::Int
#     g::SimpleGraph
#     atom_pairs::Vector{Pair{Int, Int}}
# end

struct Compound2
    name::String
    cid::Int
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
    input_url = "$(PUG_URL)/compound/name/$(cname)/record/JSON"
    verbose && @info input_url
    res = HTTP.get(input_url)
    if res.status == 200
        return JSON3.read(String(res.body))
    else
        error("Cannot Find CID of the species $cname.")
    end
end

function get_json_and_view_from_cname(cname::AbstractString; verbose=false)
    cname = HTTP.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/record/JSON"
    
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
    csym = setmetadata(csym, PubChemReactions.Compound2, PubChemReactions.Compound2(jview.Record.RecordTitle, jview.Record.RecordNumber))
    csym = setmetadata(csym, PubChemReactions.AtomBondGraph, PubChemReactions.AtomBondGraph(g, atom_pairs))
    csym = setmetadata(csym, PubChemReactions.CompoundCharge, PubChemReactions.CompoundCharge(j.PC_Compounds[1].charge))
    csym
end

rxnspecies(rxn::Reaction) = unique(reduce(vcat, (rxn.substrates, rxn.products)))
species_(rxns::Vector{<:Reaction}) = unique(reduce(vcat, map(rxnspecies, rxns)))

function parse_rhea_equation2(eq::AbstractString)
    reactants, products = PubChemReactions.rhea_to_reacts_prods(eq)
    rs = map(PubChemReactions.make_stoich_from_rhea, reactants)
    ps = map(PubChemReactions.make_stoich_from_rhea, products)
    rstoich, reactants = first.(rs), last.(rs)
    pstoich, products = first.(ps), last.(ps)

    symbolic_species_from_name.(reactants), symbolic_species_from_name.(products), rstoich, pstoich
end

"check that the element counts in substrates is equal to products"
function isbalanced(rxn)
    all(hasmetadata.(rxnspecies(rxn), Compound2)) || error("some species do not have atom graph metadata")
    all(hasmetadata.(rxnspecies(rxn), AtomBondGraph)) || error("some species do not have atom graph metadata")
    atom_counts(rxn.substrates, rxn.substoich) == atom_counts(rxn.products, rxn.prodstoich)
end

"check that the element counts in sub"
function isbalanced(rn::ReactionSystem)
    all(isbalanced.(reactions(rn)))
end

function countmap_(s)
    c = getmetadata(s, AtomBondGraph)
    aps = c.atoms
    countmap(last.(aps))
end

function atom_counts(speciess, stoichs)
    countmaps = countmap_.(speciess)

    for (stoich, cm) in zip(stoichs, countmaps)
        for (k, v) in cm 
            cm[k] = stoich * v
        end
    end

    mergewith(+, countmaps...)
end

"generate a chemical species"
macro species_str(cname)
    symbolic_species_from_name(cname)
end

function isspecies(s)
    hasmetadata(s, AtomBondGraph)
end

get_graph(s) = isspecies(s) ? getmetadata(s, AtomBondGraph) : error("no graph for var $s")
get_charge(s) = isspecies(s) ? getmetadata(s, CompoundCharge).charge : error("no charge for var $s")

# function graphplot(s::Symbolics.Symbolic)
#     g = getmetadata(s, AtomBondGraph).g
#     graphplot(g)
# end

elements(s) = unique(last.(get_graph(s).atoms))
elements(s::Vector) = Set(reduce(vcat, elements.(s)))

"""


# http://mathgene.usc.es/matlab-profs-quimica/reacciones.pdf

should i try to catch underdetermined soon, or just let LA give SingularException?


"""
function get_balanced_reaction(substrates, products)
    all_species = vcat(substrates, products)
    all(PubChemReactions.isspecies.(all_species)) || error("provide chemcial species (with graphs)")
    
    occuring_elements = collect(PubChemReactions.elements(all_species))
    atom_counts = PubChemReactions.countmap_.(all_species)
    charges = get_charge.(all_species)

    n_subs = length(substrates)
    n_prods = length(products)

    n_elems = length(occuring_elements)
    n_eqs = n_specs = length(all_species)
    
    A = zeros(Int, n_specs, n_specs)

    for i in 1:n_elems
        for j in 1:n_specs
            amt_of_i = occuring_elements[i]
            coeff = j > n_subs ? -1 : 1
            A[i, j] = haskey(atom_counts[j], amt_of_i) ? coeff * atom_counts[j][amt_of_i] : 0
        end
    end

    # check if we need another equation for charges
    if !all(charges .== 0) 
        for j in 1:n_specs # big hack, needs to be cleaned up 
            coeff = j > n_subs ? -1 : 1
            A[n_elems+1, j] = coeff * get_charge(all_species[j])
        end
    end

    A[end] = 1 # extra so not underdetermined

    b = zeros(Int, n_specs)
    b[end] = 1 # extra equation 

    x = A\b
    x .= x ./ minimum(x)
    x = round.(Int , x)

    Reaction(nothing, substrates, products, x[1:n_subs], x[n_subs+1:end])
end
