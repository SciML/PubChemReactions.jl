struct Compound2
    name::String
    cid::Int
    g::SimpleGraph
    atom_pairs::Vector{Pair{Int, Int}}
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
    bonds = compound.bonds
    bond_pairs = bonds.aid1 .=> bonds.aid2
    g = SimpleGraph(length(atom_pairs))
    for bp in bond_pairs
        add_edge!(g, bp)
    end
    g, atom_pairs
end

function symbolic_species_from_name(cname)
    j, jview = PubChemReactions.get_json_and_view_from_cname(cname) #; verbose=true);
    g, atom_pairs = compound_json_to_simplegraph(j)
    csym = Symbol(cname)
    csym = Symbolics.unwrap(first(@variables $csym(Catalyst.DEFAULT_IV)))
    csym = setmetadata(csym, PubChemReactions.Compound2, PubChemReactions.Compound2(jview.Record.RecordTitle, jview.Record.RecordNumber, g, atom_pairs))
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
    atom_counts(rxn.substrates, rxn.substoich) == atom_counts(rxn.products, rxn.prodstoich)
end

"check that the element counts in sub"
function isbalanced(rn::ReactionSystem)
    all(isbalanced.(reactions(rn)))
end

function countmap_(s)
    c = getmetadata(s, Compound2)
    aps = c.atom_pairs
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
