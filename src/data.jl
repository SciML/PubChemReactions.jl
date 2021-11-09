# houses all the data and getters
function get_cids_from_cname(cname::AbstractString; verbose=false)
    cname = HTTP.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/cids/JSON"
    verbose && @info input_url
    res = HTTP.get(input_url)
    if res.status == 200
        return JSON3.read(String(res.body))
    else
        error("Cannot Find CID of the species $cname.")
    end
end

function get_cids(cname::AbstractString)
    convert(Vector{Int}, get_cids_from_cname(cname)[:IdentifierList][:CID])
end

function get_chebi_id(csym)
    cid = getmetadata(csym, Compound).cid
    input_url = "$PUG_VIEW_URL/data/compound/$cid/JSON/?heading=Biochemical+Reactions"
    res = HTTP.get(input_url)
    if res.status == 200
        species_res = JSON3.read(String(res.body))
        chebi_id = species_res[:Record][:Reference][1][:SourceID]
        chebi_id
    else
        error("Cannot find CHEBI ID from $csym.")
    end
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

get_graph(s) = isspecies(s) ? getmetadata(s, AtomBondGraph) : error("no graph for var $s")
get_charge(s) = isspecies(s) ? getmetadata(s, CompoundCharge).charge : error("no charge for var $s")

function get_reaction(eq)
    x = parse_rhea_equation(eq)
    get_balanced_reaction(x[1], x[2])
end

function get_cid(s)
    getmetadata(s, Compound).cid
end

function get_j(s)
    getmetadata(s, Compound).json.PC_Compounds[1]
end

function get_jview(s)
    getmetadata(s, Compound).json_view.Record
end

function get_name(s)
    getmetadata(s, Compound).name
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

function get_molecular_formula(s)
    jv = get_jview(s)
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

"fix this, its misleading because it doesn't return a bool"
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