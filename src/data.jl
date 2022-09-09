pug_url_name(cname) = joinpath(PUG_URL, "compound/name/$(cname)/record/JSON")
pug_view_url_name(cname) = joinpath(PUG_VIEW_URL, "compound/name/$(cname)/JSON")

pug_url_cid(cid) = joinpath(PUG_URL, "compound/cid/$(cid)/record/JSON")
pug_view_url_cid(cid) = joinpath(PUG_VIEW_URL, "compound/cid/$(cid)/JSON")

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

function get_json_and_view_from_cname(cname::AbstractString; kwargs...)
    cname = HTTP.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/record/JSON/"#?record_type=3d" # FIX
    get_json_and_view(input_url; kwargs...)
end

function get_json_and_view_from_cid(cid; kwargs...)
    cid = HTTP.escapeuri(cid)
    input_url = "$(PUG_URL)/compound/cid/$(cid)/record/JSON/"#?record_type=3d" # FIX
    get_json_and_view(input_url; kwargs...)
end

compound_url(cid::AbstractString) = joinpath(PC_ROOT, "compound/$cid")
compound_url(s) = compound_url(string(get_cid(s)))
compound_fns(cid) = joinpath(COMPOUNDS_DIR, string(cid)) .* ("/pug.json", "/pug_view.json")

function load_json_and_view_from_cid(cid)
    fns = compound_fns(cid)
    # all(isfile.(fns))
    JSON3.read.(read.(fns))
end

function get_json_and_view(input_url; verbose=false)
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
    atom_pairs = compound.atoms.aid .=> compound.atoms.element
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
    Reaction(1, x[1], x[2])
    # balance(x[1], x[2])
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

"I may want to just compute a molecular formula from the graph"
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
val_from_sec(s) = s.Information[1].Value.StringWithMarkup[1].String

function get_smiles(s)
    jv = get_jview(s)
    for sec in jv.Section
        if sec.TOCHeading == "Names and Identifiers"
            for sec2 in sec.Section
                if sec2.TOCHeading == "Computed Descriptors"
                    for sec3 in sec2.Section
                        if sec3.TOCHeading == "Canonical SMILES"
                            return sec3.Information[1].Value.StringWithMarkup[1].String
                        end
                    end
                end
            end
        end
    end
    error("not found")
end

function get_chebi_id(csym)
    cid = get_cid(csym)
    input_url = "$PUG_VIEW_URL/data/compound/$cid/JSON/?heading=Biochemical+Reactions"
    res = HTTP.get(input_url)
    JSON3.read(String(res.body))[:Record][:Reference][1][:SourceID]
end

"fix this, its misleading because it doesn't return a bool"
function is_mass_conserved(rxn)
    netmass(rxn) == 0
end

function netmass(rxn)
    subs, prods = rxn_masses(rxn)
    prods - subs
end

function rxn_masses(rxn)
    subs = sum(rxn.substoich .* get_mass.(rxn.substrates))
    prods = sum(rxn.prodstoich .* get_mass.(rxn.products))
    subs, prods
end

pubchem_search(s) = open_in_default_browser(joinpath(PC_ROOT, "#query=$s"))
