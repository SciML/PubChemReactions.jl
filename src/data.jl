pug_url_name(cname) = joinpath(PUG_URL, "compound/name/$(cname)/record/JSON")
pug_view_url_name(cname) = joinpath(PUG_VIEW_URL, "compound/name/$(cname)/JSON")

pug_url_cid(cid) = joinpath(PUG_URL, "compound/cid/$(cid)/record/JSON")
pug_view_url_cid(cid) = joinpath(PUG_VIEW_URL, "compound/cid/$(cid)/JSON")

# houses all the data and getters
parse_json(json::AbstractString) = JSON.parse(json; dicttype = JSON.Object{Symbol, Any})
parse_json(json::AbstractVector{UInt8}) = parse_json(String(json))

function get_cids_from_cname(cname::AbstractString; verbose = false)
    cname = URIs.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/cids/JSON"
    verbose && @info input_url
    return parse_json(get_page(input_url))
end

function get_cids(cname::AbstractString)
    identifier_list = get_cids_from_cname(cname)[:IdentifierList]
    cids = identifier_list[:CID]
    return convert(Vector{Int}, cids)
end

function get_json_from_cname(cname::AbstractString; verbose = false)
    cname = URIs.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/record/JSON/" #?record_type=3d"
    verbose && @info input_url
    return parse_json(get_page(input_url))
end

function get_json_and_view_from_cname(cname::AbstractString; kwargs...)
    cname = URIs.escapeuri(cname)
    input_url = "$(PUG_URL)/compound/name/$(cname)/record/JSON/" #?record_type=3d" # FIX
    return get_json_and_view(input_url; kwargs...)
end

function get_json_and_view_from_cid(cid; kwargs...)
    cid = URIs.escapeuri(cid)
    input_url = "$(PUG_URL)/compound/cid/$(cid)/record/JSON/" #?record_type=3d" # FIX
    return get_json_and_view(input_url; kwargs...)
end

compound_url(cid::AbstractString) = joinpath(PC_ROOT, "compound/$cid")
compound_url(s) = compound_url(string(get_cid(s)))
compound_fns(cid) = joinpath(COMPOUNDS_DIR, string(cid)) .* ("/pug.json", "/pug_view.json")

function load_json_and_view_from_cid(cid)
    fns = compound_fns(cid)
    # all(isfile.(fns))
    return parse_json.(read.(fns))
end

function get_json_and_view(input_url; verbose = false)
    verbose && @info input_url
    j = parse_json(get_page(input_url))
    cid = j.PC_Compounds[1].id.id.cid
    input_url2 = "$(PUG_VIEW_URL)/data/compound/$(cid)/JSON"
    j2 = parse_json(get_page(input_url2))
    return j, j2
end

function get_json_from_cid(cid; verbose = false)
    input_url = "$(PUG_URL)/compound/cid/$(cid)/record/JSON"
    verbose && @info input_url
    return parse_json(get_page(input_url))
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
    return g, atom_pairs
end

"""
    get_graph(s)

Return the atom-bond graph metadata attached to PubChem species `s`.

# Arguments
- `s`: symbolic species created by `PubChemReactions.@species` or
  `search_compound`.

# Returns
- Atom-bond graph metadata for `s`.

# Throws
- `ErrorException`: if `s` has no PubChem atom-bond graph metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
graph = get_graph(water)
```
"""
get_graph(s) = isspecies(s) ? getmetadata(s, AtomBondGraph) : error("no graph for var $s")

"""
    get_charge(s) -> Int

Return the formal charge stored in the PubChem species metadata for `s`.

# Arguments
- `s`: symbolic species created by `PubChemReactions.@species` or
  `search_compound`.

# Returns
- `Int`: formal molecular charge reported by PubChem.

# Throws
- `ErrorException`: if `s` has no PubChem charge metadata.

# Examples
```julia
hydron = PubChemReactions.search_compound("hydron")
get_charge(hydron)
```
"""
function get_charge(s)
    return isspecies(s) ? getmetadata(s, CompoundCharge).charge : error("no charge for var $s")
end

function get_reaction(eq)
    x = parse_rhea_equation(eq)
    return Reaction(1, x[1], x[2])
    # balance(x[1], x[2])
end

"""
    get_cid(s) -> Int

Return the PubChem compound identifier attached to species `s`.

# Arguments
- `s`: symbolic species created by `PubChemReactions.@species` or
  `search_compound`.

# Returns
- `Int`: PubChem compound identifier.

# Throws
- `KeyError`: if `s` has no [`Compound`](@ref) metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
get_cid(water)
```
"""
function get_cid(s)
    return getmetadata(s, Compound).cid
end

function get_j(s)
    return getmetadata(s, Compound).json.PC_Compounds[1]
end

function get_jview(s)
    return getmetadata(s, Compound).json_view.Record
end

"""
    get_name(s) -> String

Return the PubChem compound name attached to species `s`.

# Arguments
- `s`: symbolic species created by `PubChemReactions.@species` or
  `search_compound`.

# Returns
- `String`: PubChem record title.

# Throws
- `KeyError`: if `s` has no [`Compound`](@ref) metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
get_name(water)
```
"""
function get_name(s)
    return getmetadata(s, Compound).name
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

"""
    get_molecular_formula(s) -> String

Return the molecular formula recorded in PubChem's PUG-View metadata for species `s`.

# Arguments
- `s`: PubChem species carrying [`Compound`](@ref) metadata.

# Returns
- `String`: the molecular formula reported by PubChem.

# Throws
- `ErrorException`: if the record has no molecular-formula section.
"""
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
    j = parse_json(get_page(input_url))
    record = j[:Record]
    reference = record[:Reference]
    return reference[1][:SourceID]
end

"""
    is_mass_conserved(rxn::Catalyst.Reaction) -> Bool

Return whether the exact masses of the substrates and products of `rxn` are equal.

# Arguments
- `rxn::Catalyst.Reaction`: reaction whose species carry PubChem compound metadata.

# Returns
- `Bool`: `true` when the signed mass difference is zero.
"""
function is_mass_conserved(rxn)
    return netmass(rxn) == 0
end

function netmass(rxn)
    subs, prods = rxn_masses(rxn)
    return prods - subs
end

function rxn_masses(rxn)
    subs = sum(rxn.substoich .* get_mass.(rxn.substrates))
    prods = sum(rxn.prodstoich .* get_mass.(rxn.products))
    return subs, prods
end

"""
    pubchem_search(query) -> Bool

Open a PubChem search for `query` in the system default browser.

# Arguments
- `query`: value converted to text and used as the PubChem search query.

# Returns
- `Bool`: whether the browser command was started successfully.

# Examples
```julia
pubchem_search("citric acid")
```
"""
pubchem_search(s) = open_in_default_browser(joinpath(PC_ROOT, "#query=$s"))
