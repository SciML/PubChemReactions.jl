module PubChemReactions

using JSON3, HTTP, Symbolics, CSV, DataFrames
using Catalyst
using Symbolics:variable
using StatsBase, LightGraphs
# using GraphMakie, GLMakie

const PUG_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
const PUG_VIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
const RHEA_URL = "https://www.rhea-db.org/rhea"

struct Compound
    name::String
    cids::Vector{Int}
end

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

function species_info(csym)
    meta = getmetadata(csym, Compound)
    cid = meta.cids[1]
    input_url = "$(PUG_URL)/compound/cid/$(cid)/description/JSON"
    res = HTTP.get(input_url)
    if res.status == 200
        species_res = JSON3.read(String(res.body))
        return species_res
    else
        error("Cannot find description of the species $csym")
    end
end

function get_chebi_id(csym)
    cid = getmetadata(csym, Compound).cids[1]
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

function gen_sym(cname)
    csym = Symbol(cname)
    csym = Symbolics.unwrap(first(@variables $csym(Catalyst.DEFAULT_IV)))
    csym = setmetadata(csym, Compound, Compound(cname, get_cids(cname)))
    csym
end

function get_biochem_rxns(csym, csyms...)
    chebi_ids = get_chebi_id(csym)
    for c in csyms
        chebi_id = get_chebi_id(c)
        chebi_ids = chebi_ids * "+" * chebi_id
    end

    input_url = "$RHEA_URL/?query=$(chebi_ids)&columns=rhea-id,equation,chebi-id&format=tsv"
    res = HTTP.get(input_url)
    if res.status == 200
        return CSV.read(IOBuffer(res.body), DataFrame)
    else
        error("Cannot find Biochemical reactions")
    end
end

"includes stoich values" 
function rhea_to_reacts_prods(eq::AbstractString)
    lhs, rhs = split(eq, " = ")
    split(lhs, " + "), split(rhs, " + ")
end

function parse_rhea_equation(eq::AbstractString)
    reactants, products = rhea_to_reacts_prods(eq)
    rs = map(make_stoich_from_rhea, reactants)
    ps = map(make_stoich_from_rhea, products)
    rstoich, reactants = first.(rs), last.(rs)
    pstoich, products = first.(ps), last.(ps)

    gen_sym.(reactants), gen_sym.(products), rstoich, pstoich
end

function make_stoich_from_rhea(s)
    if startswith(s, r"(\d).* ")
        ss = split(s, " ")
        parse(Int, ss[1]), ss[2]
    else 
        1, s
    end
end

include("graph.jl")

export Compound

export @species_str, AtomBondGraph, CompoundCharge

export get_graph, get_charge

end # module
