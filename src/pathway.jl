sdq_url(query) = string(SDQ_URL, "?infmt=json&outfmt=json&query=", HTTP.escapeuri(JSON3.write(query)))

"""
    pathway_record(pwacc)

Fetch the PubChem `pathway` collection record for accession `pwacc`, given exactly
as PubChem lists it (`"<source>:<external id>"`, case-sensitive, e.g.
`"Reactome:R-HSA-70171"`). Throws `ArgumentError` for malformed accessions or when
PubChem has no pathway with exactly that accession. Internal, not exported.
"""
function pathway_record(pwacc)
    pwacc isa AbstractString && occursin(r"^[^:\s]+:\S+$", pwacc) ||
        throw(ArgumentError("invalid pathway accession $(repr(pwacc))"))
    query = (select = "*", collection = "pathway", where = (ands = [(pwacc = pwacc,)],), start = 1, limit = 10)
    rows = only(JSON3.read(get_page(sdq_url(query))).SDQOutputSet).rows
    matches = filter(r -> r.pwacc == pwacc, rows)
    length(matches) == 1 ||
        throw(ArgumentError("PubChem has no pathway with accession $(repr(pwacc))"))
    return only(matches)
end

function pathway_json(pid)
    pw = pathway_record(pid)
    query = (
        download = "*", collection = "pathwayreaction",
        where = (ands = [(source = pw.source,), (externalid = pw.externalid,)],),
        order = ["relevancescore,desc"], start = 1, limit = 10_000_000,
    )
    return JSON3.read(get_page(sdq_url(query)))
end

"""
need to fix rate handling, defaulting all to 1 sucks
"""
function get_pathway(pid)
    # `pathway_json` deserializes to a JSON3.Array of reaction objects; assert it so the
    # filter/broadcast below operate on `JSON3.Object` entries rather than the wide
    # `JSON3.read` return union (which would otherwise admit String/Char/Nothing branches).
    jrxns = pathway_json(pid)::JSON3.Array
    # `Reaction` rejects entries with neither reactants nor products, so drop them.
    valid = filter(jrxns) do jr
        rcids, pcids = pc_pathway_rxn_to_rp_cids(jr)
        !(isempty(rcids) && isempty(pcids))
    end
    return pathway_reaction_to_reaction.(valid)
end

function species_(s, cid)
    cid === nothing && return species_from_cid(s)
    return species_from_cid_and_name(s, cid)
end

function parse_pathway_reaction(rxn_str)
    subs, prods = rhea_to_reacts_prods(rxn_str)
    subs2 = map(cid_from_a_tag, subs)
    prods2 = map(cid_from_a_tag, prods)
    ss = map(x -> species_(x...), subs2)
    ps = map(x -> species_(x...), prods2)
    return ss, ps
end

"""
    pathway_reaction(rxn_str)

Parse a PubChem pathway reaction string into a Catalyst `Reaction` with unit rate.

# Arguments

- `rxn_str`: A PubChem pathway reaction string containing HTML links to the
  reactant and product compounds.

# Returns

A Catalyst `Reaction` with rate `1`, with reactants and products represented as
PubChem species. Compound records are retrieved while parsing the linked
species.
"""
pathway_reaction(rxn_str) = Reaction(1, parse_pathway_reaction(rxn_str)...)

function cid_from_a_tag(str)
    r = parsehtml(str).root
    a = eachmatch(Selector("a")::Cascadia.Selector, r)
    length(a) == 0 && return (str, nothing) # error("$str has no cid")
    a = only(a)
    t = Gumbo.text(a)
    url = a.attributes["href"]
    paths = splitpath(url)
    paths[end - 1] != "compound" && return (str, nothing) # not protein
    return t, paths[end]
end

function get_page(url)
    io = IOBuffer()
    Downloads.download(url, io; timeout = 300)
    return String(take!(io))
end

function pc_pathway_rxn_to_rp_cids(jr)
    rcids = get(jr, :cidsreactant, [])
    pcids = get(jr, :cidsproduct, [])
    return to_arr.((rcids, pcids))
end

to_arr(xs) = isa(xs, AbstractArray) ? xs : [xs]
to_arr(xs::AbstractString) = isempty(xs) ? String[] : [xs]

"""
check that the reaction str `jr.reaction`, when parsed, has the same length as the rcids and pcids from the json
"""
function is_reacts_prods_cids_aligned(jr)
    rxn_str = jr.reaction
    rp = pc_pathway_rxn_to_rp_cids(jr)
    sp = rhea_to_reacts_prods(rxn_str)
    return length.(rp) == length.(sp)
end

"""
todo: fix that all reactions are unidirectional
"""
function pathway_reaction_to_reaction(jr)
    rcids, pcids = pc_pathway_rxn_to_rp_cids(jr)
    rs = species_from_cid.(rcids)
    ps = species_from_cid.(pcids)
    # Catalyst's `Reaction` rejects untyped (`Vector{Any}`) substrate/product vectors.
    # Broadcasting over an empty cid list yields `Any[]`, so narrow empty sides to the
    # concrete species element type taken from whichever side has entries.
    ET = isempty(rs) ? (isempty(ps) ? Any : eltype(ps)) : eltype(rs)
    isempty(rs) && (rs = Vector{ET}())
    isempty(ps) && (ps = Vector{ET}())
    return Reaction(1, rs, ps)
end
