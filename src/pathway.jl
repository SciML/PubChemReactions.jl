function pathway_json(pid)
    internal_pid = get_internal_pathwayid(pid)
    url = join(
        [
            PubChemReactions.RXN_TABLE_BASE_URL,
            "?infmt=json&outfmt=json&query={%22download%22:%22*%22,%22collection%22:%22pathwayreaction%22,%22where%22:{%22ands%22:[{%22pathwayid%22:%22$(internal_pid)%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22$(internal_pid)_pathwayreaction%22}",
        ]
    )
    return parse_json(get_page(url))
end

"""
    get_pathway(pid::AbstractString) -> Vector{Catalyst.Reaction}

Retrieve the PubChem pathway identified by `pid` and convert its nonempty reaction records
to unit-rate Catalyst reactions.

# Arguments
- `pid::AbstractString`: PubChem or external pathway identifier accepted by the PubChem
  pathway endpoint.

# Returns
- `Vector{Catalyst.Reaction}`: reactions whose species carry PubChem metadata.

# Throws
- `ErrorException`: if PubChem cannot retrieve the pathway or its referenced compounds.

# Examples
```julia
reactions = PubChemReactions.get_pathway("Reactome:R-HSA-70171")
```
"""
function get_pathway(pid)
    # The PubChem endpoint returns reaction objects. Drop entries without either side before
    # creating Catalyst reactions, which reject an entirely empty reaction.
    jrxns = pathway_json(pid)
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
    pathway_reaction(rxn_str) -> Catalyst.Reaction

Parse a PubChem pathway reaction string into a Catalyst `Reaction` with unit rate.

# Arguments
- `rxn_str::AbstractString`: pathway reaction expression using a supported Rhea arrow
  and optional PubChem compound-id links.

# Returns
- `Catalyst.Reaction`: unit-rate reaction built from resolved PubChem species.

# Throws
- `ErrorException`: if PubChem records cannot be retrieved for a reaction species.

# Examples
```julia
pathway_reaction("water = water")
```
"""
pathway_reaction(rxn_str) = Reaction(1, parse_pathway_reaction(rxn_str)...)

function cid_from_a_tag(str)
    r = parsehtml(str).root
    a = eachmatch(Selector("a")::Selector, r)
    length(a) == 0 && return (str, nothing) # error("$str has no cid")
    a = only(a)
    t = Gumbo.text(a)
    url = a.attributes["href"]
    paths = splitpath(url)
    paths[end - 1] != "compound" && return (str, nothing) # not protein
    return t, paths[end]
end

function get_page(url)
    return _get_page(url, _request_page, (1, 2, 4, 8))
end

_request_page(url, io) = Downloads.request(url; output = io)

function _get_page(url, request, retry_delays)
    for delay in (retry_delays..., nothing)
        io = IOBuffer()
        response = request(url, io)
        200 <= response.status < 300 && return String(take!(io))

        retryable = response.status == 408 || response.status == 429 ||
            500 <= response.status < 600
        retryable && delay !== nothing || error("HTTP $(response.status) while requesting $url")
        sleep(delay)
    end
    return
end

function get_html(url)
    p = get_page(url)
    return parsehtml(p)
end

# PubChem's pathway table endpoint requires the site's internal pathway identifier.
function get_internal_pathwayid(pid)
    url = joinpath(PC_ROOT, "pathway/$(pid)")
    h = get_html(url).root
    ms = eachmatch(Selector("meta")::Selector, h)
    m = only(
        filter(
            x -> haskey(x.attributes, "name") &&
                (x.attributes["name"] == "pubchem_uid_value"), ms
        )
    )
    return m.attributes["content"]
end

function pc_pathway_rxn_to_rp_cids(jr)
    rcids = get(jr, :cidsreactant, [])
    pcids = get(jr, :cidsproduct, [])
    return to_arr.((rcids, pcids))
end

to_arr(xs) = isa(xs, AbstractArray) ? xs : [xs]

function is_reacts_prods_cids_aligned(jr)
    rxn_str = jr.reaction
    rp = pc_pathway_rxn_to_rp_cids(jr)
    sp = rhea_to_reacts_prods(rxn_str)
    return length.(rp) == length.(sp)
end

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
