function pathway_json(pid)
    internal_pid = get_internal_pathwayid(pid)
    url = join(
        [
            PubChemReactions.RXN_TABLE_BASE_URL,
            "?infmt=json&outfmt=json&query={%22download%22:%22*%22,%22collection%22:%22pathwayreaction%22,%22where%22:{%22ands%22:[{%22pathwayid%22:%22$(internal_pid)%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22$(internal_pid)_pathwayreaction%22}",
        ]
    )
    return JSON3.read(get_page(url))
end

"""
need to fix rate handling, defaulting all to 1 sucks
"""
function get_pathway(pid)
    jrxns = pathway_json(pid)
    return pathway_reaction_to_reaction.(jrxns)
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

pathway_reaction(rxn_str) = Reaction(1, parse_pathway_reaction(rxn_str)...)

function cid_from_a_tag(str)
    r = parsehtml(str).root
    a = eachmatch(Selector("a"), r)
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
    Downloads.download(url, io)
    return String(take!(io))
end

function get_html(url)
    p = get_page(url)
    return parsehtml(p)
end

"""
Until PubChem gets back to me about how to query the reactions for a pathway, this is the only way it seems.
"""
function get_internal_pathwayid(pid)
    url = joinpath(PC_ROOT, "pathway/$(pid)")
    h = get_html(url).root
    ms = eachmatch(Selector("meta"), h)
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
    return Reaction(1, rs, ps)
end
