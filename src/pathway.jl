"need to fix rate handling, defaulting all to 1 sucks"
function get_pathway(pid)
    internal_pid = get_internal_pathwayid(pid)
    url = join([PubChemReactions.RXN_TABLE_BASE_URL, "?infmt=json&outfmt=json&query={%22download%22:%22*%22,%22collection%22:%22pathwayreaction%22,%22where%22:{%22ands%22:[{%22pathwayid%22:%22$(internal_pid)%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22$(internal_pid)_pathwayreaction%22}"])
    jrxns = JSON3.read(get_page(url))
    rxn_strs = getproperty.(jrxns, :reaction)
    rxns = Catalyst.Reaction[]
    failed = []
    for (i, rxn_str) in enumerate(rxn_strs)
        try
            subs, prods = subs_prods_from_pathway_rxn(rxn_str)

            push!(rxns, balance(subs, prods))
        catch e
            push!(failed, rxn_str)
            @info i, rxn_str, e
        end
    end
    rxns, failed
end

function subs_prods_from_pathway_rxn(rxn_str)
    subs, prods = rhea_to_reacts_prods(rxn_str)
    subs2 = map(cid_from_a_tag, subs)
    prods2 = map(cid_from_a_tag, prods)
    ss = map(x -> species_from_cid_and_name(x...), subs2)
    ps = map(x -> species_from_cid_and_name(x...), prods2)
    ss, ps
end

function cid_from_a_tag(str)
    r = parsehtml(str).root
    a = only(eachmatch(Selector("a"), r))
    # length(as) == 0 && error("$str has no cid")
    t = Gumbo.text(a)
    url = a.attributes["href"]
    paths = splitpath(url)
    @assert paths[end-1] == "compound" # not protein
    t, paths[end]
end

function get_page(url)
    io = IOBuffer()
    Downloads.download(url, io)
    String(take!(io))
end

function get_html(url)
    p = get_page(url)
    parsehtml(p)
end

"""
Until PubChem gets back to me about how to query the reactions for a pathway, this is the only way it seems.
"""
function get_internal_pathwayid(pid)
    url = joinpath(PC_ROOT, "pathway/$(pid)")
    h = get_html(url).root
    ms = eachmatch(Selector("meta"), h)
    m = only(filter(x -> haskey(x.attributes, "name") && (x.attributes["name"] == "pubchem_uid_value"), ms))
    m.attributes["content"]
end
