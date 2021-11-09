
function atomplot(s)
    # g = getmetadata(s, AtomBondGraph).g
    ag = get_graph(s)

    g = ag.g
    mg = MetaGraph(g)
    j = s.metadata[PubChemReactions.Compound2].json
    coords = j.PC_Compounds[1].coords[1]["conformers"][1]
    as = ag.atoms
    as_d = Dict(as)
    vs = vertices(mg)
    for (i, v) in enumerate(vs)
        set_prop!(mg, i, :aid, as_d[i])
        set_prop!(mg, i, :x, coords["x"][i])
        set_prop!(mg, i, :y, coords["y"][i])
        # set_prop!(mg, i, :z, coords["z"][i])
    end
    
    function mylayout(g)
       return Point.(zip(coords[:x], coords[:y])) #, coords[:z]))
    end

    graphplot(g, layout=mylayout)
end

function atomplot2(s)
    gplot(get_graph(s).g)
end

"this is the best one, it just downloads and shows the image"
function atomplot3(s)
    cid = get_cid(s)
    io = IOBuffer()
    url = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=$(cid)&t=l"
    Downloads.download(url, io)
    # df = DataFrame()
    sname = get_name(s)
    mf = get_mf(s)
    title = "$cid: $mf | $sname"
    @info (cid=cid, name=sname, formula=mf)
    p = Plots.plot(load(io))
    title!(p, title)
    display(p)
end
