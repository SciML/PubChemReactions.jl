"uses the 2D image provided by PubChem, rather than trying to use the graph data"
function atomplot(s)
    cid = get_cid(s)
    io = IOBuffer()
    url = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=$(cid)&t=l"
    Downloads.download(url, io)
    sname = get_name(s)
    mf = get_molecular_formula(s)
    title = "$cid: $mf | $sname"
    @info (cid=cid, name=sname, formula=mf)
    p = Plots.plot(load(io))
    title!(p, title)
    display(p)
end
