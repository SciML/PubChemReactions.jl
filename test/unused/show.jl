s = only(@species glucose(t))
io2 = copy(io)
p = tempname()
write(p, take!(io))

struct Species
    imgpath::Any
end
struct Species2
    s::Any
end

function Base.show(io::IO, ::MIME"text/markdown", s::Num)
    _, io_ = PubChemReactions.download_atomplot(s)
    pl = Plots.plot(load(io_); size = (200, 200))
    p = abspath(joinpath(tempname(), ".png"))
    savefig(pl, p)
    @info p
    # save(p, pl)
    return print(io, "![$(get_name(s))]($(p))")
end

function Base.show(io, ::MIME"text/markdown", s::Species)
    return print(io, "![]($(s.imgpath))")
end

function Base.show(io::IO, ::MIME"text/markdown", s2::Species2)
    s = s2.s
    cid, io_ = PubChemReactions.download_atomplot(s)
    pl = Plots.plot(load(io_); size = (200, 200), background_color = :transparent)
    p = abspath(cid, ".png")
    savefig(pl, p)
    @info p
    # save(p, pl)
    return print(io, "![]($(p))")
end

cid, io = PubChemReactions.download_atomplot(s)
# p = "glucose.png"
p = abspath("glucose.png")
using ImageIO, FileIO
pl = Plots.plot(load(io))
save(p, pl)

s2 = Species("/Users/anand/.julia/dev/PubChemReactions.jl/glucose.png")
s3 = Species(p)

s = only(@species glucose(t))
s
typeof(s)

v = @variables x
s2 = Species2(s)
s2
x

function Base.show(io::IO, ::MIME"text/markdown", s::Num)
    _, io_ = PubChemReactions.download_atomplot(s)
    pl = Plots.plot(load(io_); size = (200, 200))
    p = abspath(joinpath(tempname(), ".png"))
    savefig(pl, p)
    @info p
    # save(p, pl)
    return print(io, "![$(get_name(s))]($(p))")
end

s = only(@species glucose(t))
s
methods(show, [Any, MIME, Num])
