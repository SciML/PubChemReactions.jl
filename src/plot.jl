function download_atomplot(s)
    cid = get_cid(s)
    io = IOBuffer()
    url = joinpath(PC_ROOT, "image/imgsrv.fcgi?cid=$(cid)&t=l")
    Downloads.download(url, io)
    return cid, io
end

"""
    atomplot(s; verbose = false)

Download PubChem's 2D structure image for species `s` and display it with the
compound identifier, molecular formula, and name in the plot title.

# Arguments

- `s`: symbolic PubChem species.

# Keyword Arguments

- `verbose`: if `true`, log the compound identifier, name, and formula.

# Returns

The result of displaying the generated plot.

# Errors

Propagates PubChem download, image decoding, and metadata lookup errors.

# Examples

```julia
atomplot(species; verbose = true)
```
"""
function atomplot(s; verbose = false)
    cid, io = download_atomplot(s)
    sname = get_name(s)
    mf = get_molecular_formula(s)
    title = "$cid: $mf | $sname"
    verbose && @info(cid = cid, name = sname, formula = mf)
    p = Plots.plot(load(io))
    title!(p, title)
    return display(p)
end

# copied from https://github.com/fonsp/Pluto.jl/blob/6f5876228671f9d89d6f01cedc10221d83f012d6/src/webserver/WebServer.jl#L10-L33
function detectwsl()
    return Sys.islinux() &&
        isfile("/proc/sys/kernel/osrelease") &&
        occursin(r"Microsoft|WSL"i, read("/proc/sys/kernel/osrelease", String))
end

function open_in_default_browser(url::AbstractString)::Bool
    return try
        if Sys.isapple()
            Base.run(`open $url`)
            true
        elseif Sys.iswindows() || detectwsl()
            Base.run(`powershell.exe Start "'$url'"`)
            true
        elseif Sys.islinux()
            Base.run(`xdg-open $url`)
            true
        else
            false
        end
    catch ex
        false
    end
end

"""
    atomplot3d(s)

Open the PubChem 3D conformer viewer for species `s`.

# Arguments

- `s`: symbolic PubChem species.

# Returns

The `Bool` result from the platform browser launcher.

# Errors

Throws a metadata lookup error when `s` is not a PubChem species.

# Examples

```julia
atomplot3d(species)
```
"""
function atomplot3d(s)
    cid = string(get_cid(s))
    url = joinpath(PC_ROOT, "compound", "$cid#section=3D-Conformer&fullscreen=true")
    return open_in_default_browser(url)
end

"""
    atomplot2d(s)

Open the PubChem 2D structure viewer for species `s`.

# Arguments

- `s`: symbolic PubChem species.

# Returns

The `Bool` result from the platform browser launcher.

# Errors

Throws a metadata lookup error when `s` is not a PubChem species.

# Examples

```julia
atomplot2d(species)
```
"""
function atomplot2d(s)
    cid = string(get_cid(s))
    url = joinpath(PC_ROOT, "compound", "$cid#section=2D-Structure&fullscreen=true")
    return open_in_default_browser(url)
end
