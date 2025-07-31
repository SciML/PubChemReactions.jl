# struct SpeciesCid end
# struct SpeciesName end
# struct SpeciesSave end
# struct SpeciesLoad end

# Symbolics.option_to_metadata_type(::Val{:cid}) = SpeciesCid
# Symbolics.option_to_metadata_type(::Val{:name}) = SpeciesName
# Symbolics.option_to_metadata_type(::Val{:save}) = SpeciesSave
# Symbolics.option_to_metadata_type(::Val{:load}) = SpeciesLoad

for metadata in [:cid, :name, :save, :load]
    struct_name = Symbol("Species", uppercasefirst(string(metadata)))
    @eval struct $struct_name end
    T = Val{metadata}
    @eval Symbolics.option_to_metadata_type(::$T) = $struct_name
end

function set_species_metadata(s, j, jview)
    g, atom_pairs = compound_json_to_simplegraph(j)
    s = setmetadata(s,
        PubChemReactions.Compound,
        PubChemReactions.Compound(jview.Record.RecordTitle, jview.Record.RecordNumber, j, jview))
    s = setmetadata(s, PubChemReactions.AtomBondGraph, PubChemReactions.AtomBondGraph(g, atom_pairs))
    s = setmetadata(s, PubChemReactions.CompoundCharge,
        PubChemReactions.CompoundCharge(j.PC_Compounds[1].charge))
    s
end

"""
generate a chemical species and plot it
"""
macro species_str(cname)
    s = search_compound(cname)
    atomplot(s)
    s
end

function isspecies(s)
    hasmetadata(s, AtomBondGraph)
end

"""
generate a chemical species.

@species CO(t) [cid=281]
@species Glucose(t)

todo
if the main thing is balancing, we probably dont have to request all the data, just the counts, since the jsons are massive
"""
macro species(xs...)
    Symbolics._parse_vars(:species,
        Real,
        xs,
        tospecies
    ) |> esc
end

"""
    tospecies(s::Sym)

Maps the variable to a species.
We probably want to do this async somehow
"""
function tospecies(s; jsons = nothing)
    if s isa Symbolics.Arr
        Symbolics.wrap(tospecies(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(tospecies, s)
    elseif Symbolics.symtype(s) <: AbstractArray
        Symbolics.recurse_and_apply(tospecies, s)
    else
        if hasmetadata(s, SpeciesName)
            cname = getmetadata(s, PubChemReactions.SpeciesName)
        else
            cname = string(Symbolics.getname(s))
        end
        if jsons === nothing
            if hasmetadata(s, SpeciesCid)
                cid = getmetadata(s, PubChemReactions.SpeciesCid)
                if hasmetadata(s, SpeciesLoad) && getmetadata(s, SpeciesLoad) &&
                   all(isfile.(compound_fns(cid)))
                    j, jview = PubChemReactions.load_json_and_view_from_cid(cid)
                else
                    @info "@species JSON GET CID: $cid"
                    j, jview = PubChemReactions.get_json_and_view_from_cid(cid)
                end
            else
                @info "@species JSON GET CNAME: $cname"
                j, jview = PubChemReactions.get_json_and_view_from_cname(cname)
                cid = jview.Record.RecordNumber
            end
        else
            j, jview = jsons
        end
        s = set_species_metadata(s, j, jview)
        if hasmetadata(s, SpeciesSave) && getmetadata(s, SpeciesSave)
            save_species(s)
        end
        s
    end
end
tospecies(s::Num; kwargs...) = Num(tospecies(Symbolics.value(s); kwargs...))

"""
makes a Num with the Compound metadata
"""
function search_compound(cname)
    csym = Symbol(cname)
    csym = Symbolics.unwrap(first(@variables $csym(Catalyst.DEFAULT_IV)))
    j, jview = PubChemReactions.get_json_and_view_from_cname(cname)
    set_species_metadata(csym, j, jview)
end
# alias
species_from_name(cname) = search_compound(cname)

function species_from_cid(cid, j, jview)
    name = Symbol(jview.Record.RecordTitle)
    csym = Symbolics.unwrap(first(@variables $name(Catalyst.DEFAULT_IV)))
    set_species_metadata(csym, j, jview)
end

function species_from_cid(cid)
    j, jview = PubChemReactions.get_json_and_view_from_cid(cid)
    species_from_cid(cid, j, jview)
end

function species_from_cid_and_name(name, cid; save = true, load = true)
    name = Symbol(name)
    only(@species $name(Catalyst.DEFAULT_IV) [save = save, cid = cid, load = load])
end

function save_species(s; path = COMPOUNDS_DIR)
    # isspecies(s) || error("$s is not a PubChemReactions species")
    cid = string(PubChemReactions.get_cid(s))

    # this would handle not requiring CID, but i dont like it
    # meta_fn = joinpath(datadir, "compounds.csv")
    # if isfile(meta_fn)
    #     meta = CSV.read()

    p = joinpath(path, cid)
    isdir(p) && return "$s is already saved to $p"
    mkpath(p)

    c = getmetadata(s, Compound)
    pug_fn = joinpath(p, "pug.json")
    pug_view_fn = joinpath(p, "pug_view.json")
    JSON3.write(pug_fn, c.json)
    JSON3.write(pug_view_fn, c.json_view)
    p
end

function load_species(cid)
    _, jview = load_json_and_view_from_cid(cid)
    name = Symbol(jview.Record.RecordTitle)
    species_from_cid_and_name(name, cid)
end
load_species(cid::Integer) = load_species(string(cid))

"""
maybe a bad idea overloading open
"""
function Base.open(s::Num)
    isspecies(s) || error("$s not species")
    open_in_default_browser(compound_url(s))
    nothing
end
