struct SpeciesCid end
struct SpeciesName end
struct SpeciesSave end
struct SpeciesLoad end

Symbolics.option_to_metadata_type(::Val{:cid}) = SpeciesCid
Symbolics.option_to_metadata_type(::Val{:name}) = SpeciesName
Symbolics.option_to_metadata_type(::Val{:save}) = SpeciesSave
Symbolics.option_to_metadata_type(::Val{:load}) = SpeciesLoad

function catalyst_species(name::Symbol)
    t = Catalyst.default_t()
    return only(Catalyst.@species $name($t))
end

function set_species_metadata(s, j, jview)
    g, atom_pairs = compound_json_to_simplegraph(j)
    s = setmetadata(
        s,
        PubChemReactions.Compound,
        PubChemReactions.Compound(jview.Record.RecordTitle, jview.Record.RecordNumber, j, jview)
    )
    s = setmetadata(s, PubChemReactions.AtomBondGraph, PubChemReactions.AtomBondGraph(g, atom_pairs))
    s = setmetadata(
        s, PubChemReactions.CompoundCharge,
        PubChemReactions.CompoundCharge(j.PC_Compounds[1].charge)
    )
    return s
end

"""
    @species_str(cname)

Look up the PubChem compound named by the string literal `cname`, attach the
compound metadata to a symbolic species, and display its atom plot.

# Arguments
- `cname::String`: PubChem compound name supplied as a string literal.

# Returns
- A Symbolics variable carrying PubChem metadata.

# Examples
```julia
water = PubChemReactions.@species_str "water"
get_cid(water)
```
"""
macro species_str(cname)
    s = search_compound(cname)
    atomplot(s)
    return s
end

"""
    isspecies(s) -> Bool

Return `true` when `s` has PubChem atom-bond graph metadata.

# Arguments
- `s`: symbolic value to inspect.

# Returns
- `Bool`: `true` when `s` carries the package's atom-bond graph metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
isspecies(water)
```
"""
function isspecies(s)
    return hasmetadata(s, AtomBondGraph)
end

"""
    @species declarations

Create Catalyst species and attach PubChem compound, atom-bond graph, and charge metadata.
The declarations follow `Catalyst.@species` syntax and may use PubChem metadata options.

# Arguments
- `declarations`: one or more `Catalyst.@species` declarations. A declaration must provide a
  name that PubChem can resolve or a `cid` option.

# Metadata options
- `cid::Integer`: PubChem compound identifier used instead of a name lookup.
- `name::AbstractString`: PubChem compound name to use for the lookup.
- `save::Bool = false`: save retrieved PubChem records in the local compound cache.
- `load::Bool = false`: read records from the local compound cache when they are available.

# Returns
- Catalyst symbolic species with [`Compound`](@ref), atom-bond graph, and charge metadata.

# Throws
- `ErrorException`: if PubChem cannot resolve or retrieve a requested compound.

# Examples
```julia
using PubChemReactions: @species
using Symbolics: @variables

@variables t
@species water(t) [cid = 962]
```
"""
macro species(xs...)
    catalyst_species = Expr(
        :macrocall,
        GlobalRef(Catalyst, Symbol("@species")),
        __source__,
        xs...
    )
    return esc(_attach_pubchem_metadata(macroexpand(__module__, catalyst_species)))
end

function _attach_pubchem_metadata(expr)
    expr isa Expr && expr.head === :block || return expr
    args = map(expr.args) do arg
        if arg isa Expr && arg.head === :(=)
            lhs, rhs = arg.args
            Expr(:(=), lhs, Expr(:call, GlobalRef(PubChemReactions, :tospecies), rhs))
        else
            arg
        end
    end
    return Expr(:block, args...)
end

"""
    tospecies(s; jsons = nothing)

Attach PubChem metadata to a symbolic species during internal species construction.

# Arguments
- `s`: symbolic scalar, symbolic array, or collection to enrich with PubChem metadata.

# Keywords
- `jsons = nothing`: optional `(pug_record, pug_view_record)` tuple. When omitted,
  the record is downloaded or read from the local compound cache according to the
  metadata already attached to `s`.

# Returns
- A value with the same symbolic shape as `s` carrying PubChem compound, atom-bond,
  and charge metadata.
"""
function tospecies(s; jsons = nothing)
    return if s isa Symbolics.Arr
        Symbolics.wrap(tospecies(SymbolicUtils.unwrap(s)))
    elseif s isa AbstractArray
        map(tospecies, s)
    elseif SymbolicUtils.symtype(s) <: AbstractArray
        map(tospecies, collect(Symbolics.wrap(s)))
    else
        if hasmetadata(s, SpeciesName)
            cname = getmetadata(s, PubChemReactions.SpeciesName)
        else
            cname = string(SymbolicIndexingInterface.getname(s))
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
    search_compound(cname::AbstractString) -> Symbolics.Num

Retrieve the PubChem record for `cname` and return a symbolic species carrying its compound,
atom-bond graph, and charge metadata.

# Arguments
- `cname::AbstractString`: PubChem compound name to resolve.

# Returns
- `Symbolics.Num`: a Catalyst-compatible symbolic species with PubChem metadata.

# Throws
- `ErrorException`: if PubChem cannot resolve `cname` or return its compound records.

# Examples
```julia
water = search_compound("water")
get_cid(water)
```
"""
function search_compound(cname::AbstractString)
    csym = SymbolicUtils.unwrap(catalyst_species(Symbol(cname)))
    j, jview = PubChemReactions.get_json_and_view_from_cname(cname)
    return set_species_metadata(csym, j, jview)
end
# alias
species_from_name(cname) = search_compound(cname)

function species_from_cid(cid, j, jview)
    name = Symbol(jview.Record.RecordTitle)
    csym = SymbolicUtils.unwrap(catalyst_species(name))
    return set_species_metadata(csym, j, jview)
end

function species_from_cid(cid)
    j, jview = PubChemReactions.get_json_and_view_from_cid(cid)
    return species_from_cid(cid, j, jview)
end

function species_from_cid_and_name(name, cid; save = true, load = true)
    name = Symbol(name)
    return only(@species $name(Catalyst.default_t()) [save = save, cid = cid, load = load])
end

"""
    save_species(s; path = COMPOUNDS_DIR) -> String

Write the PubChem PUG and PUG-View JSON payloads attached to species `s`.

The payloads are stored under a compound-id subdirectory of `path`.

# Arguments
- `s`: PubChem species whose downloaded records will be saved.

# Keywords
- `path::AbstractString = COMPOUNDS_DIR`: directory that will contain the
  compound-id subdirectory.

# Returns
- `String`: newly-created compound directory, or a message identifying an existing cache.

# Throws
- `KeyError`: if `s` has no [`Compound`](@ref) metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
save_species(water; path = tempname())
```
"""
function save_species(s; path = COMPOUNDS_DIR)
    # isspecies(s) || error("$s is not a PubChemReactions species")
    cid = string(PubChemReactions.get_cid(s))

    # this would handle not requiring CID, but i dont like it
    p = joinpath(path, cid)
    isdir(p) && return "$s is already saved to $p"
    mkpath(p)

    c = getmetadata(s, Compound)
    pug_fn = joinpath(p, "pug.json")
    pug_view_fn = joinpath(p, "pug_view.json")
    open(pug_fn, "w") do io
        JSON.print(io, c.json)
    end
    open(pug_view_fn, "w") do io
        JSON.print(io, c.json_view)
    end
    return p
end

"""
    load_species(cid) -> Num

Load a saved PubChem species by compound identifier `cid`.

`cid` must identify a directory previously created by [`save_species`](@ref)
inside `COMPOUNDS_DIR`.

# Arguments
- `cid::Integer` or `cid::AbstractString`: PubChem compound identifier.

# Returns
- `Num`: a symbolic species configured to load its metadata from the local cache.

# Throws
- `SystemError`: if the expected cached JSON files do not exist.

# Examples
```julia
water = PubChemReactions.search_compound("water")
save_species(water)
cached_water = load_species(get_cid(water))
```
"""
function load_species(cid)
    _, jview = load_json_and_view_from_cid(cid)
    name = Symbol(jview.Record.RecordTitle)
    return species_from_cid_and_name(name, cid)
end
load_species(cid::Integer) = load_species(string(cid))

function open_compound(s::Num)
    isspecies(s) || error("$s not species")
    open_in_default_browser(compound_url(s))
    return nothing
end
