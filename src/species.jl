struct SpeciesCid end
struct SpeciesName end

Symbolics.option_to_metadata_type(::Val{:cid}) = SpeciesCid
Symbolics.option_to_metadata_type(::Val{:name}) = SpeciesName

function set_species_metadata(s, j, jview)
    g, atom_pairs = compound_json_to_simplegraph(j)
    s = setmetadata(s, PubChemReactions.Compound, PubChemReactions.Compound(jview.Record.RecordTitle, jview.Record.RecordNumber, j, jview))
    s = setmetadata(s, PubChemReactions.AtomBondGraph, PubChemReactions.AtomBondGraph(g, atom_pairs))
    s = setmetadata(s, PubChemReactions.CompoundCharge, PubChemReactions.CompoundCharge(j.PC_Compounds[1].charge))
    s
end

"generate a chemical species and plot it"
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
        tospecies,
    ) |> esc
end

"""
    tospecies(s::Sym)

Maps the variable to a species.
We probably want to do this async somehow
"""
function tospecies(s; jsons=nothing)
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
                j, jview = PubChemReactions.get_json_and_view_from_cid(cid) #; verbose=true);
            else
                j, jview = PubChemReactions.get_json_and_view_from_cname(cname) #; verbose=true);
            end
        else
            j, jview = jsons
        end
        s = set_species_metadata(s, j, jview)
    end
end
tospecies(s::Num; kwargs...) = Num(tospecies(Symbolics.value(s); kwargs...))

"makes a Num with the Compound metadata"
function search_compound(cname)
    csym = Symbol(cname)
    csym = Symbolics.unwrap(first(@variables $csym(Catalyst.DEFAULT_IV)))
    j, jview = PubChemReactions.get_json_and_view_from_cname(cname) #; verbose=true);
    set_species_metadata(csym, j, jview)
end
