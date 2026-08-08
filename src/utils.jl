"""
    make_at_species(sps; current_name = true) -> String

Return a Catalyst `@species` expression string for PubChem species `sps`.

# Arguments
- `sps`: PubChem species from which to derive the species declaration.

# Keywords
- `current_name::Bool = true`: use the current Symbolics variable name when `true`;
  use the PubChem record title when `false`.

# Returns
- `String`: a Catalyst `@species` declaration including the compound identifier and
  local cache options.

# Examples
```julia
water = PubChemReactions.search_compound("water")
make_at_species(water)
```
"""
function make_at_species(sps; current_name = true)
    name = current_name ? SymbolicIndexingInterface.getname(sps) :
        """var"$(get_name(sps))\""""
    return """@species $(name)(t) [cid = $(get_cid(sps)), save = true, load = true]"""
end

function eq_str_to_wl(str)
    str = replace(str, "~" => "==")
    str = replace(str, "&" => "&&")
    return str = replace(str, "|" => "||")
end

"""
    eqs_to_mathematica(eqs) -> String

Convert symbolic equations to a Wolfram Language list string.

# Arguments
- `eqs`: iterable of symbolic equations whose textual form uses `~`, `&`, or `|`.

# Returns
- `String`: braced Wolfram Language expression list with equality and boolean
  operators translated.

# Examples
```julia
using Symbolics: @variables

@variables x y
eqs_to_mathematica([x ~ y])
```
"""
function eqs_to_mathematica(eqs)
    es = string.(eqs)
    eq_strs = map(eq_str_to_wl, es)
    return join(["{", join(eq_strs, ", "), "}"])
end

function rxns_to_wl(rxns)
end

function cid_from_substance_json(x)
    for s in x.Record.Section
        s.TOCHeading == "Related Records" &&
            return splitdir(s.Section[1].Information[1].Value.StringWithMarkup[1].Markup[1].URL)[end]
    end
    return nothing
end
