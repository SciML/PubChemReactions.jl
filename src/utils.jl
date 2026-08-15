"""
    make_at_species(sps; current_name = true)

Return a Catalyst `@species` expression string for PubChem species `sps`.

# Arguments

- `sps`: A PubChem symbolic species with compound name and CID metadata.

# Keyword Arguments

- `current_name`: If `true`, use the current symbolic operation name. If
  `false`, use the PubChem compound name. The default is `true`.

# Returns

A string containing a Catalyst `@species` expression with the species CID and
`save`/`load` options enabled.
"""
function make_at_species(sps; current_name = true)
    name = current_name ? nameof(operation(Symbolics.value(sps))) :
        """var"$(get_name(sps))\""""
    return """@species $(name)(t) [cid = $(get_cid(sps)), save = true, load = true]"""
end

function eq_str_to_wl(str)
    str = replace(str, "~" => "==")
    str = replace(str, "&" => "&&")
    return str = replace(str, "|" => "||")
end

"""
    eqs_to_mathematica(eqs)

Convert symbolic equations to a Wolfram Language list string.

# Arguments

- `eqs`: An iterable collection of symbolic equations.

# Returns

A string containing the equations in a Wolfram Language list, with Julia's
`~`, `&`, and `|` operators converted to `==`, `&&`, and `||`.
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
