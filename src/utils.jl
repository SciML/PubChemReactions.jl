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
