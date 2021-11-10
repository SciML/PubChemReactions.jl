function parse_rhea_equation(eq::AbstractString)
    reactants, products = PubChemReactions.rhea_to_reacts_prods(eq)
    rs = map(PubChemReactions.make_stoich_from_rhea, reactants)
    ps = map(PubChemReactions.make_stoich_from_rhea, products)
    rstoich, reactants = first.(rs), last.(rs)
    pstoich, products = first.(ps), last.(ps)

    search_compound.(reactants), search_compound.(products), rstoich, pstoich
end

"includes stoich values" 
function rhea_to_reacts_prods(eq::AbstractString)
    eq = foldl(replace, ARROWS .=> "=", init=eq)
    lhs, rhs = split(eq, " = ")
    split(lhs, " + "), split(rhs, " + ")
end

function make_stoich_from_rhea(s)
    if startswith(s, r"(\d).* ")
        ss = split(s, " ")
        parse(Int, ss[1]), ss[2]
    else 
        1, s
    end
end

"searches the Rhea reactions DB for reactions that include the species in args"
function get_biochem_rxns(csym, csyms...)
    chebi_ids = get_chebi_id(csym)
    for c in csyms
        chebi_id = get_chebi_id(c)
        chebi_ids = chebi_ids * "+" * chebi_id
    end

    input_url = "$RHEA_URL/?query=$(chebi_ids)&columns=rhea-id,equation,chebi-id&format=tsv"
    res = HTTP.get(input_url)
    if res.status == 200
        return CSV.read(IOBuffer(res.body), DataFrame)
    else
        error("Cannot find Biochemical reactions")
    end
end
