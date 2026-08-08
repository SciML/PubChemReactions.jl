"""
    parse_rhea_equation(eq::AbstractString)

Parse a Rhea equation into PubChem species and their stoichiometric coefficients.

# Arguments
- `eq::AbstractString`: Rhea equation containing one of the supported reaction arrows.

# Returns
- `Tuple`: reactant species, product species, reactant coefficients, and product
  coefficients, in that order.

# Throws
- `ErrorException`: if a named compound cannot be retrieved from PubChem.

# Examples
```julia
reactants, products, rstoich, pstoich =
    PubChemReactions.parse_rhea_equation("water = water")
```
"""
function parse_rhea_equation(eq::AbstractString)
    reactants, products = PubChemReactions.rhea_to_reacts_prods(eq)
    rs = map(PubChemReactions.make_stoich_from_rhea, reactants)
    ps = map(PubChemReactions.make_stoich_from_rhea, products)
    rstoich, reactants = first.(rs), last.(rs)
    pstoich, products = first.(ps), last.(ps)

    return search_compound.(reactants), search_compound.(products), rstoich, pstoich
end

function rhea_to_reacts_prods(eq::AbstractString)
    eq = foldl(replace, ARROWS .=> "=", init = eq)
    lhs, rhs = split(eq, " = ")
    return strip.(split(lhs, " + ")), strip.(split(rhs, " + "))
end

function make_stoich_from_rhea(s)
    return if startswith(s, r"(\d).* ")
        ss = split(s, " ")
        parse(Int, ss[1]), ss[2]
    else
        1, s
    end
end

"""
    get_biochem_rxns(csym, csyms...) -> DataFrame

Search Rhea for biochemical reactions that include all supplied PubChem species.

# Arguments
- `csym`: PubChem species used to start the Rhea query.
- `csyms...`: additional PubChem species required in the matching reactions.

# Returns
- `DataFrame`: matching Rhea identifiers, equations, and ChEBI identifiers.

# Throws
- `ErrorException`: if a species has no usable ChEBI reference or Rhea cannot be reached.

# Examples
```julia
water = PubChemReactions.search_compound("water")
reactions = PubChemReactions.get_biochem_rxns(water)
```
"""
function get_biochem_rxns(csym, csyms...)
    chebi_ids = get_chebi_id(csym)
    for c in csyms
        chebi_id = get_chebi_id(c)
        chebi_ids = chebi_ids * "+" * chebi_id
    end

    input_url = "$RHEA_URL/?query=$(chebi_ids)&columns=rhea-id,equation,chebi-id&format=tsv"
    data, header = readdlm(IOBuffer(get_page(input_url)), '\t', header = true)
    return DataFrame(data, vec(Symbol.(header)))
end
