function atom_count_diff(rxn::Catalyst.Reaction)
    return replace_atom_counts_with_elements(mergewith(-, reverse(atom_counts(rxn))...))
end

"""
    isbalanced(substrates, products; substoich = ones(length(substrates)),
        prodstoich = ones(length(products)))

Determine whether the substrate and product collections contain the same
stoichiometric atom counts.

# Arguments

- `substrates`: Collection of PubChem species on the left-hand side. Every
  species must have the atom-bond graph metadata used by `PubChemReactions`.
- `products`: Collection of PubChem species on the right-hand side.

# Keyword Arguments

- `substoich`: Numeric coefficients for `substrates`. Its length must match
  `substrates`.
- `prodstoich`: Numeric coefficients for `products`. Its length must match
  `products`.

# Returns

`true` when the weighted atom counts match exactly, and `false` otherwise.

# Examples

```julia
using PubChemReactions

h2 = PubChemReactions.search_compound("hydrogen")
o2 = PubChemReactions.search_compound("oxygen")
h2o = PubChemReactions.search_compound("water")

isbalanced([h2, o2], [h2o]; substoich = [2, 1], prodstoich = [2])
```
"""
function isbalanced(
        substrates, products; substoich = ones(length(substrates)),
        prodstoich = ones(length(products))
    )
    return atom_counts(substrates, substoich) == atom_counts(products, prodstoich)
end

"""
    isbalanced(rxn)

Determine whether the substrate and product species in Catalyst reaction `rxn`
have matching stoichiometric atom counts.

# Arguments

- `rxn`: A Catalyst `Reaction` whose species have PubChem `Compound` and
  atom-bond graph metadata.

# Returns

`true` when the reaction is element-balanced, and `false` otherwise. An error
is thrown when a reaction species lacks the metadata required to count atoms.

The reaction's `substoich` and `prodstoich` fields are passed to the collection
method, so the check uses the reaction's current stoichiometric coefficients.
"""
function isbalanced(rxn)
    all(hasmetadata.(reaction_species(rxn), Compound)) ||
        error("some species do not have atom graph metadata")
    all(hasmetadata.(reaction_species(rxn), AtomBondGraph)) ||
        error("some species do not have atom graph metadata")
    return isbalanced(rxn.substrates, rxn.products; substoich = rxn.substoich, prodstoich = rxn.prodstoich)
end

"""
    isbalanced(rn::ReactionSystem)

Determine whether every reaction in Catalyst reaction system `rn` is
element-balanced.

# Arguments

- `rn`: A Catalyst `ReactionSystem` whose reactions contain PubChem species.

# Returns

`true` when `isbalanced` returns `true` for every reaction in `rn`, and `false`
otherwise. Metadata errors from an individual reaction are propagated.
"""
function isbalanced(rn::ReactionSystem)
    return all(isbalanced.(reactions(rn)))
end

function replace_atom_counts_with_elements(atomcounts)
    return PeriodicTable.elements[first.(atomcounts)] .=> last.(atomcounts)
end
function replace_atom_counts_with_elements(atomcounts::Dict)
    return Dict(PeriodicTable.elements[collect(keys(atomcounts))] .=> values(atomcounts))
end

"""
    element_counts(x)

Return atom counts for `x` keyed by `PeriodicTable.Element` values.

For a species or species collection, the result is a dictionary keyed by
periodic-table elements. For a `Reaction`, the result is a tuple containing
the substrate and product dictionaries.

# Arguments

- `x`: PubChem species, a collection of species, or a Catalyst `Reaction`.

# Returns

A dictionary, or a substrate/product tuple of dictionaries, containing the
weighted atom counts keyed by `PeriodicTable.Element`.
"""
element_counts(x) = replace_atom_counts_with_elements(atom_counts(x))
element_counts(x::Reaction) = replace_atom_counts_with_elements.(atom_counts(x))

"""
    atom_counts(x)

Return atom counts for a PubChem species, reaction, or species collection.

# Arguments

- `x`: A PubChem species, a collection of species, or a Catalyst `Reaction`.

# Returns

For a species or collection, a dictionary keyed by the atomic numbers stored in
the PubChem atom-bond graph. For a `Reaction`, a tuple containing separate
substrate and product dictionaries is returned. Reaction dictionaries include
the reaction's stoichiometric coefficients.
"""
function atom_counts(s::Num)
    c = getmetadata(s, AtomBondGraph)
    aps = c.atoms
    return countmap(last.(aps))
end
atom_counts(xs::Vector{Num}) = mergewith(+, atom_counts.(xs)...)
atom_counts(s::S) where {S <: SymbolicUtils.BasicSymbolic} = atom_counts(Num(s))
function atom_counts(rxn::Reaction)
    return (atom_counts(rxn.substrates, rxn.substoich), atom_counts(rxn.products, rxn.prodstoich))
end

function atom_counts(speciess, stoichs)
    countmaps = atom_counts.(speciess)

    for (stoich, cm) in zip(stoichs, countmaps)
        for (k, v) in cm
            cm[k] = stoich * v
        end
    end

    return mergewith(+, countmaps...)
end

get_elements(s) = unique(last.(get_graph(s).atoms))
get_elements(s::Vector) = Set(reduce(vcat, get_elements.(s)))

# balance(rxn; kwargs...) = balance(rxn.substrates, rxn.products; k=rxn.rate, kwargs...)

# """

# # http://mathgene.usc.es/matlab-profs-quimica/reacciones.pdf

# should i try to catch underdetermined soon, or just let LA give SingularException?

# """
# function balance(substrates, products; k=nothing, add_constraint_eq=true, force_integer_stoich=true, short_circuit=true, verbose=true)
#     # hack for now
#     k = k === nothing ? 1 : k
#     # might want an early exit
#     short_circuit && isbalanced(substrates, products) && return Reaction(k, substrates, products)

#     all_species = vcat(substrates, products)
#     all(PubChemReactions.isspecies.(all_species)) || error("provide chemical species (with graphs)")

#     occurring_elements, atomcounts, chgs, n_specs, n_subs = balance_setup(substrates, products)

#     @polyvar x[1:n_specs] # couldn't get to work with @variables
#     eqs = balance_eqs(x, occurring_elements, atomcounts, chgs, n_specs, n_subs; add_constraint_eq)
#     ts = eq_to_term.(eqs)
#     newt = groebner(ts)
#     sol = only(realsolutions(Symbolics.generic_extension_solve(newt)))

#     if force_integer_stoich
#         sol ./= minimum(real.(sol))
#         sol = 1000 * (sol) # for the love of god I want a symbolic solver that handles infinite solutions
#         sol = convert.(Int, sol)
#         sol ./= gcd(sol)
#     end
#     @assert all(>(0), sol)
#     rxn = Reaction(k, substrates, products, sol[1:n_subs], sol[n_subs+1:end])
#     @info rxn
#     @assert isbalanced(rxn)
#     rxn
# end

function balance_eqs(
        x, occurring_elements, atomcounts, chgs, n_specs, n_subs; add_constraint_eq = false
    )
    eqs = Equation[]
    for (i, e) in enumerate(occurring_elements)
        lhs = 0
        rhs = 0
        for (j, d) in enumerate(atomcounts)
            !haskey(d, e) && continue
            term = d[e] * x[j]
            j <= n_subs ? lhs += term : rhs += term
        end
        eq = lhs ~ rhs
        push!(eqs, eq)
    end
    if !all(chgs .== 0)
        lhs = 0
        rhs = 0
        for (j, d) in enumerate(atomcounts)
            term = chgs[j] * x[j]
            j <= n_subs ? lhs += term : rhs += term
        end
        push!(eqs, lhs ~ rhs)
    end
    add_constraint_eq && push!(eqs, x[end] ~ 1)
    return eqs
end

"""
refactor
"""
function balance_setup(substrates, products)
    all_species = vcat(substrates, products)
    all(PubChemReactions.isspecies.(all_species)) ||
        error("provide chemical species (with graphs)")

    occurring_elements = collect(PubChemReactions.get_elements(all_species))
    atomcounts = PubChemReactions.atom_counts.(all_species)
    charges = map(x -> PubChemReactions.get_charge.(x), (substrates, products))
    chgs = reduce(vcat, charges)

    n_subs = length(substrates)
    n_elems = length(occurring_elements)
    has_charges = !all(charges .== 0)
    n_specs = length(all_species)

    return occurring_elements, atomcounts, chgs, n_specs, n_subs
end

"""
    balance_eqs(substrates, products; add_constraint_eq = true)

Construct symbolic element and charge balance equations for collections of
PubChem substrate and product species.

# Arguments

- `substrates`: Collection of PubChem species on the left-hand side.
- `products`: Collection of PubChem species on the right-hand side.

# Keyword Arguments

- `add_constraint_eq`: If `true`, append an equation fixing the coefficient of
  the final species to `1`. This removes the scale ambiguity from the returned
  homogeneous balance equations. The default is `true`.

# Returns

A `Vector{Symbolics.Equation}` containing one symbolic coefficient for each
species in `vcat(substrates, products)`, with equations for every occurring
element and, when needed, total charge. The substrate coefficients appear on
the left-hand side and product coefficients on the right-hand side.

# Examples

```julia
using PubChemReactions

h2 = PubChemReactions.search_compound("hydrogen")
o2 = PubChemReactions.search_compound("oxygen")
h2o = PubChemReactions.search_compound("water")
eqs = balance_eqs([h2, o2], [h2o])
```
"""
function balance_eqs(substrates, products; add_constraint_eq = true)
    occurring_elements, atomcounts, chgs, n_specs,
        n_subs = balance_setup(substrates, products)
    @variables x[1:n_specs]
    return balance_eqs(x, occurring_elements, atomcounts, chgs, n_specs, n_subs; add_constraint_eq)
end

"""
    balance_eqs(rxn::Reaction; add_constraint_eq = true)

Construct symbolic element and charge balance equations for the substrate and
product species in Catalyst reaction `rxn`.

# Arguments

- `rxn`: A Catalyst `Reaction` containing PubChem species. Its current
  stoichiometric coefficients are not used; the returned equations introduce a
  new coefficient for each species.

# Keyword Arguments

- `add_constraint_eq`: If `true`, append an equation fixing the coefficient of
  the final species to `1`. The default is `true`.

# Returns

A `Vector{Symbolics.Equation}` equivalent to
`balance_eqs(rxn.substrates, rxn.products; add_constraint_eq)`.
"""
function balance_eqs(rxn::Reaction; add_constraint_eq = true)
    return balance_eqs(rxn.substrates, rxn.products; add_constraint_eq)
end
eq_to_term(eq) = eq.lhs - eq.rhs

"""
    atom_matrix(rxn::Reaction)

Return the linear coefficient matrix for the balance equations of `rxn`.

# Arguments

- `rxn`: A Catalyst `Reaction` containing PubChem species.

# Returns

A matrix whose columns correspond to the species in the reaction and whose
rows contain the coefficients of the element and charge balance equations. The
matrix is formed from `balance_eqs(rxn)` and therefore includes the
normalization equation by default.
"""
function atom_matrix(rxn::Reaction)
    eqs = balance_eqs(rxn)
    ts = eq_to_term.(eqs)
    # vars = unique(reduce(vcat, Symbolics.get_variables.(ts))) # get_variables permutes, since the order they show in eqs
    # vars = unique(reduce(vcat, x[1:4]))
    vars = only(@variables x[1:length(reaction_species(rxn))])
    vars = collect(vars)
    a, b, islinear = Symbolics.linear_expansion(ts, vars)
    return a
end
