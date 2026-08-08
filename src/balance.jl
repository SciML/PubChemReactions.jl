function atom_count_diff(rxn::Catalyst.Reaction)
    return replace_atom_counts_with_elements(mergewith(-, reverse(atom_counts(rxn))...))
end

"""
    isbalanced(
            substrates, products; substoich = ones(length(substrates)),
            prodstoich = ones(length(products))
        ) -> Bool
    isbalanced(rxn::Catalyst.Reaction) -> Bool
    isbalanced(rn::Catalyst.ReactionSystem) -> Bool

Return whether PubChem species have equal element counts on both sides of a
reaction.

# Arguments
- `substrates`: iterable of PubChem species on the reactant side.
- `products`: iterable of PubChem species on the product side.
- `rxn::Catalyst.Reaction`: reaction to check.
- `rn::Catalyst.ReactionSystem`: every reaction in this system is checked.

# Keywords
- `substoich = ones(length(substrates))`: stoichiometric coefficients for
  `substrates`.
- `prodstoich = ones(length(products))`: stoichiometric coefficients for
  `products`.

# Returns
- `Bool`: `true` only when both sides have identical element-count mappings.

# Throws
- `ErrorException`: if reaction species do not carry PubChem structural metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
isbalanced([water], [water])
```
"""
function isbalanced(
        substrates, products; substoich = ones(length(substrates)),
        prodstoich = ones(length(products))
    )
    return atom_counts(substrates, substoich) == atom_counts(products, prodstoich)
end

function isbalanced(rxn)
    all(hasmetadata.(reaction_species(rxn), Compound)) ||
        error("some species do not have atom graph metadata")
    all(hasmetadata.(reaction_species(rxn), AtomBondGraph)) ||
        error("some species do not have atom graph metadata")
    return isbalanced(rxn.substrates, rxn.products; substoich = rxn.substoich, prodstoich = rxn.prodstoich)
end

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

# Arguments
- `x`: PubChem species, reaction, or supported collection of PubChem species.

# Returns
- A mapping from `PeriodicTable.Element` values to atom counts, or a pair of
  mappings for a reaction.

# Examples
```julia
water = PubChemReactions.search_compound("water")
element_counts(water)
```
"""
element_counts(x) = replace_atom_counts_with_elements(atom_counts(x))
element_counts(x::Reaction) = replace_atom_counts_with_elements.(atom_counts(x))

"""
    atom_counts(x)

Return atom counts for a PubChem species, reaction, or stoichiometric species
collection.

# Arguments
- `x`: PubChem species, `Catalyst.Reaction`, or a vector of PubChem species.

# Returns
- `Dict`: element-number-to-count mapping for a species or collection.
- `Tuple{Dict, Dict}`: substrate and product mappings for a reaction.

# Throws
- `KeyError`: if a species has no atom-bond graph metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
atom_counts(water)
```
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

function _balance_eqs(
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
    balance_eqs(rxn::Catalyst.Reaction; add_constraint_eq = true)

Construct symbolic element and charge balance equations for PubChem species.

# Arguments
- `substrates`: iterable of PubChem species on the left side of a reaction.
- `products`: iterable of PubChem species on the right side of a reaction.
- `rxn::Catalyst.Reaction`: reaction from which substrates and products are read.

# Keywords
- `add_constraint_eq::Bool = true`: append an equation that fixes the final stoichiometric
  variable to one and removes the scale indeterminacy.

# Returns
- `Vector{Symbolics.Equation}`: homogeneous element and charge balance equations.

# Throws
- `ErrorException`: if any species has no PubChem structural metadata.

# Examples
```julia
water = PubChemReactions.search_compound("water")
balance_eqs([water], [water])
```
"""
function balance_eqs(substrates, products; add_constraint_eq = true)
    occurring_elements, atomcounts, chgs, n_specs,
        n_subs = balance_setup(substrates, products)
    @variables x[1:n_specs]
    return _balance_eqs(x, occurring_elements, atomcounts, chgs, n_specs, n_subs; add_constraint_eq)
end

function balance_eqs(rxn::Reaction; add_constraint_eq = true)
    return balance_eqs(rxn.substrates, rxn.products; add_constraint_eq)
end
eq_to_term(eq) = eq.lhs - eq.rhs

"""
    atom_matrix(rxn::Reaction) -> AbstractMatrix

Return the linear coefficient matrix for the balance equations of `rxn`.

# Arguments
- `rxn::Catalyst.Reaction`: reaction with PubChem species metadata.

# Returns
- `AbstractMatrix`: coefficient matrix for the equations produced by
  [`balance_eqs`](@ref).

# Throws
- `ErrorException`: if any species in `rxn` has no PubChem structural metadata.

# Examples
```julia
using Catalyst: Reaction

water = PubChemReactions.search_compound("water")
rxn = Reaction(1, [water], [water])
atom_matrix(rxn)
```
"""
function atom_matrix(rxn::Reaction)
    eqs = balance_eqs(rxn)
    ts = eq_to_term.(eqs)
    # vars = unique(reduce(vcat, Symbolics.get_variables.(ts))) # get_variables permutes, since the order they show in eqs
    # vars = unique(reduce(vcat, x[1:4]))
    n_species = length(rxn.substrates) + length(rxn.products)
    vars = only(@variables x[1:n_species])
    vars = collect(vars)
    a, b, islinear = Symbolics.linear_expansion(ts, vars)
    return a
end
