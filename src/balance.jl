"also probably bad idea maybe"
Base.diff(rxn::Catalyst.Reaction) = replace_atom_counts_with_elements(mergewith(-, reverse(atom_counts(rxn))...))

function isbalanced(substrates, products; substoich=ones(length(substrates)), prodstoich=ones(length(products)))
    atom_counts(substrates, substoich) == atom_counts(products, prodstoich)
end

"check that the element counts in substrates is equal to products"
function isbalanced(rxn)
    all(hasmetadata.(species(rxn), Compound)) || error("some species do not have atom graph metadata")
    all(hasmetadata.(species(rxn), AtomBondGraph)) || error("some species do not have atom graph metadata")
    isbalanced(rxn.substrates, rxn.products; substoich=rxn.substoich, prodstoich=rxn.prodstoich)
end

"check that the element counts in sub"
function isbalanced(rn::ReactionSystem)
    all(isbalanced.(reactions(rn)))
end

replace_atom_counts_with_elements(atomcounts) = PeriodicTable.elements[first.(atomcounts)] .=> last.(atomcounts)
replace_atom_counts_with_elements(atomcounts::Dict) = Dict(PeriodicTable.elements[collect(keys(atomcounts))] .=> values(atomcounts))
element_counts(x) = replace_atom_counts_with_elements(atom_counts(x))
element_counts(x::Reaction) = replace_atom_counts_with_elements.(atom_counts(x))

function atom_counts(s::Num)
    c = getmetadata(s, AtomBondGraph)
    aps = c.atoms
    countmap(last.(aps))
end
atom_counts(xs::Vector{Num}) = mergewith(+, atom_counts.(xs)...)
atom_counts(s::S) where {S<:SymbolicUtils.Symbolic} = atom_counts(Num(s))
atom_counts(rxn::Reaction) = (atom_counts(rxn.substrates, rxn.substoich), atom_counts(rxn.products, rxn.prodstoich))

function atom_counts(speciess, stoichs)
    countmaps = atom_counts.(speciess)

    for (stoich, cm) in zip(stoichs, countmaps)
        for (k, v) in cm
            cm[k] = stoich * v
        end
    end

    mergewith(+, countmaps...)
end

get_elements(s) = unique(last.(get_graph(s).atoms))
get_elements(s::Vector) = Set(reduce(vcat, get_elements.(s)))

balance(rxn; kwargs...) = balance(rxn.substrates, rxn.products; k=rxn.rate, kwargs...)

"""


# http://mathgene.usc.es/matlab-profs-quimica/reacciones.pdf

should i try to catch underdetermined soon, or just let LA give SingularException?


"""
function balance(substrates, products; k=nothing, add_constraint_eq=true, force_integer_stoich=true, short_circuit=true, verbose=true)
    # hack for now 
    k = k === nothing ? 1 : k
    # might want an early exit
    short_circuit && isbalanced(substrates, products) && return Reaction(k, substrates, products)

    all_species = vcat(substrates, products)
    all(PubChemReactions.isspecies.(all_species)) || error("provide chemcial species (with graphs)")

    occuring_elements, atomcounts, chgs, n_specs, n_subs = balance_setup(substrates, products)

    @polyvar x[1:n_specs] # couldn't get to work with @variables 
    eqs = balance_eqs(x, occuring_elements, atomcounts, chgs, n_specs, n_subs; add_constraint_eq)
    ts = eq_to_term.(eqs)
    newt = groebner(ts)
    sol = only(realsolutions(Symbolics.generic_extension_solve(newt)))

    if force_integer_stoich
        sol ./= minimum(real.(sol))
        sol = 1000 * (sol) # for the love of god I want a symbolic solver that handles infinite solutions
        sol = convert.(Int, sol)
        sol ./= gcd(sol)
    end
    @assert all(>(0), sol)
    rxn = Reaction(k, substrates, products, sol[1:n_subs], sol[n_subs+1:end])
    @info rxn
    @assert isbalanced(rxn)
    rxn
end

function balance_eqs(x, occuring_elements, atomcounts, chgs, n_specs, n_subs; add_constraint_eq=false)
    eqs = Equation[]
    for (i, e) in enumerate(occuring_elements)
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
    eqs
end

"refactor"
function balance_setup(substrates, products)
    all_species = vcat(substrates, products)
    all(PubChemReactions.isspecies.(all_species)) || error("provide chemcial species (with graphs)")

    occuring_elements = collect(PubChemReactions.get_elements(all_species))
    atomcounts = PubChemReactions.atom_counts.(all_species)
    charges = map(x -> PubChemReactions.get_charge.(x), (substrates, products))
    chgs = reduce(vcat, charges)

    n_subs = length(substrates)
    n_elems = length(occuring_elements)
    has_charges = !all(charges .== 0)
    n_specs = length(all_species)

    occuring_elements, atomcounts, chgs, n_specs, n_subs
end

function balance_eqs(substrates, products; add_constraint_eq=true)
    occuring_elements, atomcounts, chgs, n_specs, n_subs = balance_setup(substrates, products)
    @variables x[1:n_specs]
    balance_eqs(x, occuring_elements, atomcounts, chgs, n_specs, n_subs; add_constraint_eq)
end

balance_eqs(rxn::Reaction; add_constraint_eq=true) = balance_eqs(rxn.substrates, rxn.products; add_constraint_eq)
eq_to_term(eq) = eq.lhs - eq.rhs

function atom_matrix(rxn::Reaction)
    eqs = balance_eqs(rxn)
    ts = eq_to_term.(eqs)
    # vars = unique(reduce(vcat, Symbolics.get_variables.(ts))) # get_variables permutes, since the order they show in eqs
    # vars = unique(reduce(vcat, x[1:4]))
    vars = only(@variables x[1:length(species(rxn))])
    vars = collect(vars)
    a, b, islinear = Symbolics.linear_expansion(ts, vars)
    a
end
