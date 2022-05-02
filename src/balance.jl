
"check that the element counts in substrates is equal to products"
function isbalanced(rxn)
    all(hasmetadata.(species(rxn), Compound)) || error("some species do not have atom graph metadata")
    all(hasmetadata.(species(rxn), AtomBondGraph)) || error("some species do not have atom graph metadata")
    atom_counts(rxn.substrates, rxn.substoich) == atom_counts(rxn.products, rxn.prodstoich)
end

"check that the element counts in sub"
function isbalanced(rn::ReactionSystem)
    all(isbalanced.(reactions(rn)))
end

replace_atom_counts_with_elements(atomcounts) = PeriodicTable.elements[first.(atomcounts)] .=> last.(atomcounts)
replace_atom_counts_with_elements(atomcounts::Dict) = Dict(PeriodicTable.elements[collect(keys(atomcounts))] .=> values(atomcounts))

function atom_counts(s::Num)
    c = getmetadata(s, AtomBondGraph)
    aps = c.atoms
    countmap(last.(aps))
end

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
function balance(substrates, products; k=nothing, verbose=true)
    all_species = vcat(substrates, products)
    all(PubChemReactions.isspecies.(all_species)) || error("provide chemcial species (with graphs)")

    occuring_elements = collect(get_elements(all_species))
    atom_counts = PubChemReactions.atom_counts.(all_species)
    charges = map(x -> get_charge.(x), (substrates, products))
    n_subs = length(substrates)
    n_prods = length(products)

    n_elems = length(occuring_elements)
    n_eqs = n_specs = length(all_species)

    A = zeros(Int, n_specs, n_specs)

    for i in 1:n_elems
        for j in 1:n_specs
            amt_of_i = occuring_elements[i]
            coeff = j > n_subs ? -1 : 1
            A[i, j] = haskey(atom_counts[j], amt_of_i) ? coeff * atom_counts[j][amt_of_i] : 0
        end
    end

    # check if we need another equation for charges
    if !all(charges .== 0)
        for j in 1:n_specs # big hack, needs to be cleaned up 
            coeff = j > n_subs ? -1 : 1
            A[n_elems+1, j] = coeff * get_charge(all_species[j])
        end
    end

    A[end] = 1 # extra so not underdetermined

    b = zeros(Int, n_specs)
    b[end] = 1 # extra equation 

    x = A \ b
    x .= x ./ minimum(x)
    x = round.(Int, x)

    # need a better way to set rate. wolfram doesn't include rate in the Reaction type, just subs, prods, and stoichs
    k = k === nothing ? 1 : k
    Reaction(k, substrates, products, x[1:n_subs], x[n_subs+1:end])
end
