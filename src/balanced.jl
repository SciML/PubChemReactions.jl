
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

function atom_countmap(s)
    c = getmetadata(s, AtomBondGraph)
    aps = c.atoms
    countmap(last.(aps))
end

function atom_counts(speciess, stoichs)
    countmaps = atom_countmap.(speciess)

    for (stoich, cm) in zip(stoichs, countmaps)
        for (k, v) in cm 
            cm[k] = stoich * v
        end
    end

    mergewith(+, countmaps...)
end

elements(s) = unique(last.(get_graph(s).atoms))
elements(s::Vector) = Set(reduce(vcat, elements.(s)))

"""


# http://mathgene.usc.es/matlab-profs-quimica/reacciones.pdf

should i try to catch underdetermined soon, or just let LA give SingularException?


"""
function get_balanced_reaction(substrates, products; verbose=true)
    all_species = vcat(substrates, products)
    all(PubChemReactions.isspecies.(all_species)) || error("provide chemcial species (with graphs)")
    
    occuring_elements = collect(PubChemReactions.elements(all_species))
    atom_counts = PubChemReactions.atom_countmap.(all_species)
    charges = map(x->get_charge.(x), (substrates, products))
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

    x = A\b
    x .= x ./ minimum(x)
    x = round.(Int , x)

    # need a better way to set rate. wolfram doesn't include rate in the Reaction type, just subs, prods, and stoichs
    Reaction(1, substrates, products, x[1:n_subs], x[n_subs+1:end])
end
