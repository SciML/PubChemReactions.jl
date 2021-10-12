using PubChemReactions, Catalyst, Test
using OrdinaryDiffEq
using MetaGraphs, LightGraphs, NetworkLayout
using Downloads, FileIO, Images
using PhysicalConstants.CODATA2018
using Unitful
using Symbolics

# note biostructures.jl has a lot of this https://github.com/BioJulia/BioStructures.jl/blob/master/src/model.jl#L102
abstract type AbstractAtom end 

struct Compound{T}
    "pubchem metadata"
    num_metadata # but we dont have the actual Nums? i guess thats fine
    node_atom_pairs::Vector{Pair{T, Atom}}
    g::MultiGraph{T}
end

struct Atom <: AbstractAtom 
    protons
    neutrons
    electrons
end

proton() = Atom(1, 0, 0)
neutron() = Atom(0, 1, 0)
electron() = Atom(0, 0, 1)
# Base.:+(a1::Atom, a2::Atom) = Atom(a1.protons + a2.protons, a1.neutrons + a2.neutrons, a1.electrons + a2.electrons)
p = proton()

# mass(a) = 

# "modifies both "
# function covalent_bond!(g, a1, a2) end


# valence(a::Atom)

# valence(c::Compound, a::Atom)

charge(a::Atom) = (a.protons - a.electrons) * ElementaryCharge
charge_num(a::Atom) = (a.protons - a.electrons)
# charge(g, a::Atom) = a.protons - a.electrons + sum(mul.(edges(g, e))

