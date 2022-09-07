module PubChemReactions

using JSON3, HTTP, Symbolics, CSV, DataFrames
using Catalyst, Graphs
using StatsBase
using ImageIO, FileIO, Plots
using Downloads
# using URIs maybe
using PeriodicTable
using Cascadia, Gumbo # grumble grumble
# using Groebner, DynamicPolynomials
# using PolynomialRoots, PolynomialFactors

const PC_ROOT = "https://pubchem.ncbi.nlm.nih.gov"
const PUG_URL = joinpath(PC_ROOT, "rest/pug")
const PUG_VIEW_URL = joinpath(PC_ROOT, "rest/pug_view")
const RXN_TABLE_BASE_URL = joinpath(PC_ROOT, "sdq/sdqagent.cgi")

const RHEA_URL = "https://www.rhea-db.org/rhea"

const ARROWS = ["<->", "->", "<-", "<=>", "=>", "⟶"] # lame hack for parsing string reactions

const DATADIR = joinpath(@__DIR__, "../data/")
const COMPOUNDS_DIR = joinpath(DATADIR, "compounds")

include("data.jl")
include("rhea.jl")
include("graph.jl")
include("plot.jl")
include("species.jl")
# include("balance.jl")
include("pathway.jl")
include("utils.jl")

export Compound
export @species_str, @species
export get_cid, get_name, get_charge, get_graph
export atom_counts, element_counts, atom_matrix
export save_species, load_species, isspecies
export pubchem_search
# export balance, balance_eqs, isbalanced
export atomplot, atomplot2d, atomplot3d
export make_at_species, eqs_to_mathematica
export pathway_reaction

pc() = open_in_default_browser(PC_ROOT)
export pc

end # module
