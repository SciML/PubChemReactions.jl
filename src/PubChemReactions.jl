module PubChemReactions

using JSON3, HTTP, Symbolics, CSV, DataFrames
using Catalyst, Graphs
using StatsBase
using Images, FileIO, Plots
using Downloads

const PUG_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
const PUG_VIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
const RHEA_URL = "https://www.rhea-db.org/rhea"
const ARROWS = ["<->", "->", "<-", "<=>", "=>"] # lame hack for parsing string reactions

include("data.jl")
include("rhea.jl")
include("graph.jl")
include("plot.jl")
include("balanced.jl")

export Compound
export @species_str

end # module
