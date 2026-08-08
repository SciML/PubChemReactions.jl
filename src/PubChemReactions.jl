module PubChemReactions

import Cascadia, Downloads, Gumbo, ImageIO, JSON, PeriodicTable, SymbolicIndexingInterface,
    SymbolicUtils, URIs
import Catalyst, Symbolics
using Catalyst: Reaction, ReactionSystem, reactions
using Cascadia: Selector
using DataFrames: DataFrame
using DelimitedFiles: readdlm
using FileIO: load
using Graphs: SimpleGraph, add_edge!
using Gumbo: parsehtml
import Plots
using StatsBase: countmap
using SymbolicUtils: getmetadata, hasmetadata, setmetadata
using Symbolics: Equation, Num, @variables

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
include("balance.jl")
include("pathway.jl")
include("utils.jl")

export Compound
# `@species` is intentionally not exported because it clashes with `Catalyst.@species`.
# Use `PubChemReactions.@species` or import it explicitly.
export @species_str
export get_cid, get_name, get_charge, get_graph
export atom_counts, element_counts, atom_matrix
export save_species, load_species, isspecies
export pubchem_search
export balance_eqs, isbalanced # , balance
export atomplot, atomplot2d, atomplot3d
export make_at_species, eqs_to_mathematica
export pathway_reaction

"""
    pc() -> Bool

Open the PubChem homepage in the system default browser.

# Returns
- `Bool`: whether the browser command was started successfully.

# Examples
```julia
pc()
```
"""
pc() = open_in_default_browser(PC_ROOT)
export pc

end # module
