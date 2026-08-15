module PubChemReactions

import CSV, Cascadia, Catalyst, Downloads, Graphs, Gumbo, HTTP, ImageIO, JSON3,
    PeriodicTable, Plots, SymbolicIndexingInterface, Symbolics, SymbolicUtils
using Cascadia: Selector
using Catalyst: Reaction, ReactionSystem, reactions
using DataFrames: DataFrame
using FileIO: load
using Graphs: SimpleGraph, add_edge!
using Gumbo: parsehtml
using Plots: title!
using StatsBase: countmap
using Symbolics: Equation, Num, @variables
using SymbolicUtils: getmetadata, hasmetadata, operation, setmetadata
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
include("balance.jl")
include("pathway.jl")
include("utils.jl")

export Compound
# `@species` is intentionally not exported: it would clash with `Catalyst.@species`
# when both packages are `using`-imported. Call `PubChemReactions.@species` to use the
# PubChem-aware version (resolves `cid`/`name`/`save`/`load` metadata); `Catalyst.@species`
# remains available unqualified via `using Catalyst`.
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
    pc()

Open the PubChem homepage in the system default browser.

# Returns

The `Bool` result from the platform browser launcher.

# Examples

```julia
pc()
```
"""
pc() = open_in_default_browser(PC_ROOT)
export pc

end # module
