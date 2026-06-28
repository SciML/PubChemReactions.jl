using SciMLTesting, PubChemReactions, Test
using JET

run_qa(
    PubChemReactions;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    ei_kwargs = (;
        # `escapeuri` is owned by URIs (an external dep) but accessed via HTTP, which
        # re-exports it non-publicly; HTTP is the natural import here, so keep the ignore.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :escapeuri,   # owner URIs, accessed via HTTP
            ),
        ),
        # Qualified accesses to non-SciML deps' currently-non-public names; ignore until
        # those packages mark them public (or declare them with `public`). Verified still
        # non-public against the released graph (JSON3 1.14, Catalyst 16.2, HTTP 1.11,
        # CSV 0.10); none are SciMLBase/DiffEqBase-owned.
        all_qualified_accesses_are_public = (;
            ignore = (
                :Array,            # JSON3
                :DEFAULT_IV,       # Catalyst
                :Object,           # JSON3
                :VariableSpecies,  # Catalyst
                :escapeuri,        # HTTP
                :get,              # HTTP
                :read,             # JSON3 and CSV
                :write,            # JSON3
            ),
        ),
    ),
    # The module relies on many implicit imports from its heavy `using` deps (CSV,
    # Cascadia, Catalyst, Symbolics, SymbolicUtils, DataFrames, Downloads, FileIO,
    # Graphs, Gumbo, HTTP, ImageIO, JSON3, PeriodicTable, Plots, StatsBase). Making
    # them explicit is a large, risky refactor; tracked rather than mass-rewritten.
    # https://github.com/SciML/PubChemReactions.jl/issues/62
    ei_broken = (:no_implicit_imports,),
)
