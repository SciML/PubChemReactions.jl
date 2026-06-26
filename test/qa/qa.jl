using SciMLTesting, PubChemReactions, Test
using JET

run_qa(
    PubChemReactions;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    ei_kwargs = (;
        # Names re-exported by a non-owner module (they resolve, just via a re-export);
        # they go owner-correct as the base libs settle their re-export surface.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :escapeuri,   # owner URIs, accessed via HTTP
                :getname,     # owner SymbolicIndexingInterface, accessed via Symbolics
                :scalarize,   # owner SymbolicUtils, accessed via Symbolics
                :symtype,     # owner SymbolicUtils, accessed via Symbolics
                :unwrap,      # owner SymbolicUtils, accessed via Symbolics
            ),
        ),
        # Qualified accesses to other packages' currently-non-public names; ignore until
        # those packages mark them public (or declare them with `public`).
        all_qualified_accesses_are_public = (;
            ignore = (
                :Array,            # JSON3
                :DEFAULT_IV,       # Catalyst
                :Object,           # JSON3
                :VariableSpecies,  # Catalyst
                :escapeuri,        # HTTP
                :get,              # HTTP
                :getname,          # Symbolics (re-export of a SymbolicIndexingInterface internal)
                :read,             # JSON3 and CSV
                :symtype,          # Symbolics (re-export of a SymbolicUtils internal)
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
