using SciMLTesting, PubChemReactions, Test
using JET

run_qa(
    PubChemReactions;
    jet_kwargs = (; target_defined_modules = true),
    ei_kwargs = (;
        # `escapeuri` is owned by URIs but accessed via HTTP. `scalarize` and
        # `unwrap` are public Symbolics entry points whose methods are owned by
        # SymbolicUtils, so Base.which reports the implementation owner.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :escapeuri, # owner URIs, accessed via HTTP
                :scalarize, # owner SymbolicUtils, accessed via Symbolics
                :unwrap,    # owner SymbolicUtils, accessed via Symbolics
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
)
