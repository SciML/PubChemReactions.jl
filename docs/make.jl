using Documenter
using PubChemReactions

makedocs(;
    modules = [PubChemReactions],
    sitename = "PubChemReactions.jl",
    pages = [
        "Public API" => "index.md",
    ],
    doctest = true,
    checkdocs = :exports,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
)
