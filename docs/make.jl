# Standard stuff
cd(@__DIR__)
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using CairoMakie, Documenter, Literate
using DocumenterTools: Themes
using DocumenterCitations
ENV["JULIA_DEBUG"] = "Documenter"

# Packages specific to these docs
using Hydrology

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)

Literate.markdown("src/examples/Kazmierczak2024.jl", "src/examples"; credit = false)

example_pages = [
    "examples/Kazmierczak2024.md",
]

ref_pages = ["references.md"]
# %% Build docs
PAGES = [
    "index.md",
    "Examples" => example_pages,
    "References" => ref_pages,
]

include("style.jl")

makedocs(
    modules = [Hydrology],
    format = Documenter.HTML(
        prettyurls = CI,
        assets = [
            asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
        collapselevel = 2,
        ),
    sitename = "Hydrology.jl",
    authors = "Takis Angelides",
    pages = PAGES,
    doctest = CI,
    draft = false,
    plugins = [bib],
    checkdocs = :none,
    warnonly = true,
)

deploydocs(
    repo = "github.com/TakisAngelides/Hydrology.jl",
    devbranch = "main"
)
