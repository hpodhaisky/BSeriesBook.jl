# Build the documentation:
#     julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#     julia --project=docs docs/make.jl
# The result is in docs/build/index.html.

using Documenter
using BSeriesBook

DocMeta.setdocmeta!(BSeriesBook, :DocTestSetup, :(using BSeriesBook);
                    recursive = true)

makedocs(;
    modules = [BSeriesBook],
    sitename = "BSeriesBook.jl",
    authors = "John C. Butcher, Helmut Podhaisky",
    repo = Remotes.GitHub("hpodhaisky", "BSeriesBook.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://hpodhaisky.github.io/BSeriesBook.jl",
        edit_link = "main",
        assets = String[],
        example_size_threshold = 32_000,   # the tree table on trees.md
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Standard numbering" => "trees.md",
        "Conventions" => "conventions.md",
        "Underlying one-step method" => "glm.md",
        "Benchmarks" => "benchmarks.md",
        "Reference" => "reference.md",
        "Literature" => "literature.md",
    ],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/hpodhaisky/BSeriesBook.jl",
    devbranch = "main",
    push_preview = true,
)
