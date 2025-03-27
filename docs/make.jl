using GenomicBreedingPlots
using Documenter

DocMeta.setdocmeta!(GenomicBreedingPlots, :DocTestSetup, :(using GenomicBreedingPlots); recursive = true)

makedocs(;
    modules = [GenomicBreedingPlots],
    authors = "jeffersonparil@gmail.com",
    sitename = "GenomicBreedingPlots.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GenomicBreedingPlots.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GenomicBreedingPlots.jl", devbranch = "main")
