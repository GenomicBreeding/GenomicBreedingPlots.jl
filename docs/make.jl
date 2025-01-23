using Pkg
# Pkg.develop(PackageSpec(path=pwd()))
Pkg.activate(".")
using GBPlots
using Documenter

DocMeta.setdocmeta!(GBPlots, :DocTestSetup, :(using GBPlots); recursive = true)

makedocs(;
    modules = [GBPlots],
    authors = "jeffersonparil@gmail.com",
    sitename = "GBPlots.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GBPlots.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GBPlots.jl", devbranch = "main")
