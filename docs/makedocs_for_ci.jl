using Pkg
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()
include("make.jl")