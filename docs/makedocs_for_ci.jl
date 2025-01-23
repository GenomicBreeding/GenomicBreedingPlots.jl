using Pkg
cd("..")
# Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()
include("make.jl")