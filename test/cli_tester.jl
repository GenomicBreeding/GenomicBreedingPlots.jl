using Pkg
Pkg.activate(".")
using JuliaFormatter
Pkg.activate(".")
# Format
format(".")
# Test
include("runtests.jl")
