using Pkg
Pkg.activate(".")
try
    Pkg.add(url = "https://github.com/GenomicBreeding/GBCore.jl")
catch
    nothing
end
using GBPlots
using GBCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra
using DataFrames
using Plots, StatsPlots, Distances, Clustering, Measures
using GLMakie, ColorSchemes
