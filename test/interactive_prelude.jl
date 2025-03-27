using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using GenomicBreedingPlots
using GenomicBreedingCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra
using DataFrames, Random
using Distances, Clustering, Measures
using CairoMakie, ColorSchemes
