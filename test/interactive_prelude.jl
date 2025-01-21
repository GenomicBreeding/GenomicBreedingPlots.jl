using Pkg
Pkg.activate(".")
Pkg.add(url = "https://github.com/GenomicBreeding/GBCore.jl")
using GBPlots
using GBCore
using StatsBase, Distributions, LinearAlgebra
using DataFrames
using Plots, StatsPlots, Distances, Clustering, Measures
