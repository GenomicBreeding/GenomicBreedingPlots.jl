module GBPlots

using GBCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra
using DataFrames
using Distances, Clustering, Measures
using CairoMakie, ColorSchemes
# using PrecompileTools: @compile_workload

include("commons.jl")
# include("genomes.jl")
include("phenomes.jl")
# include("trials.jl")
# include("tebv.jl")
include("cv.jl")

export PlotsGB, DistributionPlots, ViolinPlots, CorHeatPlots, TreePlots, BarPlots, BoxPlots
export checkdims, labeltofname, saveplots
export plot

# @compile_workload begin
#     true
# end

end
