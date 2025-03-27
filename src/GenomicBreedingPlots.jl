module GenomicBreedingPlots

using GenomicBreedingCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra
using DataFrames, Random
using Distances, Clustering, Measures
using CairoMakie, ColorSchemes

include("commons.jl")
include("genomes.jl")
include("phenomes.jl")
include("cv.jl")
# include("trials.jl")
# include("tebv.jl")

export PlotsGB, DistributionPlots, ViolinPlots, CorHeatPlots, TreePlots, BarPlots, BoxPlots, PCBiPlots
export checkdims, labeltofname, saveplots
export plot

end
