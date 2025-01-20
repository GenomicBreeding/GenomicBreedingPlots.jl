module GBPlots

using GBCore
using StatsBase, Distributions, LinearAlgebra
using DataFrames
using Plots, StatsPlots, Distances, Clustering, Measures
using PlotlyJS, PlotlyBase, Genie, GenieFramework
using PrecompileTools: @compile_workload

include("commons.jl")
# include("genomes.jl")
include("phenomes.jl")
# include("trials.jl")
# include("tebv.jl")

export PlotsGB, DistributionPlots, ViolinPlots, CorHeatPlots, TreePlots
export checkdims, labeltofname, saveplots
export plotstatic, plotinteractive

@compile_workload begin
    true
end

end
