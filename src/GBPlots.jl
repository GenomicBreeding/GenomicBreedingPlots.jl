module GBPlots

using GBCore
using StatsBase, Distributions, LinearAlgebra
using DataFrames
using Plots, StatsPlots
using Genie
using PrecompileTools: @compile_workload

include("commons.jl")
# include("genomes.jl")
include("phenomes.jl")
# include("trials.jl")
# include("tebv.jl")

export PlotsGB, DistributionPlots, ViolinPlots
export checkdims, labeltofname, saveplots
export plotstatic, plotinteractive

@compile_workload begin
    true
end

end
