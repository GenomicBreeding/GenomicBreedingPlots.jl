using Pkg
Pkg.activate(".")
using GenomicBreedingPlots
using GenomicBreedingCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra
using DataFrames, Random
using Distances, Clustering, Measures
using CairoMakie, ColorSchemes
