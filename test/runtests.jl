using GBPlots
using GBCore
using StatsBase
using Test
using Documenter

Documenter.doctest(GBPlots)

@testset "GBPlots.jl" begin
    # Phenomes
    genomes = GBCore.simulategenomes(n = 300, verbose = false)
    genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace = true)
    trials, _ = GBCore.simulatetrials(
        genomes = genomes,
        n_years = 1,
        n_seasons = 1,
        n_harvests = 1,
        n_sites = 1,
        n_replications = 1,
        verbose = false,
    )
    phenomes = extractphenomes(trials)
    phenomes.phenotypes[1, 1] = missing
    plots = GBPlots.plot(DistributionPlots, phenomes)
    fnames = saveplots(plots)
    @test length(fnames) == length(plots.plots)
    rm.(fnames)
    plots = GBPlots.plot(ViolinPlots, phenomes)
    fnames = saveplots(plots)
    @test length(fnames) == length(plots.plots)
    rm.(fnames)
    plots = GBPlots.plot(CorHeatPlots, phenomes)
    fnames = saveplots(plots)
    @test length(fnames) == length(plots.plots)
    rm.(fnames)
    plots = GBPlots.plot(TreePlots, phenomes)
    fnames = saveplots(plots)
    @test length(fnames) == length(plots.plots)
    rm.(fnames)
    # CV
    cvs::Vector{CV} = []
    for m = 1:3
        for t = 1:5
            for p = 1:4
                n = 100
                populations = if rand() < 0.5
                    string.("population_", sample(1:4, n, replace = true))
                else
                    string.("population_", repeat([p], n))
                end
                entries = string.("entry_", sample(1:1_000, n, replace = true))
                for r = 1:5
                    for f = 1:5
                        fit = Fit(n = 10, l = 1_000)
                        fit.model = string("model_", m)
                        fit.populations .= string("population_", p)
                        fit.trait = string("trait_", t)
                        fit.metrics = Dict("cor" => rand() / maximum([1, 5 * rand()]), "rmse" => rand())
                        cv = CV(
                            string("replication_", r),
                            string("fold_", f),
                            fit,
                            populations,
                            entries,
                            rand(n),
                            rand(n),
                            fit.metrics,
                        )
                        push!(cvs, cv)
                    end
                end
            end
        end
    end
    plots = GBPlots.plot(BarPlots, cvs)
    fnames = saveplots(plots)
    @test length(fnames) == length(plots.plots)
    rm.(fnames)
    plots = GBPlots.plot(BoxPlots, cvs)
    fnames = saveplots(plots)
    @test length(fnames) == length(plots.plots)
    rm.(fnames)

end
