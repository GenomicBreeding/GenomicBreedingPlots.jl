"""
phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations .= "pop_1"; phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
x = plotstatic(DistributionPlots, phenomes)
fnames = saveplots(x)
rm.(fnames)
"""
function plotstatic(type::Type{T}, phenomes::Phenomes)::T where {T<:DistributionPlots}
    # type = DistributionPlots; phenomes = Phenomes(n=10, t=3); phenomes.entries = string.("entry_", 1:10); phenomes.populations .= "pop_1"; phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(10,3); phenomes.phenotypes[1, 1] = missing;
    df = GBCore.tabularise(phenomes)
    populations::Vector{String} = unique(df.populations)
    labels = Vector{String}(undef, length(phenomes.traits) + (length(phenomes.traits) * length(populations)))
    plots = Vector{Plots.Plot}(undef, length(labels))
    i::Int64 = 0
    # Across all populations
    for trait in Symbol.(phenomes.traits)
        # trait = Symbol.(phenomes.traits)[1]
        i += 1
        idx = findall(.!ismissing.(df[!, trait]))
        labels[i] = string("All populations\nTrait: ", trait, " (n=", length(idx), ")")
        plots[i] = StatsPlots.@df df[idx, :] StatsPlots.density(
            cols(trait),
            title = labels[i],
            xlabel = string(trait),
            ylabel = "Density",
            legend = false,
        )
    end
    # Per population
    for trait in Symbol.(phenomes.traits)
        for pop in populations
            # trait = Symbol.(phenomes.traits)[1]; pop = populations[1];
            i += 1
            idx = findall((.!ismissing.(df[!, trait])) .&& (phenomes.populations .== pop))
            labels[i] = string("Population: ", pop, "\nTrait: ", trait, " (n=", length(idx), ")")
            plots[i] = StatsPlots.@df df[idx, :] StatsPlots.density(
                cols(trait),
                title = labels[i],
                xlabel = string(trait),
                ylabel = "Density",
                legend = false,
            )
        end
    end
    out::type = DistributionPlots(labels, plots)
    out
end


"""
phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
x = plotstatic(ViolinPlots, phenomes)
fnames = saveplots(x)
rm.(fnames)
"""
function plotstatic(type::Type{T}, phenomes::Phenomes)::T where {T<:ViolinPlots}
    # type = DistributionPlots; phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
    if !GBCore.checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    df = GBCore.tabularise(phenomes)
    labels = Vector{String}(undef, 1 + length(phenomes.traits))
    plots = Vector{Plots.Plot}(undef, length(labels))
    i::Int64 = 0
    # Across all populations
    i += 1
    df_reshaped = DataFrames.stack(df, propertynames(df)[(end-(length(phenomes.traits)-1)):end])
    idx = findall(.!ismissing.(df_reshaped[!, :value]))
    lengths = combine(groupby(df_reshaped[idx, :], :variable), :value => length)
    for x in eachrow(lengths)
        # x = eachrow(lengths)[1]
        df_reshaped.variable[df_reshaped.variable.==x[1]] .= string(x[1], "\n(n=", x[2], ")")
    end
    labels[i] = "Trait distributions"
    p = StatsPlots.@df df_reshaped[idx, :] StatsPlots.violin(
        :variable,
        :value,
        title = labels[i],
        xlabel = "Traits",
        ylabel = "Phenotype values",
        legend = false,
    )
    StatsPlots.@df df_reshaped[idx, :] StatsPlots.boxplot!(:variable, :value, fill = (0.5))
    # StatsPlots.@df df_reshaped[idx, :] StatsPlots.dotplot!(:variable, :value, fill=(0.5)); gui(p)
    plots[i] = p
    # Per trait
    for trait in Symbol.(phenomes.traits)
        # trait = Symbol.(phenomes.traits)[1]
        i += 1
        df.pop = deepcopy(df.populations)
        idx = findall(.!ismissing.(df[!, trait]))
        lengths = combine(groupby(df[idx, :], :populations), trait => length)
        for x in eachrow(lengths)
            # x = eachrow(lengths)[1]
            df.pop[df.populations.==x[1]] .= string(x[1], "\n(n=", x[2], ")")
        end
        labels[i] = string("Trait: ", trait)
        p = StatsPlots.@df df[idx, :] StatsPlots.violin(
            :pop,
            cols(trait),
            title = labels[i],
            xlabel = "Populations",
            ylabel = string(trait),
            legend = false,
        )
        StatsPlots.@df df[idx, :] StatsPlots.boxplot!(:pop, cols(trait), fill = (0.5))
        # StatsPlots.@df df[idx, :] StatsPlots.dotplot!(:pop, cols(trait), fill=(0.5)); gui(p)
        plots[i] = p
    end
    out::type = ViolinPlots(labels, plots)
    out
end


function plotinteractive(phenomes::Phenomes)
    # Genie Hello World!
    # As simple as Hello
    route("/hello") do
        "Hello Genie!"
    end

    up(8888)

    down()

end
