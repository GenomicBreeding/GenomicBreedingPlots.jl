"""
    plotstatic(type::Type{T}, phenomes::Phenomes)::T where {T<:DistributionPlots}

Plot distributions of each trait across populations and per population

# Examples
```julia
phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations .= "pop_1"; phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
x = plotstatic(DistributionPlots, phenomes)
fnames = saveplots(x)
rm.(fnames)
```
"""
function plotstatic(type::Type{T}, phenomes::Phenomes)::T where {T<:DistributionPlots}
    # type = DistributionPlots; phenomes = Phenomes(n=10, t=3); phenomes.entries = string.("entry_", 1:10); phenomes.populations .= "pop_1"; phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(10,3); phenomes.phenotypes[1, 1] = missing;
    df = tabularise(phenomes)
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
    plotstatic(type::Type{T}, phenomes::Phenomes)::T where {T<:ViolinPlots}

Violin plots across populations and per population across traits

# Examples
```julia
phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
x = plotstatic(ViolinPlots, phenomes)
fnames = saveplots(x)
rm.(fnames)
```
"""
function plotstatic(type::Type{T}, phenomes::Phenomes)::T where {T<:ViolinPlots}
    # type = DistributionPlots; phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    df = tabularise(phenomes)
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

"""
    plotstatic(type::Type{T}, phenomes::Phenomes; color_scheme::Symbol=:YlGn)::T where {T<:CorHeatPlots}

Trait correlation heatmaps across populations and per population

# Examples
```julia
phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
x = plotstatic(CorHeatPlots, phenomes, color_scheme=:viridis)
fnames = saveplots(x)
rm.(fnames)
```
"""
function plotstatic(type::Type{T}, phenomes::Phenomes; color_scheme::Symbol = :YlGn)::T where {T<:CorHeatPlots}
    # type = CorHeatPlots; phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["trait_1", "trait_2", "long_trait_name number 3"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing; color_scheme::Symbol=:YlGn
    # type = CorHeatPlots; phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["trait_1", "trait_2", "trait_3"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing; color_scheme::Symbol=:YlGn
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    df = tabularise(phenomes)
    populations = unique(phenomes.populations)
    labels = Vector{String}(undef, 1 + length(populations))
    plots = Vector{Plots.Plot}(undef, length(labels))
    # Plot across populations and per population
    i::Int64 = 0
    for pop in vcat("All populations", populations)
        # pop = populations[1]
        i += 1
        if pop == "All populations"
            idx = collect(1:size(df, 1))
        else
            idx = findall(df.populations .== pop)
        end
        X = Matrix(df[idx, Symbol.(phenomes.traits)])
        idx = findall(.!ismissing.(sum(X, dims = 2)[:, 1]))
        labels[i] = string("Population: ", pop, "\n(n=", length(idx), ")")
        C = cor(X[idx, :])
        C = C[reverse(collect(1:end)), :]
        p = Plots.heatmap(
            C,
            color = color_scheme,
            title = labels[i],
            titlefont = font(10),
            clim = (-1, 1),
            xmirror = true,
            xrotation = -45,
            yrotation = 0,
            xticks = (collect(1:length(phenomes.traits)), phenomes.traits),
            yticks = (reverse(collect(1:length(phenomes.traits))), phenomes.traits),
            top_margin = (maximum(length.(phenomes.traits)))mm,
        )
        for i in axes(C, 1)
            for j in axes(C, 2)
                if C[i, j] <= -0.5
                    annotate!(p, (j, i, text(round(C[i, j], digits = 3), 8, :red)))
                elseif C[i, j] < 0.0
                    annotate!(p, (j, i, text(round(C[i, j], digits = 3), 8, :purple)))
                elseif C[i, j] < 0.5
                    annotate!(p, (j, i, text(round(C[i, j], digits = 3), 8, :black)))
                else
                    annotate!(p, (j, i, text(round(C[i, j], digits = 3), 8, :white)))
                end
            end
        end
        # gui(p)
        plots[i] = p
    end
    out::type = CorHeatPlots(labels, plots)
    out
end


"""
    plotstatic(type::Type{T}, phenomes::Phenomes; color_scheme::Symbol=:YlGn)::T where {T<:TreePlots}

Plot tree diagrams showing the relationships of each entry using trait information.

- Distance metric: Euclidean
- Grouping/linkage: Ward's distance
- Branch order: Optimal (Bar-Joseph et al, 2001. Bionformatics.)

# Examples
```julia
phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; 
Ψ::Matrix{Float64} = phenomes.phenotypes;
groups = Clustering.kmeans(Ψ', 5);
phenomes.populations = string.("pop_", groups.assignments); 
phenomes.phenotypes[1, 1] = missing;
x = plotstatic(TreePlots, phenomes)
fnames = saveplots(x)
rm.(fnames)
```
"""
function plotstatic(type::Type{T}, phenomes::Phenomes; color_scheme::Symbol = :default)::T where {T<:TreePlots}
    # type = TreePlots; phenomes = Phenomes(n=10, t=3); phenomes.entries = string.("entry_", 1:10); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 10, replace=true); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 10)'; phenomes.phenotypes[1, 1] = missing; color_scheme::Symbol=:default;
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    df = tabularise(phenomes)
    populations = unique(phenomes.populations)
    labels::Vector{String} = fill("", 1 + length(populations))
    plots = Vector{Plots.Plot}(undef, length(labels))
    # Plot across populations and per population
    colors_per_population = palette(color_scheme)[1:length(populations)]
    colors = Vector{RGB{Float64}}(undef, length(phenomes.populations))
    for i in eachindex(phenomes.populations)
        # i, p = 1, phenomes.populations[1]
        colors[i] = colors_per_population[findall(populations .== phenomes.populations[i])[1]]
    end
    i::Int64 = 0
    for pop in vcat("All populations", populations)
        # pop = "All populations"
        # pop = populations[1]
        # println(pop)
        i += 1
        if pop == "All populations"
            idx_entries_1 = collect(1:size(df, 1))
        else
            idx_entries_1 = findall(df.populations .== pop)
        end
        X = Matrix(df[idx_entries_1, Symbol.(phenomes.traits)])
        # Omit entries with at least one missing phenotype
        idx_entries_2 = findall(.!ismissing.(sum(X, dims = 2)[:, 1]))
        if length(idx_entries_2) < 2
            continue
        end
        X = X[idx_entries_2, :]
        # Omit traits with zero variance
        idx_traits = findall(std(X, dims = 1)[1, :] .> 1e-20)
        if length(idx_traits) < 2
            continue
        end
        X = X[:, idx_traits]
        # Normalise
        X = (X .- mean(X, dims = 1)) ./ std(X, dims = 1)
        # Plot
        labels[i] = string("Population: ", pop, "\n(n=", size(X, 1), "; t=", size(X, 2), ")")
        distances = Distances.pairwise(Distances.Euclidean(), X, dims = 1)
        clusters = Clustering.hclust(distances, linkage = :ward, branchorder = :optimal)
        p = Plots.plot(
            clusters,
            xaxis = false,
            xticks = (1:length(idx_entries_2), phenomes.entries[idx_entries_1][idx_entries_2][clusters.order]),
            yrotation = 0,
            yflip = true,
            title = labels[i],
            ylabel = "Euclidean distance",
            permute = (:x, :y),
            size = (400 + 10 * maximum(length.(phenomes.entries)), 600 + 10 * length(phenomes.entries)),
            top_margin = 10mm,
            right_margin = 1.75 * maximum(length.(phenomes.entries))mm,
        )
        cols = colors[idx_entries_1][idx_entries_2][clusters.order]
        yt, yl = Plots.yticks(p)[1]
        y0 = zeros(length(yl)) .- 0.05
        yticks!(yt, fill(" ", size(phenomes.entries)))
        for (yi, xi, li, ci) in zip(yt, y0, yl, cols)
            annotate!(xi, yi, text(li, 8, ci, :left, rotation = 0))
        end
        # Plots.gui(p)
        plots[i] = p
    end
    # Remove skipped populations
    idx = findall(labels .!= "")
    # Output
    out::type = TreePlots(labels[idx], plots[idx])
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
