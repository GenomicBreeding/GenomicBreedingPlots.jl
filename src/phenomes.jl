"""
    plot(
        type::Type{T},
        phenomes::Phenomes;
        plot_size::Tuple{Int64,Int64} = (600, 450),
    )::T where {T<:DistributionPlots}

Plot distributions of each trait across populations

# Examples
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;

julia> dplots = GBPlots.plot(DistributionPlots, phenomes);

julia> fnames = saveplots(dplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    phenomes::Phenomes;
    plot_size::Tuple{Int64,Int64} = (600, 450),
)::T where {T<:DistributionPlots}
    # type = DistributionPlots
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    traits::Vector{String} = sort(unique(phenomes.traits))
    labels = Vector{String}(undef, length(traits))
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    df = tabularise(phenomes)
    # Per trait across populations
    for (i, trait) in enumerate(traits)
        # i = 1; trait = traits[i]
        y = df[:, trait]
        idx = findall(.!ismissing.(y) .&& .!isnan.(y) .&& .!isinf.(y))
        fig = if length(idx) > 3
            y = y[idx]
            labels[i] = string("Trait: ", trait, " (n=", length(idx), ")")
            fig = CairoMakie.Figure(size = plot_size)
            axs = CairoMakie.Axis(fig[1, 1], title = labels[i])
            CairoMakie.density!(axs, y)
            fig
        end
        plots[i] = fig
    end
    # Output
    out::type = DistributionPlots(labels, plots)
    out
end


"""
    plot(type::Type{T}, phenomes::Phenomes; plot_size::Tuple{Int64, Int64} = (600, 450))::T where {T<:ViolinPlots}

Violin plots per trait per population

# Examples
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;

julia> vplots = GBPlots.plot(ViolinPlots, phenomes);

julia> fnames = saveplots(vplots)

julia> rm.(fnames);
```
"""
function plot(type::Type{T}, phenomes::Phenomes; plot_size::Tuple{Int64,Int64} = (600, 450))::T where {T<:ViolinPlots}
    # type = ViolinPlots
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    traits::Vector{String} = sort(unique(phenomes.traits))
    populations::Vector{String} = sort(unique(phenomes.populations))
    labels = Vector{String}(undef, length(traits))
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    df = tabularise(phenomes)
    # Reshape so that we have the following fields: entries, populations, traits, y, reversed population ID, and colours per population (needing because all x and y values need to be numeric) for violin plotting
    colours = repeat(colorschemes[:tol_light][1:end], Int(ceil(length(populations) / length(colorschemes[:tol_light]))))
    df_reshaped = begin
        df_reshaped = DataFrames.stack(df, propertynames(df)[(end-(length(phenomes.traits)-1)):end])
        rename!(df_reshaped, :variable => :traits)
        rename!(df_reshaped, :value => :y)
        df_reshaped.popid_reversed = [findall(reverse(populations) .== pop)[1] for pop in df_reshaped.populations]
        df_reshaped.popid_colors =
            [colours[findall(reverse(populations) .== pop)[1]] for pop in df_reshaped.populations]
        idx = findall(.!ismissing.(df_reshaped.y) .&& .!isnan.(df_reshaped.y) .&& .!isinf.(df_reshaped.y))
        df_reshaped[idx, :]
    end
    for (i, trait) in enumerate(traits)
        # i = 1; trait = traits[i]
        idx = findall(df_reshaped.traits .== trait)
        fig = if length(idx) > 3
            counts::Vector{Int64} = [sum(df_reshaped.populations .== population) for population in populations]
            labels[i] = string("Trait: ", trait, " (per population)")
            fig = CairoMakie.Figure(size = plot_size)
            axs = CairoMakie.Axis(
                fig[1, 1],
                title = labels[i],
                yticks = (1:length(populations), reverse(string.(populations, " (n=", counts, ")"))),
            )
            CairoMakie.violin!(
                axs,
                df_reshaped.popid_reversed[idx],
                df_reshaped.y[idx],
                color = df_reshaped.popid_colors[idx],
                orientation = :horizontal,
                transparency = 0.75,
                scale = :width,
            )
            CairoMakie.boxplot!(
                axs,
                df_reshaped.popid_reversed[idx],
                df_reshaped.y[idx],
                orientation = :horizontal,
                strokewidth = 1,
                outliercolor = :black,
                width = 1 / length(populations),
            )
            fig
        end
        plots[i] = fig
    end
    # Output
    out::type = ViolinPlots(labels, plots)
    out
end

"""
    plot(
        type::Type{T},
        phenomes::Phenomes;
        plot_size::Tuple{Int64, Int64} = (600, 450),
        color_scheme::Symbol = :viridis,
        rev_label_colors::Bool = false,
    )::T where {T<:CorHeatPlots}

Correlation heatmaps:
- between traits across populations
- between entries across populations
- between traits per populations
- between entries per populations

# Examples
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;

julia> hplots = GBPlots.plot(CorHeatPlots, phenomes);

julia> fnames = saveplots(hplots, format="png")

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    phenomes::Phenomes;
    plot_size::Tuple{Int64,Int64} = (600, 450),
    color_scheme::Symbol = :viridis,
    rev_label_colors::Bool = false,
    n_threshold_to_show_text::Int64 = 1_000,
)::T where {T<:CorHeatPlots}
    # type = CorHeatPlots
    # genomes = GBCore.simulategenomes(n=100, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing
    # plot_size = (700, 500); color_scheme = :viridis; rev_label_colors = false; n_threshold_to_show_text = 1_000
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    traits::Vector{String} = sort(unique(phenomes.traits))
    populations::Vector{String} = sort(unique(phenomes.populations))
    labels = Vector{String}(undef, 2 * (1 + length(populations)))
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    # Instantiate the vectors of correlation matrices, counts per pairwise correlation, labels and grouping names
    correlations = Vector{Matrix{Float64}}(undef, length(labels))
    counts = Vector{Matrix{Int64}}(undef, length(labels))
    labellings = Vector{Vector{String}}(undef, length(labels))
    groupings = Vector{String}(undef, length(labels))
    # Population the vectors above
    traits, entries, dist = distances(phenomes, distance_metrics = ["correlation"], standardise_traits = true)
    correlations[1:2] = [dist["traits|correlation"], dist["entries|correlation"]]
    counts[1:2] = [Int.(dist["traits|counts"]), Int.(dist["entries|counts"])]
    labellings[1:2] = [traits, entries]
    groupings[1:2] = ["between traits across all populations", "between entries across all populations"]
    i = 3
    for population in populations
        idx_entries = findall(phenomes.populations .== population)
        traits_per_pop, entries_per_pop, dist_per_pop =
            distances(slice(phenomes, idx_entries = idx_entries), distance_metrics = ["correlation"])
        correlations[i:(i+1)] = [dist_per_pop["traits|correlation"], dist_per_pop["entries|correlation"]]
        counts[i:(i+1)] = [Int.(dist_per_pop["traits|counts"]), Int.(dist_per_pop["entries|counts"])]
        labellings[i:(i+1)] = [traits_per_pop, entries_per_pop]
        groupings[i:(i+1)] = [string("between traits (", population, ")"), string("between entries (", population, ")")]
        i += 2
    end
    # Plot heatmaps
    for (i, (C, N, lab, grp)) in enumerate(zip(correlations, counts, labellings, groupings))
        # i = 1; C=correlations[i]; N=counts[i]; lab=labellings[i]; grp=groupings[i];
        n = size(C, 1)
        font_size_ticks = minimum([20, 15 / (0.1 * n)])
        font_size_labels = minimum([20, 9 / (0.1 * n)])
        x = repeat(1:n, inner = n)
        y = repeat(1:n, outer = n)
        c = string.(round.(reshape(C, n * n, 1)[:, 1], digits = 2), "\n(n=", reshape(N, n * n, 1)[:, 1], ")")
        col = if !rev_label_colors
            [x > 0.5 ? :black : :white for x in reshape(C, n * n, 1)[:, 1]]
        else
            [x < 0.5 ? :black : :white for x in reshape(C, n * n, 1)[:, 1]]
        end
        labels[i] = string("Correlations ", grp)
        fig = CairoMakie.Figure(size = plot_size)
        axs = CairoMakie.Axis(
            fig[1, 1],
            xticks = (1:n, lab),
            yticks = (1:n, lab),
            xticklabelrotation = deg2rad(90),
            yreversed = true,
            xaxisposition = :top,
            yaxisposition = :right,
            xticklabelsize = font_size_ticks,
            yticklabelsize = font_size_ticks,
        )
        CairoMakie.Label(fig[2, 1], labels[i])
        CairoMakie.heatmap!(axs, 1:n, 1:n, C, colorrange = (-1.0, 1.0), colormap = color_scheme)
        if n^2 < n_threshold_to_show_text
            CairoMakie.text!(axs, x, y, text = c, align = (:center, :center), color = col, fontsize = font_size_labels)
        end
        plots[i] = fig
    end
    # Output
    out::type = CorHeatPlots(labels, plots)
    out
end


"""
    plot(
        type::Type{T},
        phenomes::Phenomes;
        plot_size::Tuple{Int64,Int64} = (600, 450),
        horizontal::Bool = true,
        standardise_traits::Bool = true,
    )::T where {T<:TreePlots}

Plot tree diagrams showing the relationships of each entry using trait information.

- Distance metric: Euclidean
- Grouping/linkage: Ward's distance
- Branch order: Optimal (Bar-Joseph et al, 2001. Bionformatics.)

# Examples
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;

julia> tplots = GBPlots.plot(TreePlots, phenomes);

julia> fnames = saveplots(tplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    phenomes::Phenomes;
    plot_size::Tuple{Int64,Int64} = (600, 450),
    horizontal::Bool = true,
    standardise_traits::Bool = true,
)::T where {T<:TreePlots}
    # type = TreePlots
    # genomes = GBCore.simulategenomes(n=100, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;
    # plot_size = (700, 500); horizontal = true; standardise_traits = true;
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    traits::Vector{String} = sort(unique(phenomes.traits))
    populations::Vector{String} = sort(unique(phenomes.populations))
    labels = Vector{String}(undef, 2)
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    # Estimate pairwise Euclidean distances
    traits, entries, dist = distances(phenomes, distance_metrics = ["euclidean"])
    clusters = [
        Clustering.hclust(dist["traits|euclidean"], linkage = :ward, branchorder = :optimal),
        Clustering.hclust(dist["entries|euclidean"], linkage = :ward, branchorder = :optimal),
    ]
    leaves = [phenomes.traits, phenomes.entries]
    titles = ["Traits dendrogram", "Entries dendrogram"]
    for (i, (clust, tips, title)) in enumerate(zip(clusters, leaves, titles))
        # i = 2; clust = clusters[i]; tips = leaves[i]; title = titles[i];
        n = length(tips)
        m = length(clust.heights)
        font_size_ticks = minimum([20, 50 / (0.1 * n)])
        fig = CairoMakie.Figure(size = plot_size)
        axs = if !horizontal
            CairoMakie.Axis(
                fig[1, 1],
                title = title,
                xgridvisible = false,
                ylabel = "Euclidean Distance",
                xticklabelsize = font_size_ticks,
                xticks = (1:n, tips[clust.order]),
                xticklabelrotation = deg2rad(90),
            )
        else
            CairoMakie.Axis(
                fig[1, 1],
                title = title,
                ygridvisible = false,
                xlabel = "Euclidean Distance",
                yticklabelsize = font_size_ticks,
                yticks = (1:n, tips[clust.order]),
            )
        end
        colours = repeat(colorschemes[:tol_light][1:end], Int(ceil(m / length(colorschemes[:tol_light]))))
        midpoints = fill(0.0, m)
        for i = 1:m
            # i = 1
            # println(i)
            x1, x2 = (clust.merges[i, 1], clust.merges[i, 2])
            d = clust.heights[i]
            dx1 = x1 < 0 ? 0.0 : clust.heights[x1]
            dx2 = x2 < 0 ? 0.0 : clust.heights[x2]
            x1 = x1 < 0 ? findall(clust.order .== -x1)[1] : midpoints[x1]
            x2 = x2 < 0 ? findall(clust.order .== -x2)[1] : midpoints[x2]
            midpoints[i] = mean(abs.([x1, x2]))
            if !horizontal
                CairoMakie.lines!(axs, [x1, x1], [dx1, d], color = colours[i])
                CairoMakie.lines!(axs, [x2, x2], [dx2, d], color = colours[i])
                CairoMakie.lines!(axs, [x1, x2], [d, d], color = colours[i])
            else
                CairoMakie.lines!(axs, [dx1, d], [x1, x1], color = colours[i])
                CairoMakie.lines!(axs, [dx2, d], [x2, x2], color = colours[i])
                CairoMakie.lines!(axs, [d, d], [x1, x2], color = colours[i])
            end
        end
        labels[i] = title
        plots[i] = fig
    end
    # Output
    out::type = TreePlots(labels, plots)
    out
end
