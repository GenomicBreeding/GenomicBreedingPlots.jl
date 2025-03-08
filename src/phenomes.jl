"""
    plot(type::Type{T}, phenomes::Phenomes; plot_size::Tuple{Int64,Int64}=(600, 450))::T where {T<:DistributionPlots}

Create density plots for each trait in the provided phenomes data.

# Arguments
- `type::Type{T}`: The type of distribution plot to create (must be a subtype of `DistributionPlots`).
- `phenomes::Phenomes`: A Phenomes struct containing the phenotypic data to plot.
- `plot_size::Tuple{Int64,Int64}`: Optional tuple specifying the width and height of each plot in pixels. 
  Defaults to (600, 450).

# Returns
- `T<:DistributionPlots`: A DistributionPlots object containing the generated density plots and their labels.

# Details
The function creates density plots for each trait in the phenomes data, excluding any missing, NaN, or infinite 
values. Only traits with more than 3 valid observations will be plotted. Each plot includes the trait name and 
the number of observations in its title.

# Throws
- `ArgumentError`: If the dimensions of the Phenomes struct are invalid or corrupted.

# Example
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
        if length(idx) > 3
            y = y[idx]
            labels[i] = string("Trait: ", trait, " (n=", length(idx), ")")
            fig = CairoMakie.Figure(size = plot_size)
            axs = CairoMakie.Axis(fig[1, 1], title = labels[i])
            CairoMakie.density!(axs, y)
            plots[i] = fig
        end
    end
    # Output
    out::type = DistributionPlots(labels, plots)
    out
end



"""
    plot(type::Type{T}, phenomes::Phenomes; plot_size::Tuple{Int64,Int64}=(600, 450), colour_scheme::Symbol=:viridis)::T where {T<:ViolinPlots}

Create violin plots for phenotypic data across different populations.

# Arguments
- `type::Type{T}`: The type of plot to create (must be a subtype of ViolinPlots)
- `phenomes::Phenomes`: A Phenomes struct containing the phenotypic data to visualize
- `plot_size::Tuple{Int64,Int64}`: Size of the plot in pixels (width, height)
- `colour_scheme::Symbol`: Color scheme for the violin plots (default: :viridis)

# Returns
- `T<:ViolinPlots`: A ViolinPlots object containing the generated plots and their labels

# Details
The function creates violin plots for each trait in the phenomes data, showing the distribution
of phenotypic values across different populations. Each plot includes:
- Violin plots showing the distribution density
- Box plots overlaid on the violin plots
- Population counts in the y-axis labels
- Horizontal orientation for better readability

Missing, NaN, and infinite values are automatically filtered out.

# Throws
- `ArgumentError`: If the Phenomes struct dimensions are corrupted

# Example
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;

julia> vplots = GBPlots.plot(ViolinPlots, phenomes);

julia> fnames = saveplots(vplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    phenomes::Phenomes;
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :viridis,
)::T where {T<:ViolinPlots}
    # type = ViolinPlots
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing
    # plot_size = (600, 450); colour_scheme = :viridis;
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    traits::Vector{String} = sort(unique(phenomes.traits))
    populations::Vector{String} = sort(unique(phenomes.populations))
    labels = Vector{String}(undef, length(traits))
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    df = tabularise(phenomes)
    # Reshape so that we have the following fields: entries, populations, traits, y, reversed population ID, and colours per population (needing because all x and y values need to be numeric) for violin plotting
    colours =
        repeat(colorschemes[colour_scheme][1:end], Int(ceil(length(populations) / length(colorschemes[colour_scheme]))))
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
        if length(idx) > 3
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
            plots[i] = fig
        end
    end
    # Output
    out::type = ViolinPlots(labels, plots)
    out
end


"""
    plot(type::Type{T}, phenomes::Phenomes; 
         plot_size::Tuple{Int64,Int64}=(600, 450),
         colour_scheme::Symbol=:viridis,
         rev_label_colors::Bool=false,
         n_threshold_to_show_text::Int64=1_000)::T where {T<:CorHeatPlots}

Generate correlation heatmap plots for phenotypic data.

Creates multiple correlation heatmaps showing relationships between traits and entries
across all populations and within individual populations.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be a subtype of CorHeatPlots)
- `phenomes::Phenomes`: A Phenomes struct containing the phenotypic data to plot

# Keywords
- `plot_size::Tuple{Int64,Int64}`: Size of the plot in pixels (width, height)
- `colour_scheme::Symbol`: Color scheme for the heatmap (e.g., :viridis)
- `rev_label_colors::Bool`: If true, reverses the text color threshold for labels
- `n_threshold_to_show_text::Int64`: Maximum number of cells before text labels are hidden

# Returns
- `CorHeatPlots`: A CorHeatPlots object containing correlation heatmaps for:
  * Traits across all populations
  * Entries across all populations
  * Traits within each population
  * Entries within each population

# Notes
- Correlation values range from -1.0 to 1.0
- Text labels show correlation values and sample sizes (when sample sizes vary)
- Text color automatically adjusts for readability against the background
- Labels are hidden when the number of cells exceeds `n_threshold_to_show_text`

# Throws
- `ArgumentError`: If the Phenomes struct dimensions are corrupted

# Example
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
    colour_scheme::Symbol = :viridis,
    rev_label_colors::Bool = false,
    n_threshold_to_show_text::Int64 = 1_000,
)::T where {T<:CorHeatPlots}
    # type = CorHeatPlots
    # genomes = GBCore.simulategenomes(n=100, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing
    # plot_size = (700, 500); colour_scheme = :viridis; rev_label_colors = false; n_threshold_to_show_text = 1_000
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
    # Population vectors above
    traits, entries, dist = distances(phenomes, distance_metrics = ["correlation"], standardise_traits = true)
    correlations[1:2] = [dist["traits|correlation"], dist["entries|correlation"]]
    counts[1:2] = [Int.(dist["traits|counts"]), Int.(dist["entries|counts"])]
    labellings[1:2] = [traits, entries]
    groupings[1:2] = ["between traits across all populations", "between entries across all populations"]
    i = 3
    for population in populations
        idx_entries = findall(phenomes.populations .== population)
        traits_per_pop, entries_per_pop, dist_per_pop = try
            distances(slice(phenomes, idx_entries = idx_entries), distance_metrics = ["correlation"])
        catch
            continue
        end
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
        c = if length(unique(N)) > 1
            string.(round.(reshape(C, n * n, 1)[:, 1], digits = 2), "\n(n=", reshape(N, n * n, 1)[:, 1], ")")
        else
            string.(round.(reshape(C, n * n, 1)[:, 1], digits = 2))
        end
        col = if !rev_label_colors
            [x > 0.5 ? :black : :white for x in reshape(C, n * n, 1)[:, 1]]
        else
            [x < 0.5 ? :black : :white for x in reshape(C, n * n, 1)[:, 1]]
        end
        labels[i] = string("Correlations (phenotypes) ", grp)
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
        CairoMakie.heatmap!(axs, 1:n, 1:n, C, colorrange = (-1.0, 1.0), colormap = colour_scheme)
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
    plot(type::Type{T}, phenomes::Phenomes; plot_size::Tuple{Int64,Int64}=(600, 450), 
         colour_scheme::Symbol=:tol_light, horizontal::Bool=true, 
         standardise_traits::Bool=true)::T where {T<:TreePlots}

Generate hierarchical clustering dendrograms for traits and entries based on phenotypic data.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be a subtype of TreePlots)
- `phenomes::Phenomes`: A Phenomes struct containing phenotypic data
- `plot_size::Tuple{Int64,Int64}`: Size of the output plots in pixels (width, height)
- `colour_scheme::Symbol`: Color scheme to use for the dendrogram branches
- `horizontal::Bool`: If true, plots dendrograms horizontally; if false, vertically
- `standardise_traits::Bool`: If true, standardizes trait values before distance calculation

# Returns
- A `TreePlots` object containing two dendrograms:
  1. Traits dendrogram showing relationships between phenotypic traits
  2. Entries dendrogram showing relationships between genetic entries

# Details
The function performs the following operations:
- Calculates pairwise Euclidean distances between traits and entries
- Performs hierarchical clustering using Ward's method
- Creates dendrograms with optimal branch ordering
- Adjusts font sizes automatically based on the number of labels
- Uses color schemes from ColorSchemes.jl

# Throws
- `ArgumentError`: If the dimensions in the Phenomes struct are corrupted

# Example
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
    colour_scheme::Symbol = :tol_light,
    horizontal::Bool = true,
    standardise_traits::Bool = true,
)::T where {T<:TreePlots}
    # type = TreePlots
    # genomes = GBCore.simulategenomes(n=100, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;
    # plot_size = (700, 500); colour_scheme = :tol_light; horizontal = true; standardise_traits = true;
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    traits::Vector{String} = sort(unique(phenomes.traits))
    populations::Vector{String} = sort(unique(phenomes.populations))
    labels = Vector{String}(undef, 2)
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    # Estimate pairwise Euclidean distances
    traits, entries, dist =
        distances(phenomes, distance_metrics = ["euclidean"], standardise_traits = standardise_traits)
    clusters = [
        Clustering.hclust(dist["traits|euclidean"], linkage = :ward, branchorder = :optimal),
        Clustering.hclust(dist["entries|euclidean"], linkage = :ward, branchorder = :optimal),
    ]
    leaves = [phenomes.traits, phenomes.entries]
    titles = ["Traits dendrogram (phenotypes)", "Entries dendrogram (phenotypes)"]
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
        colours = repeat(colorschemes[colour_scheme][1:end], Int(ceil(m / length(colorschemes[colour_scheme]))))
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

"""
    plot(type::Type{T}, phenomes::Phenomes; plot_size::Tuple{Int64,Int64}=(600, 450), colour_scheme::Symbol=:tol_muted)::T where {T<:PCBiPlots}

Generate Principal Component Analysis (PCA) biplots for phenotypic data.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be a subtype of PCBiPlots)
- `phenomes::Phenomes`: A Phenomes struct containing phenotypic data, populations, and traits
- `plot_size::Tuple{Int64,Int64}`: Size of the output plots in pixels (width, height)
- `colour_scheme::Symbol`: Color scheme to use for the plots (default: :tol_muted)

# Returns
- `PCBiPlots`: A PCBiPlots object containing two figures:
    1. PCA biplot of entries colored by population
    2. PCA biplot of traits

# Details
The function performs the following operations:
- Handles missing, NaN, and Inf values in the phenotypic data
- Standardizes the data (centers and scales)
- Performs PCA when more than 2 dimensions are present
- Creates direct trait comparison plots when only 2 dimensions are present
- Adds legends when multiple populations or traits are present

# Throws
- `ArgumentError`: If the Phenomes struct dimensions are corrupted
- `ErrorException`: If the phenotypic data is too sparse

# Examples
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> pplots = GBPlots.plot(PCBiPlots, genomes);

julia> fnames = saveplots(pplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    phenomes::Phenomes;
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :tol_muted,
)::T where {T<:PCBiPlots}
    # type = PCBiPlots
    # genomes = GBCore.simulategenomes(n=90, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=2, n_seasons=2, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;
    # plot_size = (700, 500); colour_scheme::Symbol = :tol_muted;
    if !checkdims(phenomes)
        throw(ArgumentError("Phenomes struct is corrupted."))
    end
    populations::Vector{String} = sort(unique(phenomes.populations))
    traits::Vector{String} = sort(phenomes.traits)
    labels = Vector{String}(undef, 2)
    # Prepare the allele frequencie matrix    
    Y, traits = begin
        n, t = size(phenomes.phenotypes)
        y1 = mean(phenomes.phenotypes, dims = 2)[:, 1]
        y2 = mean(phenomes.phenotypes, dims = 1)[1, :]
        idx_rows = findall(.!ismissing.(y1) .&& .!isnan.(y1) .&& .!isinf.(y1))
        idx_cols = findall(.!ismissing.(y2) .&& .!isnan.(y2) .&& .!isinf.(y2))
        A, traits = if (length(idx_cols) == 0) || (length(idx_rows) / n >= length(idx_cols) / t)
            A = phenomes.phenotypes[idx_rows, :]
            a = mean(A, dims = 1)[1, :]
            idx_cols = findall(.!ismissing.(a) .&& .!isnan.(a) .&& .!isinf.(a))
            (A[:, idx_cols], phenomes.traits[idx_cols])
        elseif (length(idx_rows) == 0) || (length(idx_rows) / n <= length(idx_cols) / t)
            A = phenomes.phenotypes[:, idx_cols]
            a = mean(A, dims = 2)[:, 1]
            idx_rows = findall(.!ismissing.(a) .&& .!isnan.(a) .&& .!isinf.(a))
            (A[idx_rows, :], phenomes.traits[idx_cols])
        else
            throw(ErrorException("Phenomes struct too sparse"))
        end
        # Remove columns with zero variance
        μ = mean(A, dims = 1)
        σ = std(A, dims = 1)
        idx = findall(σ[1, :] .> 0.0)
        A = A[:, idx]
        μ = mean(A, dims = 1)
        σ = std(A, dims = 1)
        Y::Matrix{Float64} = (A .- μ) ./ σ
        (Y, traits)
    end
    labels = ["PCA (phenotypes) biplot of entries", "PCA (phenotypes) biplot of traits"]
    fig_entries = begin
        colours = [findall(populations .== pop)[1] for pop in phenomes.populations]
        fig = CairoMakie.Figure(size = plot_size)
        plt = if size(Y, 2) > 2
            M = MultivariateStats.fit(PCA, Y)
            pc1 = M.proj[:, 1]
            pc2 = M.proj[:, 2]
            variances_explained = M.prinvars ./ sum(M.prinvars)
            axs = CairoMakie.Axis(
                fig[1, 1],
                title = labels[1],
                xlabel = string("PC1 (", round(100 * variances_explained[1], digits = 2), "%)"),
                ylabel = string("PC2 (", round(100 * variances_explained[2], digits = 2), "%)"),
            )
            CairoMakie.scatter!(
                axs,
                pc1,
                pc2,
                color = colours,
                colormap = colour_scheme,
                colorrange = (1, maximum([2, length(populations)])),
            )
        else
            labels[1] = string(traits[1], " vs ", traits[2])
            x = Y[:, 1]
            y = Y[:, 2]
            axs = CairoMakie.Axis(
                fig[1, 1],
                title = labels[1],
                xlabel = phenomes.traits[idx_cols[1]],
                ylabel = phenomes.traits[idx_cols[2]],
            )
            CairoMakie.scatter!(
                axs,
                x,
                y,
                color = colours,
                colormap = colour_scheme,
                colorrange = (1, maximum([2, length(populations)])),
            )
        end
        if length(populations) > 1
            colourmap = getproperty(ColorSchemes, plt.colormap[])
            colours = colourmap[range(start = 0.0, stop = 1.0; length = length(populations))]
            elems = [
                [MarkerElement(color = col, marker = :circle, markersize = 15, strokecolor = :black)] for col in colours
            ]
            CairoMakie.Legend(fig[1, 2], elems, populations)
        end
        fig
    end
    fig_traits = begin
        colours = [findall(traits .== trait)[1] for trait in phenomes.traits]
        fig = CairoMakie.Figure(size = plot_size)
        plt = if size(Y, 1) > 2
            M = MultivariateStats.fit(PCA, Y')
            pc1 = M.proj[:, 1]
            pc2 = M.proj[:, 2]
            variances_explained = M.prinvars ./ sum(M.prinvars)
            axs = CairoMakie.Axis(
                fig[1, 1],
                title = labels[1],
                xlabel = string("PC1 (", round(100 * variances_explained[1], digits = 2), "%)"),
                ylabel = string("PC2 (", round(100 * variances_explained[2], digits = 2), "%)"),
            )
            CairoMakie.scatter!(
                axs,
                pc1,
                pc2,
                color = colours,
                colormap = colour_scheme,
                colorrange = (1, maximum([2, length(traits)])),
            )
        else
            labels[1] = string(phenomes.traits[idx_cols[1]], " vs ", phenomes.traits[idx_cols[2]])
            x = phenomes.phenotypes[:, idx_cols[1]]
            y = phenomes.phenotypes[:, idx_cols[2]]
            axs = CairoMakie.Axis(
                fig[1, 1],
                title = labels[1],
                xlabel = phenomes.traits[idx_cols[1]],
                ylabel = phenomes.traits[idx_cols[2]],
            )
            CairoMakie.scatter!(
                axs,
                x,
                y,
                color = colours,
                colormap = colour_scheme,
                colorrange = (1, maximum([2, length(traits)])),
            )
        end
        if length(traits) > 1
            colourmap = getproperty(ColorSchemes, plt.colormap[])
            colours = colourmap[range(start = 0.0, stop = 1.0; length = length(traits))]
            elems = [
                [MarkerElement(color = col, marker = :circle, markersize = 15, strokecolor = :black)] for col in colours
            ]
            CairoMakie.Legend(fig[1, 2], elems, traits)
        end
        fig
    end
    # Output
    out::type = PCBiPlots(labels, [fig_entries, fig_traits])
    out
end
