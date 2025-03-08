"""
    plot(type::Type{T}, genomes::Genomes; n_loci_alleles::Int64 = 1_000, seed::Int64 = 42, plot_size::Tuple{Int64,Int64} = (600, 450))::T where {T<:DistributionPlots}

Generate a density plot of allele frequencies distribution from genomic data.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be a subtype of DistributionPlots)
- `genomes::Genomes`: A Genomes struct containing the genomic data to plot
- `n_loci_alleles::Int64`: Number of loci to sample for plotting (default: 1000)
- `seed::Int64`: Random seed for reproducibility (default: 42)
- `plot_size::Tuple{Int64,Int64}`: Size of the output plot in pixels (width, height) (default: (600, 450))

# Returns
- `T<:DistributionPlots`: A DistributionPlots object containing the generated plot

# Description
This function creates a density plot showing the distribution of allele frequencies
across the sampled loci. It randomly samples a specified number of loci from the
genomic data and plots their frequency distribution using CairoMakie.

# Throws
- `ArgumentError`: If the Genomes struct dimensions are invalid or corrupted

# Example
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> dplots = GBPlots.plot(DistributionPlots, genomes);

julia> fnames = saveplots(dplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    genomes::Genomes;
    n_loci_alleles::Int64 = 1_000,
    seed::Int64 = 42,
    plot_size::Tuple{Int64,Int64} = (600, 450),
)::T where {T<:DistributionPlots}
    # type = DistributionPlots
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # n_loci_alleles = 1_000; seed = 42; plot_size = (600, 450);
    if !checkdims(genomes)
        throw(ArgumentError("Genomes struct is corrupted."))
    end
    rng::TaskLocalRNG = Random.seed!(seed)
    n, p = size(genomes.allele_frequencies)
    idx_loci_alleles::Vector{Int64} = sort(sample(rng, 1:p, minimum([n_loci_alleles, p]), replace = false))
    labels = ["Distribution of allele frequencies"]
    fig = CairoMakie.Figure(size = plot_size)
    axs = CairoMakie.Axis(fig[1, 1], title = labels[1])
    q::Vector{Float64} = []
    for i = 1:n
        for j in idx_loci_alleles
            qij = genomes.allele_frequencies[i, j]
            if !ismissing(qij) && !isnan(qij) && !isinf(qij)
                append!(q, qij)
            end
        end
    end
    CairoMakie.density!(axs, q)
    # Output
    out::type = DistributionPlots(labels, [fig])
    out
end


"""
    plot(type::Type{T}, genomes::Genomes; kwargs...)::T where {T<:ViolinPlots}

Generate a violin plot visualizing allele frequencies across different populations in genomic data.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be a subtype of ViolinPlots)
- `genomes::Genomes`: A Genomes struct containing genomic data and population information

# Keywords
- `n_loci_alleles::Int64=1_000`: Number of loci alleles to sample for plotting
- `seed::Int64=42`: Random seed for reproducibility
- `plot_size::Tuple{Int64,Int64}=(600, 450)`: Size of the output plot in pixels
- `colour_scheme::Symbol=:viridis`: Color scheme to use for the violin plots

# Returns
- `::T where {T<:ViolinPlots}`: A ViolinPlots object containing the generated plot

# Description
Creates a horizontal violin plot showing the distribution of allele frequencies for each population
in the genomic data. The plot includes population labels with sample sizes and uses different
colors for each population.

# Throws
- `ArgumentError`: If the Genomes struct dimensions are invalid or corrupted

# Example
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> genomes = extractgenomes(trials); genomes.phenotypes[1,1] = missing;

julia> vplots = GBPlots.plot(ViolinPlots, genomes);

julia> fnames = saveplots(vplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    genomes::Genomes;
    n_loci_alleles::Int64 = 1_000,
    seed::Int64 = 42,
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :viridis,
)::T where {T<:ViolinPlots}
    # type = ViolinPlots
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # n_loci_alleles = 1_000; seed = 42; plot_size = (600, 450); colour_scheme = :viridis;
    if !checkdims(genomes)
        throw(ArgumentError("Genomes struct is corrupted."))
    end
    rng::TaskLocalRNG = Random.seed!(seed)
    n, p = size(genomes.allele_frequencies)
    idx_loci_alleles::Vector{Int64} = sort(sample(rng, 1:p, minimum([n_loci_alleles, p]), replace = false))
    populations = sort(unique(genomes.populations))
    counts::Vector{Int64} = [sum(genomes.populations .== population) for population in populations]
    colours =
        repeat(colorschemes[colour_scheme][1:end], Int(ceil(length(populations) / length(colorschemes[colour_scheme]))))
    pops::Vector{Int64} = []
    q::Vector{Float64} = []
    for i = 1:n
        pop = findall(populations .== genomes.populations[i])[1]
        for j in idx_loci_alleles
            qij = genomes.allele_frequencies[i, j]
            if !ismissing(qij) && !isnan(qij) && !isinf(qij)
                append!(pops, pop)
                append!(q, qij)
            end
        end
    end
    labels = ["Violinplot of allele frequencies per population"]
    fig = CairoMakie.Figure(size = plot_size)
    axs = CairoMakie.Axis(
        fig[1, 1],
        title = labels[1],
        yticks = (1:length(populations), reverse(string.(populations, " (n=", counts, ")"))),
    )
    CairoMakie.violin!(
        axs,
        reverse(pops),
        reverse(q),
        color = colours[reverse(pops)],
        orientation = :horizontal,
        transparency = 0.75,
        scale = :width,
    )
    # CairoMakie.boxplot!(
    #     axs,
    #     reverse(pops),
    #     reverse(q),
    #     orientation = :horizontal,
    #     strokewidth = 1,
    #     outliercolor = :black,
    #     width = 1 / length(populations),
    # )
    # fig
    # Output
    out::type = ViolinPlots(labels, [fig])
    out
end

"""
    plot(type::Type{T}, genomes::Genomes; kwargs...)::T where {T<:CorHeatPlots}

Generate correlation heatmap plots for genome data.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be a subtype of CorHeatPlots)
- `genomes::Genomes`: The genome data structure to analyze

# Keywords
- `n_loci_alleles::Int64 = 1_000`: Maximum number of loci alleles to include in the plot
- `seed::Int64 = 42`: Random seed for reproducibility
- `plot_size::Tuple{Int64,Int64} = (600, 450)`: Size of each plot in pixels
- `colour_scheme::Symbol = :viridis`: Color scheme for the heatmap
- `rev_label_colors::Bool = false`: If true, reverses the text color threshold for labels
- `n_threshold_to_show_text::Int64 = 1_000`: Maximum number of cells before text labels are hidden

# Returns
- `::T`: A CorHeatPlots object containing the generated correlation heatmaps

# Description
Creates correlation heatmaps showing relationships between:
- Loci alleles across all populations
- Entries across all populations
- Loci alleles within each population
- Entries within each population

Each heatmap includes correlation values and sample sizes for each pair, with automatic
text sizing and visibility based on the number of elements being displayed.

# Throws
- `ArgumentError`: If the Genomes struct dimensions are invalid

# Example
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> hplots = GBPlots.plot(CorHeatPlots, genomes);

julia> fnames = saveplots(hplots, format="png")

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    genomes::Genomes;
    n_loci_alleles::Int64 = 1_000,
    seed::Int64 = 42,
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :viridis,
    rev_label_colors::Bool = false,
    n_threshold_to_show_text::Int64 = 1_000,
)::T where {T<:CorHeatPlots}
    # type = CorHeatPlots
    # genomes = GBCore.simulategenomes(n=100, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # seed = 42; n_loci_alleles = 1_000; plot_size = (700, 500); colour_scheme = :viridis; rev_label_colors = false; n_threshold_to_show_text = 1_000
    if !checkdims(genomes)
        throw(ArgumentError("Genomes struct is corrupted."))
    end
    rng::TaskLocalRNG = Random.seed!(seed)
    n, p = size(genomes.allele_frequencies)
    idx_loci_alleles::Vector{Int64} = sort(sample(rng, 1:p, minimum([n_loci_alleles, p]), replace = false))
    loci::Vector{String} = sort(genomes.loci_alleles[idx_loci_alleles])
    populations::Vector{String} = sort(unique(genomes.populations))
    labels = Vector{String}(undef, 2 * (1 + length(populations)))
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    # Instantiate the vectors of correlation matrices, counts per pairwise correlation, labels and grouping names
    correlations = Vector{Matrix{Float64}}(undef, length(labels))
    counts = Vector{Matrix{Int64}}(undef, length(labels))
    labellings = Vector{Vector{String}}(undef, length(labels))
    groupings = Vector{String}(undef, length(labels))
    # Population the vectors above
    loci_alleles, entries, dist =
        distances(genomes, distance_metrics = ["correlation"], idx_loci_alleles = idx_loci_alleles)
    correlations[1:2] = [dist["loci_alleles|correlation"], dist["entries|correlation"]]
    counts[1:2] = [Int.(dist["loci_alleles|counts"]), Int.(dist["entries|counts"])]
    labellings[1:2] = [loci_alleles, entries]
    groupings[1:2] = ["between loci_alleles across all populations", "between entries across all populations"]
    i = 3
    for population in populations
        idx_entries = findall(genomes.populations .== population)
        loci_alleles_per_pop, entries_per_pop, dist_per_pop = try
            distances(slice(genomes, idx_entries = idx_entries), distance_metrics = ["correlation"])
        catch
            continue
        end
        correlations[i:(i+1)] = [dist_per_pop["loci_alleles|correlation"], dist_per_pop["entries|correlation"]]
        counts[i:(i+1)] = [Int.(dist_per_pop["loci_alleles|counts"]), Int.(dist_per_pop["entries|counts"])]
        labellings[i:(i+1)] = [loci_alleles_per_pop, entries_per_pop]
        groupings[i:(i+1)] =
            [string("between loci_alleles (", population, ")"), string("between entries (", population, ")")]
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
        labels[i] = string("Correlations (genotypes) ", grp)
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
    plot(type::Type{T}, genomes::Genomes; kwargs...)::T where {T<:TreePlots}

Generate hierarchical clustering dendrograms for a Genomes object.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be subtype of TreePlots)
- `genomes::Genomes`: A Genomes struct containing genetic data

# Keyword Arguments
- `n_loci_alleles::Int64=1_000`: Number of loci-alleles to sample for plotting
- `seed::Int64=42`: Random seed for reproducibility
- `plot_size::Tuple{Int64,Int64}=(600, 450)`: Size of the output plots in pixels
- `colour_scheme::Symbol=:tol_muted`: Color scheme to use for the dendrograms
- `horizontal::Bool=true`: If true, plot dendrograms horizontally; if false, vertically

# Returns
- `TreePlots`: A TreePlots object containing two dendrograms:
  1. Clustering of loci-alleles based on genotypes
  2. Clustering of entries (samples) based on genotypes

# Details
Uses Ward's minimum variance method for hierarchical clustering based on Euclidean distances.
The dendrograms are colored using the specified color scheme and optimally ordered for visualization.

# Throws
- `ArgumentError`: If the Genomes struct dimensions are invalid

# Example
```julia
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> tplots = GBPlots.plot(TreePlots, genomes);

julia> fnames = saveplots(tplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    genomes::Genomes;
    n_loci_alleles::Int64 = 1_000,
    seed::Int64 = 42,
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :tol_muted,
    horizontal::Bool = true,
)::T where {T<:TreePlots}
    # type = TreePlots
    # genomes = GBCore.simulategenomes(n=100, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # seed = 42; n_loci_alleles = 1_000; plot_size = (700, 500); colour_scheme::Symbol = :tol_muted; horizontal= true
    if !checkdims(genomes)
        throw(ArgumentError("Genomes struct is corrupted."))
    end
    rng::TaskLocalRNG = Random.seed!(seed)
    n, p = size(genomes.allele_frequencies)
    idx_loci_alleles::Vector{Int64} = sort(sample(rng, 1:p, minimum([n_loci_alleles, p]), replace = false))
    loci::Vector{String} = sort(genomes.loci_alleles[idx_loci_alleles])
    populations::Vector{String} = sort(unique(genomes.populations))
    labels = Vector{String}(undef, 2)
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    # Estimate pairwise Euclidean distances
    loci_alleles, entries, dist =
        distances(genomes, distance_metrics = ["euclidean"], idx_loci_alleles = idx_loci_alleles)
    clusters = [
        Clustering.hclust(dist["loci_alleles|euclidean"], linkage = :ward, branchorder = :optimal),
        Clustering.hclust(dist["entries|euclidean"], linkage = :ward, branchorder = :optimal),
    ]
    leaves = [genomes.loci_alleles[idx_loci_alleles], genomes.entries]
    titles = ["Loci-alleles dendrogram (genotypes)", "Entries dendrogram (genotypes)"]
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
        # Draw the dendrogram
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
    plot(type::Type{T}, genomes::Genomes; kwargs...)::T where {T<:PCBiPlots}

Generate Principal Component Analysis (PCA) biplots for genomic data.

# Arguments
- `type::Type{T}`: The type of plot to generate (must be a subtype of PCBiPlots)
- `genomes::Genomes`: A Genomes struct containing the genomic data to plot

# Keywords
- `n_loci_alleles::Int64=1_000`: Maximum number of loci-alleles to include in the analysis
- `seed::Int64=42`: Random seed for reproducibility
- `plot_size::Tuple{Int64,Int64}=(600, 450)`: Size of the output plots in pixels
- `colour_scheme::Symbol=:tol_muted`: Color scheme to use for the plots

# Returns
- `::T where {T<:PCBiPlots}`: A PCBiPlots object containing two plots:
  1. PCA biplot of entries colored by population
  2. PCA biplot of loci colored by chromosome

# Details
The function performs the following operations:
1. Validates input dimensions
2. Samples loci-alleles if necessary
3. Processes allele frequencies matrix
4. Removes columns with zero variance
5. Performs PCA and generates two biplots:
   - One showing the relationship between entries/populations
   - One showing the relationship between loci/chromosomes

# Throws
- `ArgumentError`: If the Genomes struct is corrupted
- `ErrorException`: If the Genomes struct is too sparse

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
    genomes::Genomes;
    n_loci_alleles::Int64 = 1_000,
    seed::Int64 = 42,
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :tol_muted,
)::T where {T<:PCBiPlots}
    # type = PCBiPlots
    # genomes = GBCore.simulategenomes(n=90, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # seed = 42; n_loci_alleles = 1_000; plot_size = (700, 500); colour_scheme::Symbol = :tol_muted;
    if !checkdims(genomes)
        throw(ArgumentError("Genomes struct is corrupted."))
    end
    rng::TaskLocalRNG = Random.seed!(seed)
    n, p = size(genomes.allele_frequencies)
    idx_loci_alleles::Vector{Int64} = sort(sample(rng, 1:p, minimum([n_loci_alleles, p]), replace = false))
    genomes = slice(genomes, idx_loci_alleles = idx_loci_alleles)
    loci::Vector{String} = sort(unique(genomes.loci_alleles))
    populations::Vector{String} = sort(unique(genomes.populations))
    chromosomes, positions, alleles = loci_alleles(genomes)
    labels = Vector{String}(undef, 2)
    # Prepare the allele frequencie matrix    
    q = mean(genomes.allele_frequencies, dims = 1)[1, :]
    idx_cols = findall(.!ismissing.(q) .&& .!isnan.(q) .&& .!isinf.(q))
    if length(idx_cols) == 0
        throw(ErrorException("Genomes struct too sparse"))
    end
    G::Matrix{Float64} = genomes.allele_frequencies[:, idx_cols]
    # Remove columns with zero variance
    μ = mean(G, dims = 1)
    σ = std(G, dims = 1)
    idx = findall(σ[1, :] .> 0.0)
    G = G[:, idx]
    μ = mean(G, dims = 1)
    σ = std(G, dims = 1)
    G = (G .- μ) ./ σ
    labels = ["PCA (genotypes) biplot of entries", "PCA (genotypes) biplot of loci"]
    fig_entries = begin
        colours = [findall(populations .== pop)[1] for pop in genomes.populations]
        fig = CairoMakie.Figure(size = plot_size)
        plt = if size(G, 2) > 2
            M = fit(PCA, G)
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
            labels[1] = string(genomes.loci_alleles[idx_cols[1]], " vs ", genomes.loci_alleles[idx_cols[2]])
            x = G[:, 1]
            y = G[:, 2]
            axs = CairoMakie.Axis(
                fig[1, 1],
                title = labels[1],
                xlabel = genomes.loci_alleles[idx_cols[1]],
                ylabel = genomes.loci_alleles[idx_cols[2]],
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
    fig_loci_alleles = begin
        colours = [findall(sort(unique(chromosomes)) .== chr)[1] for chr in chromosomes]
        fig = CairoMakie.Figure(size = plot_size)
        plt = if size(G, 1) > 2
            M = fit(PCA, G')
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
                colorrange = (1, maximum([2, length(unique(chromosomes))])),
            )
        else
            labels[1] = string(genomes.loci_alleles[idx_cols[1]], " vs ", genomes.loci_alleles[idx_cols[2]])
            x = G[:, 1]
            y = G[:, 2]
            axs = CairoMakie.Axis(
                fig[1, 1],
                title = labels[1],
                xlabel = genomes.loci_alleles[idx_cols[1]],
                ylabel = genomes.loci_alleles[idx_cols[2]],
            )
            CairoMakie.scatter!(
                axs,
                x,
                y,
                color = colours,
                colormap = colour_scheme,
                colorrange = (1, maximum([2, length(unique(chromosomes))])),
            )
        end
        if length(unique(chromosomes)) > 1
            colourmap = getproperty(ColorSchemes, plt.colormap[])
            colours = colourmap[range(start = 0.0, stop = 1.0; length = length(unique(chromosomes)))]
            elems = [
                [MarkerElement(color = col, marker = :circle, markersize = 15, strokecolor = :black)] for col in colours
            ]
            CairoMakie.Legend(fig[1, 2], elems, sort(unique(chromosomes)))
        end
        fig
    end
    # Output
    out::type = PCBiPlots(labels, [fig_entries, fig_loci_alleles])
    out
end
