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
    # type = ViolinPlots; phenomes = Phenomes(n=100, t=3); phenomes.entries = string.("entry_", 1:100); phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace=true); phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(Distributions.MvNormal([1,2,3], LinearAlgebra.I), 100)'; phenomes.phenotypes[1, 1] = missing;
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
    idx = findall(
        .!ismissing.(df_reshaped[!, :value]) .&& .!isnan.(df_reshaped[!, :value]) .&& .!isinf.(df_reshaped[!, :value]),
    )
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
        permute = (:x, :y),
    )
    StatsPlots.@df df_reshaped[idx, :] StatsPlots.boxplot!(:variable, :value, fill = (0.5), permute = (:x, :y))
    # StatsPlots.@df df_reshaped[idx, :] StatsPlots.dotplot!(:variable, :value, fill=(0.5)); gui(p)
    plots[i] = p
    # Per trait
    for trait in Symbol.(phenomes.traits)
        # trait = Symbol.(phenomes.traits)[1]
        i += 1
        println(i)
        df.pop = deepcopy(df.populations)
        idx = findall(.!ismissing.(df[!, trait]) .&& .!isnan.(df[!, trait]) .&& .!isinf.(df[!, trait]))
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










"""

    plotinteractive2d(
        phenomes::Phenomes; 
        idx_entries::Union{Nothing, Vector{Int64}} = nothing,
        idx_traits::Union{Nothing, Vector{Int64}} = nothing,
        ditch_some_entries_to_keep_all_traits::Bool = true,
    )::Figure

Interactive biplot with histogram and correlation heatmap

# Examples
```julia
phenomes = Phenomes(n = 100, t = 3)
phenomes.entries = string.("entry_", 1:100)
phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace = true)
phenomes.traits = ["trait_1", "trait_2", "long_trait_name number 3"]
phenomes.phenotypes = rand(Distributions.MvNormal([1, 2, 3], LinearAlgebra.I), 100)'
phenomes.phenotypes[1, 1] = missing
f = plotinteractive2d(phenomes)
```
"""
function GBPlots.plotinteractive2d(
    phenomes::Phenomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    idx_traits::Union{Nothing,Vector{Int64}} = nothing,
    ditch_some_entries_to_keep_all_traits::Bool = true,
)::Figure
    # phenomes = Phenomes(n = 100, t = 3)
    # phenomes.entries = string.("entry_", 1:100)
    # phenomes.populations = StatsBase.sample(string.("pop_", 1:5), 100, replace = true)
    # phenomes.traits = ["trait_1", "trait_2", "long_trait_name number 3"]
    # phenomes.phenotypes = rand(Distributions.MvNormal([1, 2, 3], LinearAlgebra.I), 100)'
    # phenomes.phenotypes[1, 1] = missing
    # idx_entries = nothing; idx_traits = nothing; ditch_some_entries_to_keep_all_traits = true
    # Check arguments
    if !checkdims(phenomes)
        throw(ArgumentError("The phenomes struct is corrupted."))
    end
    if isnothing(idx_entries)
        idx_entries = collect(1:length(phenomes.entries))
    else
        if (minimum(idx_entries) < 1) .|| maximum(idx_entries) > length(phenomes.entries)
            throw(
                ArgumentError(
                    "The indexes of the entries, `idx_entries` are out of bounds. Expected range: from 1 to " *
                    string(length(phenomes.entries)) *
                    " while the supplied range is from " *
                    string(minimum(idx_entries)) *
                    " to " *
                    string(maximum(idx_entries)) *
                    ".",
                ),
            )
        end
    end
    if isnothing(idx_traits)
        idx_traits = collect(1:length(phenomes.traits))
    else
        if (minimum(idx_traits) < 1) .|| maximum(idx_traits) > length(phenomes.traits)
            throw(
                ArgumentError(
                    "The indexes of the traits, `idx_traits` are out of bounds. Expected range: from 1 to " *
                    string(length(phenomes.traits)) *
                    " while the supplied range is from " *
                    string(minimum(idx_traits)) *
                    " to " *
                    string(maximum(idx_traits)) *
                    ".",
                ),
            )
        end
    end
    # Extract phenotypes and tabularise
    phenomes = slice(phenomes, idx_entries = idx_entries, idx_traits = idx_traits)
    df = tabularise(phenomes)
    traits = sort(unique(phenomes.traits))
    populations = sort(unique(phenomes.populations))
    # Remove entries with at least 1 missing/NaN/Inf trait or remove trait/s with missing/NaN/Inf to keep all traits
    idx_rows, idx_cols = if ditch_some_entries_to_keep_all_traits
        idx_rows = findall(
            mean(Matrix(.!ismissing.(df[:, 4:end]) .&& .!isnan.(df[:, 4:end]) .&& .!isinf.(df[:, 4:end])), dims = 2)[
                :,
                1,
            ] .== 1,
        )
        idx_cols = collect(1:size(df, 2))
        if length(idx_rows) == 0
            throw(
                ArgumentError(
                    "All entries have at least 1 missing trait value. Please consider setting `ditch_some_entries_to_keep_all_traits` to false or imputing missing phenotypes.",
                ),
            )
        end
        idx_rows, idx_cols
    else
        idx_rows = collect(1:size(df, 1))
        idx_cols = findall(
            mean(Matrix(.!ismissing.(df[:, 4:end]) .&& .!isnan.(df[:, 4:end]) .&& .!isinf.(df[:, 4:end])), dims = 1)[
                1,
                :,
            ] .== 1,
        )
        if length(idx_cols) == 0
            throw(
                ArgumentError(
                    "All traits have at least 1 entry with missing data. Please consider setting `ditch_some_entries_to_keep_all_traits` to true or imputing missing phenotypes.",
                ),
            )
        end
        idx_rows, idx_cols
    end
    df = df[idx_rows, idx_cols]
    # Extract the first 2 principal components
    traits = names(df)[4:end]
    t = length(traits)
    if t <= 2
        throw(ArgumentError("No need to do PCA as you only have at most 2 traits."))
    end
    A::Matrix{Float64} = Matrix(df[:, 4:end])
    A = (A .- mean(A, dims = 1)) ./ std(A, dims = 1)
    # Remove traits with no variation
    v = var(A, dims = 1)[1, :]
    idx_cols = findall((abs.(v .- 1) .< 0.00001) .&& .!isnan.(v) .&& .!ismissing.(v) .&& .!isinf.(v))
    A = A[:, idx_cols]
    M = fit(PCA, A; maxoutdim = 2)
    vec_missing::Vector{Union{Missing,Float64}} = repeat([missing], size(df, 1))
    df.pc1 = deepcopy(vec_missing)
    df.pc2 = deepcopy(vec_missing)
    df[!, :pc1] = M.proj[:, 1]
    df[!, :pc2] = M.proj[:, 2]
    # Include the 2 PCs in the list of traits
    traits = names(df)[4:end]
    traits = traits[sortperm(uppercase.(traits))]
    # Activate Makie using the OpenGL plotting backend
    GLMakie.activate!()
    GLMakie.closeall() # close any open screen
    # Set plot size
    height, width = 1_200, 800
    # Define the entire plot/figure
    fig = begin
        fig = Figure(size = (height, width))
        # Set the trait 1 selector
        default_trait_x = if length(traits) > 2
            "pc1"
        else
            traits[1]
        end
        menu_trait_x = Menu(fig, options = traits, default = default_trait_x)
        tb_x = Textbox(fig)
        # Set the trait 2 selector
        default_trait_y = if length(traits) > 2
            "pc2"
        else
            traits[2]
        end
        menu_trait_y = Menu(fig, options = traits, default = default_trait_y)
        tb_y = Textbox(fig)
        # Instantiate Pearson's correlation value per pair of traits
        ρ = Observable(string(round(cor(df[!, default_trait_x], df[!, default_trait_y]), digits = 4)))
        # Place these trait menus in a left sidebar
        fig[1, 1] = vgrid!(
            Label(fig, "Search trait 1", width = nothing),
            tb_x,
            menu_trait_x,
            Label(fig, "Search trait 2", width = nothing),
            tb_y,
            menu_trait_y,
            Label(fig, @lift("ρ = $($ρ)"));
            tellheight = false,
        )
        # Place the main scatter plot with histograms for the 2 traits
        fig_main = fig[1:2, 2] = GridLayout()
        plot_hist_x = Axis(fig_main[1, 1])
        plot_scatter = Axis(fig_main[2, 1], xlabel = default_trait_x, ylabel = default_trait_y)
        plot_hist_y = Axis(fig_main[2, 2])

        plot_heatmap = Axis(
            fig[2, 1],
            title = "Trait correlations",
            xticks = (1:length(traits), traits),
            yticks = (1:length(traits), reverse(traits)),
            xaxisposition = :top,
            xticklabelrotation = deg2rad(90),
        )

        C = cor(Matrix(df[!, traits]))
        C = C[:, reverse(collect(1:end))]
        GLMakie.heatmap!(
            plot_heatmap,
            1:length(traits),
            1:length(traits),
            C,
            colorrange = (-1, 1),
            inspector_label = (self, (i, j), p) ->
                string("ρ = ", round(C[i, j], digits = 4), "\n", traits[i], "\n", reverse(traits)[j]),
        )


        X = Observable{Any}(df[!, default_trait_x])
        Y = Observable{Any}(df[!, default_trait_y])

        for pop in populations
            idx = findall(df.populations .== pop)

            x = @lift($X[idx])
            y = @lift($Y[idx])

            GLMakie.hist!(plot_hist_x, x)
            GLMakie.hist!(plot_hist_y, y, direction = :x)
            GLMakie.scatter!(
                plot_scatter,
                x,
                y,
                label = pop,
                inspector_label = (self, i, p) -> string(df.entries[idx][i], "\n(", df.populations[idx][i], ")"),
            )
        end
        connect!(ρ, @lift(rpad(round(cor($X, $Y), digits = 4), 2 * maximum(length.(traits)), " ")))

        leg = Legend(fig_main[1, 2], plot_scatter)
        GLMakie.hidedecorations!(plot_hist_x, grid = false)
        GLMakie.hidedecorations!(plot_hist_y, grid = false)
        leg.tellheight = true


        # Reactivity stuff
        # Search-ish functionality-ish
        on(tb_x.stored_string) do s
            bool_matches = .!isnothing.(match.(Regex(s, "i"), traits)) # case-insensitive with the Regex flag "i"
            if sum(bool_matches) > 0
                traits_reordered = vcat(traits[bool_matches], traits[.!bool_matches])
                menu_trait_x.is_open = true
                menu_trait_x.selection = traits_reordered[1]
                menu_trait_x.i_selected = 1
                menu_trait_x.options = traits_reordered
            else
                menu_trait_x.options = traits
            end
        end
        on(tb_y.stored_string) do s
            bool_matches = .!isnothing.(match.(Regex(s, "i"), traits)) # case-insensitive with the Regex flag "i"
            if sum(bool_matches) > 0
                traits_reordered = vcat(traits[bool_matches], traits[.!bool_matches])
                menu_trait_y.is_open = true
                menu_trait_y.selection = traits_reordered[1]
                menu_trait_y.i_selected = 1
                menu_trait_y.options = traits_reordered
            else
                menu_trait_y.options = traits
            end
        end
        # Drop down menus
        on(menu_trait_x.selection) do s
            X[] = df[!, s]
            plot_scatter.xlabel = s
            autolimits!(plot_scatter)
            autolimits!(plot_hist_x)
            autolimits!(plot_hist_y)
        end

        on(menu_trait_y.selection) do s
            # trait_y[] = s
            Y[] = df[!, s]
            plot_scatter.ylabel = s
            autolimits!(plot_scatter)
            autolimits!(plot_hist_x)
            autolimits!(plot_hist_y)
        end

        DataInspector(fig)
        fig
    end
    # Output
    fig
end
