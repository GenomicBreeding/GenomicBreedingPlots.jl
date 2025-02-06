"""
    plot(
        type::Type{T},
        cvs::Vector{CV};
        plot_size::Tuple{Int64,Int64} = (600, 450),
        color_scheme::Symbol = :viridis,
    )::T where {T <: BarPlots}

Bar plots summarising the results of a k-fold cross-validation

# Examples
```julia
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;

julia> cvs, notes = cvperpopulation(genomes=genomes, phenomes=phenomes, models = [ridge, bayesa], n_replications = 2, n_folds = 3, verbose = false);

julia> bplots = plot(BarPlots, cvs);

julia> fnames = saveplots(bplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    cvs::Vector{CV};
    metric::String = "cor"
    plot_size::Tuple{Int64,Int64} = (600, 450),
    color_scheme::Symbol = :viridis,
)::T where {T <: BarPlots}
    # cvs::Vector{CV} = []
    # for m in 1:3
    #     for t in 1:5
    #         for p in 1:4
    #             n = 100
    #             populations = if rand() < 0.5
    #                 string.("population_", sample(1:4, n, replace=true))
    #             else
    #                 string.("population_", repeat([p], n))
    #             end
    #             entries = string.("entry_", sample(1:1_000, n, replace=true))
    #             for r in 1:5
    #                 for f in 1:5
    #                     fit = Fit(n = 10, l = 1_000); fit.model = string("model_", m); fit.populations .= string("population_", p); fit.trait = string("trait_", t); fit.metrics = Dict("cor" => rand()/maximum([1, 5*rand()]), "rmse" => rand())
    #                     cv = CV(
    #                         string("replication_", r), 
    #                         string("fold_", f), 
    #                         fit, 
    #                         populations, 
    #                         entries, 
    #                         rand(n), 
    #                         rand(n), 
    #                         fit.metrics
    #                     )
    #                     push!(cvs, cv)
    #                 end
    #             end
    #         end
    #     end
    # end
    # metric = "cor"; plot_size = (600, 450); color_scheme = :viridis;
    # Check arguments
    for (i, cv) in enumerate(cvs)
        if !checkdims(cv)
            throw(ArgumentError("The element number " * string(i) * " in the vector of CV structs is corrupted."))
        end
    end
    valid_metrics = string.(keys(cvs[1].fit.metrics))
    if !(metric ∈ valid_metrics)
        throw(ArgumentError("The supplied genomic prediction accuracy metric: `" * metric * "` does not exist in the Fit struct. Please choose from:\n\t‣ " * join(valid_metrics, "\n\t‣ ")))
    end
    # Summarise
    df_summary, df_summary_per_entry = summarise(cvs)
    # Instantiate output vectors of labels and plots
    traits = sort(unique(df_summary.trait))
    populations = sort(unique(df_summary.training_population))
    models = sort(unique(df_summary.model))
    labels = Vector{String}(undef, length(traits))
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    # Plots

    # Within population CV results
    i = 1
    for model in models
        # model = "model_3"
        idx = findall(
            (df_summary.training_population .== df_summary.validation_population) .&&
            (df_summary.model .== model)
        )
        df_sub = df_summary[idx, :]
        populations = sort(unique(df_sub.validation_population))
        traits = sort(unique(df_sub.trait))
        n = length(populations) * length(traits)
        font_size_labels = minimum([20, 20 / (0.1 * n)])
        
        fig = CairoMakie.Figure(size=plot_size);
        plt = if model !== "bulk_populations"
            x = [findall(populations .== pop)[1] for pop in df_sub.validation_population]
            y = df_sub[!, string(metric, "_mean")]
            y_std = df_sub[!, string(metric, "_std")]
            z = [findall(traits .== trait)[1] for trait in df_sub.trait]
            axs = CairoMakie.Axis(fig[1,1], xlabel = "GEBV Accuracy", limits = ((minimum([0.0, minimum(y)]), 1.5 * maximum(y)), nothing), yticks = (1:length(populations), populations), yreversed = true)
            CairoMakie.barplot!(
                axs,
                x,
                y,
                dodge = z,
                color = z,
                colorrange = (1, length(traits)),
                colormap = color_scheme,
                bar_labels = collect(1:length(y)),
                label_formatter = i -> string(round(y[i], digits=2), " (±", round(y_std[i], digits=2), ")"),
                label_size = font_size_labels,
                label = [label => (; color = i) for (i, label) in enumerate(traits)],
                direction = :x,
            )
        else
            df_sub_agg = combine(groupby(df_sub, :trait), [Symbol(string(metric, "_mean")) => mean, Symbol(string(metric, "_mean")) => std])
            x = [findall(traits .== trait)[1] for trait in df_sub_agg.trait]
            y = df_sub_agg[!, string(metric, "_mean_mean")]
            y_std = df_sub_agg[!, string(metric, "_mean_std")]
            axs = CairoMakie.Axis(fig[1,1], xlabel = "GEBV Accuracy", limits = ((minimum([0.0, minimum(y)]), 1.5 * maximum(y)), nothing), yticks = (1:length(traits), traits), yreversed = true)
            CairoMakie.barplot!(
                axs,
                x,
                y,
                color = x,
                colorrange = (1, length(traits)),
                colormap = color_scheme,
                bar_labels = collect(1:length(y)),
                label_formatter = i -> string(round(y[i], digits=2), " (±", round(y_std[i], digits=2), ")"),
                label_size = font_size_labels,
                label = [label => (; color = i) for (i, label) in enumerate(traits)],
                direction = :x,
            )
        end
        CairoMakie.Legend(fig[1,2], axs)
        plots[i] = fig
        i += 1
    end
end




"""
    plot(
        type::Type{T},
        cvs::Vector{CV};
        plot_size::Tuple{Int64,Int64} = (600, 450),
        color_scheme::Symbol = :viridis,
    )::T where {T <: BoxPlots}

Bar plots summarising the results of a k-fold cross-validation

# Examples
```julia
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials); phenomes.phenotypes[1,1] = missing;

julia> cvs, notes = cvperpopulation(genomes=genomes, phenomes=phenomes, models = [ridge, bayesa], n_replications = 2, n_folds = 3, verbose = false);

julia> bplots = plot(BoxPlots, cvs);

julia> fnames = saveplots(bplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    cvs::Vector{CV};
    plot_size::Tuple{Int64,Int64} = (600, 450),
    color_scheme::Symbol = :viridis,
)::T where {T <: BoxPlots}
    # cvs::Vector{CV} = []
    # for m in 1:3
    #     for t in 1:5
    #         for p in 1:4
    #             n = 100
    #             populations = if rand() < 0.5
    #                 string.("population_", sample(1:4, n, replace=true))
    #             else
    #                 string.("population_", repeat([p], n))
    #             end
    #             entries = string.("entry_", sample(1:1_000, n, replace=true))
    #             for r in 1:5
    #                 for f in 1:5
    #                     fit = Fit(n = 10, l = 1_000); fit.model = string("model_", m); fit.populations .= string("population_", p); fit.trait = string("trait_", t); fit.metrics = Dict("cor" => rand()/maximum([1, 5*rand()]), "rmse" => rand())
    #                     cv = CV(
    #                         string("replication_", r), 
    #                         string("fold_", f), 
    #                         fit, 
    #                         populations, 
    #                         entries, 
    #                         rand(n), 
    #                         rand(n), 
    #                         fit.metrics
    #                     )
    #                     push!(cvs, cv)
    #                 end
    #             end
    #         end
    #     end
    # end
    # plot_size = (600, 450); color_scheme = :viridis;


    # Check arguments
    for (i, cv) in enumerate(cvs)
        if !checkdims(cv)
            throw(ArgumentError("The element number " * string(i) * " in the vector of CV structs is corrupted."))
        end
    end
    # Tabularise
    df_metrics, df_ys = tabularise(cvs)
    # Instantiate output vectors of labels and plots
    traits = sort(unique(df_metrics.trait))
    populations = sort(unique(df_metrics.training_population))
    labels = Vector{String}(undef, length(traits))
    plots = Vector{CairoMakie.Figure}(undef, length(labels))
    # Plots

    # Within population CV results
    model = "model_1"
    idx = findall(
        (df_metrics.training_population .== df_metrics.validation_population) .&&
        (df_metrics.model .== model)
    )
    df_sub = df_metrics[idx, :]
    populations = sort(unique(df_sub.validation_population))
    traits = sort(unique(df_sub.trait))
    df_sub.numeric_validation_population = [findall(populations .== pop)[1] for pop in df_sub.validation_population]
    df_sub.numeric_trait = [findall(traits .== trait)[1] for trait in df_sub.trait]
    fig = CairoMakie.Figure(size=plot_size);
    axs = CairoMakie.Axis(fig[1,1], xticks = (1:length(populations), populations))
    plt = CairoMakie.boxplot!(
        axs,
        df_sub.numeric_validation_population,
        df_sub.cor,
        dodge = df_sub.numeric_trait,
        color = df_sub.numeric_trait,
        colorrange = (1, length(traits)),
        colormap = color_scheme,
        label = [label => (; color = i) for (i, label) in enumerate(traits)]
    )
    colourmap = getproperty(ColorSchemes, plt.colormap[])
	colours = colourmap[range(start=0.0, stop=1.0; length=length(traits))]
	elems = [[MarkerElement(color = col, marker=:circle, markersize = 15, strokecolor = :black)] for col in colours]
	CairoMakie.Legend(fig[1,2], elems, traits)
    fig
end