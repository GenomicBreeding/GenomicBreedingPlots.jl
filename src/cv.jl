"""
    plot(
        type::Type{T},
        cvs::Vector{CV};
        metric::String = "cor",
        plot_size::Tuple{Int64,Int64} = (600, 450),
        colour_scheme::Symbol = :viridis,
    )::T where {T <: BarPlots}

Bar plots summarising the results of a k-fold cross-validation

# Examples
```julia
julia> cvs::Vector{CV} = [];

julia> for m in 1:3
            for t in 1:5
                for p in 1:4
                    n = 100
                    populations = if rand() < 0.5
                        string.("population_", sample(1:4, n, replace=true))
                    else
                        string.("population_", repeat([p], n))
                    end
                    entries = string.("entry_", sample(1:1_000, n, replace=true))
                    for r in 1:5
                        for f in 1:5
                            fit = Fit(n = 10, l = 1_000); fit.model = string("model_", m); fit.populations .= string("population_", p); fit.trait = string("trait_", t); fit.metrics = Dict("cor" => rand()/maximum([1, 5*rand()]), "rmse" => rand())
                            cv = CV(
                                string("replication_", r), 
                                string("fold_", f), 
                                fit, 
                                populations, 
                                entries, 
                                rand(n), 
                                rand(n), 
                                fit.metrics
                            )
                            push!(cvs, cv)
                        end
                    end
                end
            end
        end;

julia> bplots = GBPlots.plot(BarPlots, cvs);

julia> fnames = saveplots(bplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    cvs::Vector{CV};
    metric::String = "cor",
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :viridis,
)::T where {T<:BarPlots}
    # type = BarPlots
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
    # metric = "cor"; plot_size = (600, 450); colour_scheme = :viridis;
    # Check arguments
    for (i, cv) in enumerate(cvs)
        if !checkdims(cv)
            throw(ArgumentError("The element number " * string(i) * " in the vector of CV structs is corrupted."))
        end
    end
    valid_metrics = string.(keys(cvs[1].fit.metrics))
    if !(metric ∈ valid_metrics)
        throw(
            ArgumentError(
                "The supplied genomic prediction accuracy metric: `" *
                metric *
                "` does not exist in the Fit struct. Please choose from:\n\t‣ " *
                join(valid_metrics, "\n\t‣ "),
            ),
        )
    end
    # Tabularise, i.e. extract the raw metrics per trait, model, population, replication and fold
    df_metrics_all, _ = tabularise(cvs)
    # Instantiate output vectors of labels and plots
    labels::Vector{String} = []
    plots::Vector{CairoMakie.Figure} = []
    # Plots
    # i = 1
    for cv_type in ["Within population", "Across populations"]
        # cv_type = "Across populations"
        df_metrics, x_names, z_names = if cv_type == "Within population"
            # Within population CV results
            idx = findall(
                isnothing.(match.(Regex(";"), df_metrics_all.training_population)) .&&
                df_metrics_all.training_population .== df_metrics_all.validation_population,
            )
            df_metrics = df_metrics_all[idx, :]
            x_names = ["model", "trait", "validation_population"]
            z_names = vcat("bulk", x_names)
            (df_metrics, x_names, z_names)
        else
            # Across population CV results
            idx = findall(df_metrics_all.training_population .!= df_metrics_all.validation_population)
            df_metrics = df_metrics_all[idx, :]
            x_names = ["validation_population"]
            z_names = ["model", "trait"]
            (df_metrics, x_names, z_names)
        end
        for x_name in x_names
            # x_name = x_names[1]
            for z_name in z_names
                # z_name = z_names[1]
                if x_name == z_name
                    continue
                end
                w_name, w_levels = if cv_type == "Within population"
                    # Within population CV
                    w_name = x_names[findall([sum([x_name, z_name] .== a) == 0 for a in x_names])[1]]
                    w_levels = if z_name == "bulk"
                        [nothing]
                    else
                        sort(unique(df_metrics[!, Symbol(w_name)]))
                    end
                    (w_name, w_levels)
                else
                    # Across populations CV
                    w_name = ("training_population", z_names[z_names.!=z_name][1])
                    w_levels = []
                    for pop in sort(unique(df_metrics[!, Symbol(w_name[1])]))
                        for trait_or_model in sort(unique(df_metrics[!, Symbol(w_name[2])]))
                            push!(w_levels, (pop, trait_or_model))
                        end
                    end
                    (w_name, w_levels)
                end
                for w_level in w_levels
                    # w_level = w_levels[6]
                    # w_level = ("population_1", "trait_1")
                    # println(string("cv_type = \"", cv_type, "\"; x_name = \"", x_name, "\"; z_name = \"", z_name, "\"; w_level = \"", w_level, "\""))
                    df_metrics_sub = if isnothing(w_level)
                        df_metrics
                    else
                        idx = if typeof(w_name) == Tuple{String,String}
                            findall((df_metrics[:, w_name[1]] .== w_level[1]) .&& (df_metrics[:, w_name[2]] .== w_level[2]))
                        else
                            findall(df_metrics[:, w_name] .== w_level)
                        end
                        if length(idx) == 0
                            continue
                        end
                        df_metrics_sub = df_metrics[idx, :]
                    end
                    df, z_levels, z, colours = if z_name == "bulk"
                        df = combine(
                            groupby(df_metrics_sub, Symbol(x_name)),
                            [Symbol(metric) => mean, Symbol(metric) => std],
                        )
                        if nrow(df) == 0
                            continue
                        end
                        x_levels = sort(unique(df[!, Symbol(x_name)]))
                        x = [findall(x_levels .== a)[1] for a in df[!, x_name]]
                        z_levels = 1
                        z = 1
                        colours = x
                        (df, z_levels, z, colours)
                    else
                        df = combine(
                            groupby(df_metrics_sub, [Symbol(x_name), Symbol(z_name)]),
                            [Symbol(metric) => mean, Symbol(metric) => std],
                        )
                        if nrow(df) == 0
                            continue
                        end
                        x_levels = sort(unique(df[!, Symbol(x_name)]))
                        x = [findall(x_levels .== a)[1] for a in df[!, x_name]]
                        z_levels = sort(unique(df[!, Symbol(z_name)]))
                        z, colours = if length(z_levels) > 1
                            z = [findall(z_levels .== a)[1] for a in df[!, z_name]]
                            z, z
                        else
                            1, x
                        end
                        (df, z_levels, z, colours)
                    end
                    # Rename level combinations to fit in the plot area
                    x_levels = begin
                        x_levels_renamed = deepcopy(x_levels)
                        for (i, x_level) in enumerate(x_levels)
                            if !isnothing(match(Regex(";"), x_level))
                                x_levels_renamed[i] = string("Bulk[", replace.(x_level, ";" => "\n"), "]")
                            end
                        end
                        x_levels_renamed
                    end
                    y = df[!, string(metric, "_mean")]
                    y_std = df[!, string(metric, "_std")]
                    n = length(x_levels) * length(z_levels)
                    font_size_labels = minimum([14, 20 / (0.1 * n)])
                    title, label = if isnothing(w_level)
                        across_names = x_names[[sum([x_name, z_name] .== a) == 0 for a in x_names]]
                        title = string("Across ", across_names[1], " and ", across_names[2])
                        label =
                            string(cv_type, " | x: ", x_name, " | y: ", metric, " | z: ", z_name, " | subset: none")
                        (title, label)
                    else
                        (
                            if typeof(w_name) == Tuple{String,String}
                                title = join(string.(w_name, ": ", w_level), " | ")
                                label = string(
                                    cv_type,
                                    " | x: ",
                                    x_name,
                                    " | y: ",
                                    metric,
                                    " | z: ",
                                    z_name,
                                    " | subset: ",
                                    join(w_level, "&"),
                                )
                                (title, label)
                            else
                                title = string(w_name, ": ", w_level)
                                label = string(
                                    cv_type,
                                    " | x: ",
                                    x_name,
                                    " | y: ",
                                    metric,
                                    " | z: ",
                                    z_name,
                                    " | subset: ",
                                    w_level,
                                )
                                (title, label)
                            end
                        )
                    end
                    # Plot
                    fig = CairoMakie.Figure(size = plot_size)
                    axs = CairoMakie.Axis(
                        fig[1, 1],
                        title = title,
                        xlabel = "GEBV Accuracy",
                        ylabel = x_name,
                        limits = ((minimum([0.0, minimum(y)]), 1.5 * maximum(y)), nothing),
                        yticks = (1:length(x_levels), x_levels),
                        yreversed = true,
                    )
                    plt = CairoMakie.barplot!(
                        axs,
                        x,
                        y,
                        dodge = z,
                        color = colours,
                        colorrange = (1, maximum([2, length(unique(colours))])),
                        colormap = colour_scheme,
                        bar_labels = collect(1:length(y)),
                        label_formatter = i ->
                            !isnan(y_std[i]) ?
                            string(round(y[i], digits = 2), " (±", round(y_std[i], digits = 2), ")") :
                            string(round(y[i], digits = 2)),
                        label_size = font_size_labels,
                        label = [label => (; color = i) for (i, label) in enumerate(z_levels)],
                        direction = :x,
                    )
                    if length(z) > 1
                        CairoMakie.Legend(fig[1, 2], axs)
                    end
                    # fig
                    # println(i); i += 1
                    # println(string("x_name = \"", x_name, "\"; z_name = \"", z_name, "\"; w_level = \"", w_level, "\""))
                    push!(labels, label)
                    push!(plots, fig)
                end
            end
        end
    end
    # Output
    out::type = BarPlots(labels, plots)
    out

end

"""
    plot(
        type::Type{T},
        cvs::Vector{CV};
        metric::String = "cor",
        plot_size::Tuple{Int64,Int64} = (600, 450),
        colour_scheme::Symbol = :viridis,
    )::T where {T <: BoxPlots}

Box plots summarising the results of a k-fold cross-validation

# Examples
```julia
julia> cvs::Vector{CV} = [];

julia> for m in 1:3
            for t in 1:5
                for p in 1:4
                    n = 100
                    populations = if rand() < 0.5
                        string.("population_", sample(1:4, n, replace=true))
                    else
                        string.("population_", repeat([p], n))
                    end
                    entries = string.("entry_", sample(1:1_000, n, replace=true))
                    for r in 1:5
                        for f in 1:5
                            fit = Fit(n = 10, l = 1_000); fit.model = string("model_", m); fit.populations .= string("population_", p); fit.trait = string("trait_", t); fit.metrics = Dict("cor" => rand()/maximum([1, 5*rand()]), "rmse" => rand())
                            cv = CV(
                                string("replication_", r), 
                                string("fold_", f), 
                                fit, 
                                populations, 
                                entries, 
                                rand(n), 
                                rand(n), 
                                fit.metrics
                            )
                            push!(cvs, cv)
                        end
                    end
                end
            end
        end;

julia> bplots = GBPlots.plot(BoxPlots, cvs);

julia> fnames = saveplots(bplots)

julia> rm.(fnames);
```
"""
function plot(
    type::Type{T},
    cvs::Vector{CV};
    metric::String = "cor",
    plot_size::Tuple{Int64,Int64} = (600, 450),
    colour_scheme::Symbol = :viridis,
)::T where {T<:BoxPlots}
    # type = BoxPlots
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
    # metric = "cor"; plot_size = (600, 450); colour_scheme = :viridis;
    # Check arguments
    for (i, cv) in enumerate(cvs)
        if !checkdims(cv)
            throw(ArgumentError("The element number " * string(i) * " in the vector of CV structs is corrupted."))
        end
    end
    valid_metrics = string.(keys(cvs[1].fit.metrics))
    if !(metric ∈ valid_metrics)
        throw(
            ArgumentError(
                "The supplied genomic prediction accuracy metric: `" *
                metric *
                "` does not exist in the Fit struct. Please choose from:\n\t‣ " *
                join(valid_metrics, "\n\t‣ "),
            ),
        )
    end
    # Tabularise, i.e. extract the raw metrics per trait, model, population, replication and fold
    df_metrics_all, _ = tabularise(cvs)
    # Instantiate output vectors of labels and plots
    labels::Vector{String} = []
    plots::Vector{CairoMakie.Figure} = []
    # Plots
    # i = 1
    for cv_type in ["Within population", "Across populations"]
        # cv_type = "Across populations"
        df_metrics, x_names, z_names = if cv_type == "Within population"
            # Within population CV results
            idx = findall(
                isnothing.(match.(Regex(";"), df_metrics_all.training_population)) .&&
                df_metrics_all.training_population .== df_metrics_all.validation_population,
            )
            df_metrics = df_metrics_all[idx, :]
            x_names = ["model", "trait", "validation_population"]
            z_names = vcat("bulk", x_names)
            (df_metrics, x_names, z_names)
        else
            # Across population CV results
            idx = findall(df_metrics_all.training_population .!= df_metrics_all.validation_population)
            df_metrics = df_metrics_all[idx, :]
            x_names = ["validation_population"]
            z_names = ["model", "trait"]
            (df_metrics, x_names, z_names)
        end
        for x_name in x_names
            # x_name = x_names[1]
            for z_name in z_names
                # z_name = z_names[1]
                if x_name == z_name
                    continue
                end
                w_name, w_levels = if cv_type == "Within population"
                    # Within population CV
                    w_name = x_names[findall([sum([x_name, z_name] .== a) == 0 for a in x_names])[1]]
                    w_levels = if z_name == "bulk"
                        [nothing]
                    else
                        sort(unique(df_metrics[!, Symbol(w_name)]))
                    end
                    (w_name, w_levels)
                else
                    # Across populations CV
                    w_name = ("training_population", z_names[z_names.!=z_name][1])
                    w_levels = []
                    for pop in sort(unique(df_metrics[!, Symbol(w_name[1])]))
                        for trait_or_model in sort(unique(df_metrics[!, Symbol(w_name[2])]))
                            push!(w_levels, (pop, trait_or_model))
                        end
                    end
                    (w_name, w_levels)
                end
                for w_level in w_levels
                    # w_level = ("population_1", "model_1")
                    # println(i); i += 1
                    # println(string("cv_type = \"", cv_type, "\"; x_name = \"", x_name, "\"; z_name = \"", z_name, "\"; w_level = \"", w_level, "\""))
                    df_metrics_sub = if isnothing(w_level)
                        df_metrics
                    else
                        idx = if typeof(w_name) == Tuple{String,String}
                            findall((df_metrics[:, w_name[1]] .== w_level[1]) .&& (df_metrics[:, w_name[2]] .== w_level[2]))
                        else
                            findall(df_metrics[:, w_name] .== w_level)
                        end
                        if length(idx) == 0
                            continue
                        end
                        df_metrics_sub = df_metrics[idx, :]
                    end
                    if z_name == "bulk"
                        x_levels = sort(unique(df_metrics[!, Symbol(x_name)]))
                        df_metrics_sub.__x__ = [findall(x_levels .== a)[1] for a in df_metrics[!, x_name]]
                        z_levels = 1
                        df_metrics_sub.__z__ .= 1
                        df_metrics_sub.__colours__ = df_metrics_sub.__x__
                    else
                        x_levels = sort(unique(df_metrics_sub[!, Symbol(x_name)]))
                        df_metrics_sub.__x__ = [findall(x_levels .== a)[1] for a in df_metrics_sub[!, x_name]]
                        z_levels = sort(unique(df_metrics_sub[!, Symbol(z_name)]))
                        df_metrics_sub.__z__ = [findall(z_levels .== a)[1] for a in df_metrics_sub[!, z_name]]
                        df_metrics_sub.__colours__ = df_metrics_sub.__z__
                    end
                    # Rename level combinations to fit in th plot area
                    x_levels = begin
                        x_levels_renamed = deepcopy(x_levels)
                        for (i, x_level) in enumerate(x_levels)
                            if !isnothing(match(Regex(";"), x_level))
                                x_levels_renamed[i] = string("Bulk[", replace.(x_level, ";" => "\n"), "]")
                            end
                        end
                        x_levels_renamed
                    end
                    df_metrics_sub.__y__ = df_metrics_sub[!, metric]
                    n = length(x_levels) * length(z_levels)
                    font_size_labels = minimum([14, 20 / (0.1 * n)])
                    title, label = if isnothing(w_level)
                        across_names = x_names[[sum([x_name, z_name] .== a) == 0 for a in x_names]]
                        title = string("Across ", across_names[1], " and ", across_names[2])
                        label =
                            string(cv_type, " | x: ", x_name, " | y: ", metric, " | z: ", z_name, " | subset: none")
                        (title, label)
                    else
                        (
                            if typeof(w_name) == Tuple{String,String}
                                title = join(string.(w_name, ": ", w_level), " | ")
                                label = string(
                                    cv_type,
                                    " | x: ",
                                    x_name,
                                    " | y: ",
                                    metric,
                                    " | z: ",
                                    z_name,
                                    " | subset: ",
                                    join(w_level, "&"),
                                )
                                (title, label)
                            else
                                title = string(w_name, ": ", w_level)
                                label = string(
                                    cv_type,
                                    " | x: ",
                                    x_name,
                                    " | y: ",
                                    metric,
                                    " | z: ",
                                    z_name,
                                    " | subset: ",
                                    w_level,
                                )
                                (title, label)
                            end
                        )
                    end
                    # Plot
                    fig = CairoMakie.Figure(size = plot_size)
                    axs = CairoMakie.Axis(
                        fig[1, 1],
                        title = title,
                        xlabel = x_name,
                        ylabel = "GEBV Accuracy",
                        limits = (
                            nothing,
                            (minimum([0.0, minimum(df_metrics_sub.__y__)]), 1.5 * maximum(df_metrics_sub.__y__)),
                        ),
                        xticks = (1:length(x_levels), x_levels),
                    )
                    plt = if length(z_levels) == 1
                        CairoMakie.boxplot!(axs, df_metrics_sub.__x__, df_metrics_sub.__y__, dodge = df_metrics_sub.__z__)
                    else
                        plt = CairoMakie.boxplot!(
                            axs,
                            df_metrics_sub.__x__,
                            df_metrics_sub.__y__,
                            dodge = df_metrics_sub.__z__,
                            color = df_metrics_sub.__colours__,
                            colorrange = (1, maximum([2, length(unique(df_metrics_sub.__colours__))])),
                            colormap = colour_scheme,
                            label = [label => (; color = i) for (i, label) in enumerate(z_levels)],
                        )
                        colourmap = getproperty(ColorSchemes, plt.colormap[])
                        colours = colourmap[range(
                            start = 0.0,
                            stop = 1.0;
                            length = length(unique(df_metrics_sub.__colours__)),
                        )]
                        elems = [
                            [MarkerElement(color = col, marker = :circle, markersize = 15, strokecolor = :black)] for col in colours
                        ]
                        CairoMakie.Legend(fig[1, 2], elems, z_levels)
                        plt
                    end
                    fig
                    push!(labels, label)
                    push!(plots, fig)
                end
            end
        end
    end
    # Output
    out::type = BoxPlots(labels, plots)
    out
end
