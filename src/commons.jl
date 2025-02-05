abstract type PlotsGB end

mutable struct DistributionPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

mutable struct ViolinPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

mutable struct CorHeatPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

mutable struct TreePlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

mutable struct BarPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

function GBCore.checkdims(x::PlotsGB)::Bool
    length(x.labels) == length(x.plots)
end

function labeltofname(; label::String, prefix::String, suffix::String)::String
    label = replace(label, " " => "_")
    label = replace(label, "\n" => "_")
    label = replace(label, "\t" => "_")
    label = replace(label, "(" => "_")
    label = replace(label, ")" => "_")
    label = replace(label, ":" => "_")
    label = replace(label, "=" => "_")
    fname = join([prefix, label], "-") * "." * suffix
    fname = replace(fname, ".." => ".")
    fname = replace(fname, "__" => "_")
    fname = replace(fname, "--" => "-")
    fname = replace(fname, "_." => ".")
    fname = replace(fname, "._" => "_")
    fname = replace(fname, "-." => ".")
    fname = replace(fname, ".-" => "-")
    fname
end

function saveplots(
    plots::PlotsGB;
    idx::Vector{Int64} = [0],
    format::String = "svg",
    prefix::String = "",
    use_labels::Bool = true,
)::Vector{String}
    # phenomes = Phenomes(n=10, t=3); phenomes.entries = string.("entry_", 1:10); phenomes.populations .= "pop_1"; phenomes.traits = ["A", "B", "C"]; phenomes.phenotypes = rand(10,3);
    # plots = GBPlots.plot(DistributionPlots, phenomes); idx = [1, 3, 5]; format = "svg"; prefix = ""; use_labels = false;
    # Check arguments
    if !checkdims(plots)
        throw(ArgumentError("The plots::DistributionPlots is corrupted."))
    end
    if idx == [0]
        idx = collect(1:length(plots.labels))
    end
    if (maximum(idx) > length(plots.labels)) || (minimum(idx) < 1)
        throw(ArgumentError("The provided indexes (idx) are out of bounds."))
    end
    if (format != "svg") && (format != "png") && (format != "pdf")
        throw(ArgumentError("we accept the following file formats: 'svg', 'png' and 'pdf'."))
    end
    if prefix == ""
        prefix = string(typeof(plots))
    end
    # Save the plots
    fnames = Vector{String}(undef, length(idx))
    for (i, j) in enumerate(idx)
        # i, j = 1, idx[1]
        if use_labels
            fnames[i] = labeltofname(label = plots.labels[j], prefix = prefix, suffix = format)
        else
            fnames[i] = labeltofname(label = string(i), prefix = prefix, suffix = format)
        end
        if isfile(fnames[i])
            throw(ErrorException("Cannot overwrite exisiting file: `" * fnames[i] * "`."))
        end
        CairoMakie.save(fnames[i], plots.plots[j])
    end
    fnames
end
