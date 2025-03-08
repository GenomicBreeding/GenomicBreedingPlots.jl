"""
    abstract type PlotsGB end

Abstract type representing the base type for various plotting structures in GBPlots.
"""
abstract type PlotsGB end

"""
    DistributionPlots <: PlotsGB

Mutable struct for storing distribution plots.

# Fields
- `labels::Vector{String}`: Vector containing labels for the plots
- `plots::Vector{CairoMakie.Figure}`: Vector containing the distribution plot figures
"""
mutable struct DistributionPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

"""
    ViolinPlots <: PlotsGB

Mutable struct for storing violin plots.

# Fields
- `labels::Vector{String}`: Vector containing labels for the plots
- `plots::Vector{CairoMakie.Figure}`: Vector containing the violin plot figures
"""
mutable struct ViolinPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

"""
    CorHeatPlots <: PlotsGB

Mutable struct for storing correlation heatmap plots.

# Fields
- `labels::Vector{String}`: Vector containing labels for the plots
- `plots::Vector{CairoMakie.Figure}`: Vector containing the correlation heatmap figures
"""
mutable struct CorHeatPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

"""
    TreePlots <: PlotsGB

Mutable struct for storing tree visualization plots.

# Fields
- `labels::Vector{String}`: Vector containing labels for the plots
- `plots::Vector{CairoMakie.Figure}`: Vector containing the tree plot figures
"""
mutable struct TreePlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

"""
    BarPlots <: PlotsGB

Mutable struct for storing bar plots.

# Fields
- `labels::Vector{String}`: Vector containing labels for the plots
- `plots::Vector{CairoMakie.Figure}`: Vector containing the bar plot figures
"""
mutable struct BarPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

"""
    BoxPlots <: PlotsGB

Mutable struct for storing box plots.

# Fields
- `labels::Vector{String}`: Vector containing labels for the plots
- `plots::Vector{CairoMakie.Figure}`: Vector containing the box plot figures
"""
mutable struct BoxPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

"""
    PCBiPlots <: PlotsGB

Mutable struct for storing principal component or biplot visualizations.

# Fields
- `labels::Vector{String}`: Vector containing labels for the plots
- `plots::Vector{CairoMakie.Figure}`: Vector containing the PC/biplot figures
"""
mutable struct PCBiPlots <: PlotsGB
    labels::Vector{String}
    plots::Vector{CairoMakie.Figure}
end

"""
    checkdims(x::PlotsGB)::Bool

Check if dimensions of labels and plots match in a PlotsGB object.

# Arguments
- `x::PlotsGB`: A PlotsGB object containing labels and plots

# Returns
- `Bool`: `true` if the number of labels equals the number of plots, `false` otherwise

"""
function GBCore.checkdims(x::PlotsGB)::Bool
    length(x.labels) == length(x.plots)
end

"""
    labeltofname(; label::String, prefix::String, suffix::String)::String

Convert a label string into a valid filename by replacing special characters.

This function takes a label string and converts it into a filename-safe string by:
1. Replacing common symbols with underscores
2. Adding a prefix and suffix
3. Cleaning up repeated separators

# Arguments
- `label::String`: The input label to be converted
- `prefix::String`: String to prepend to the filename
- `suffix::String`: The file extension (without the dot)

# Returns
- `String`: A cleaned filename string with format "prefix-label.suffix"
"""
function labeltofname(; label::String, prefix::String, suffix::String)::String
    symbol_strings::Vector{String} = [" ", "\n", "\t", "(", ")", "&", "|", ":", "=", "+", "-", "*", "/", "%"]
    for s in symbol_strings
        label = replace(label, s => "_")
    end
    fname = join([prefix, label], "-") * "." * suffix
    repeat_strings = ["..", "__", "--", "_.", "._", "-.", ".-"]
    for s in repeat_strings
        fname = replace(fname, s => ".")
    end
    fname
end

"""
    saveplots(plots::PlotsGB; idx::Vector{Int64}=[0], format::String="svg", 
              prefix::String="", use_labels::Bool=true, overwrite::Bool=false)::Vector{String}

Save plots from a PlotsGB object to files in the specified format.

# Arguments
- `plots::PlotsGB`: A PlotsGB object containing the plots to be saved
- `idx::Vector{Int64}`: Indices of plots to save. Default `[0]` saves all plots
- `format::String`: Output file format, one of "svg", "png", or "pdf". Default "svg"
- `prefix::String`: Prefix for output filenames. Default uses type name of plots
- `use_labels::Bool`: If true, use plot labels in filenames; if false, use numeric indices. Default true
- `overwrite::Bool`: If true, overwrite existing files; if false, throw error. Default false

# Returns
- `Vector{String}`: Vector of filenames where plots were saved

# Throws
- `ArgumentError`: If plots object is corrupted, indices are invalid, or format is unsupported
- `ErrorException`: If attempting to overwrite existing files when overwrite=false
"""
function saveplots(
    plots::PlotsGB;
    idx::Vector{Int64} = [0],
    format::String = "svg",
    prefix::String = "",
    use_labels::Bool = true,
    overwrite::Bool = false,
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
            if overwrite
                rm(fnames[i])
            else
                throw(ErrorException("Cannot overwrite exisiting file: `" * fnames[i] * "`."))
            end
        end
        try
            CairoMakie.save(fnames[i], plots.plots[j])
        catch
            continue
        end
    end
    fnames
end
