module Plotting

export plot_branch, plot_special_points

include("BornModel.jl")
import .BornModel as BM
using Plots
using BifurcationKit

function convert_to_dimensional(values, quantity)
    """
    Converts an array/value of quantity in nondimensional form back to dimensional form. 
    Used in plot_special_points and plot_branch
    """
    ops = Dict(
        :S2 => BM.dimensional_S,
        :μ4 => BM.dimensional_F,
        :T2 => BM.dimensional_T,
        :M  => BM.dimensional_M,
        :Us => BM.dimensional_U,
        :Ud => BM.dimensional_U,
        :Δσ31 => BM.dimensional_σ,
        :Δσ13 => BM.dimensional_σ,
        :Δσ21 => BM.dimensional_σ,
        :Δσ12 => BM.dimensional_σ,
        :Δσ43 => BM.dimensional_σ,
        :Δσ34 => BM.dimensional_σ
    )

    conversion_function = get(ops, quantity, "unknown quantity")

    return conversion_function(values)
end

function get_marker(specialpoint_type)
    ops = Dict(
        :bp => :circle,
        :hopf => :dtriangle
    )

    return get(ops, specialpoint_type, :auto)
end

function plot_special_points(branch, param, var; make_dimensional = true, plot_legend = true)
    labels = []
    for i in 1:length(branch.specialpoint) 
        if branch.specialpoint[i].type != :endpoint # | branch.specialpoint[i].type != :hopf # was eerst &&?
            lab = string(branch.specialpoint[i].type)

            if lab == "bp"
                lab = "saddle node bifurcation"
            end

            if lab == "hopf"
                lab = "Hopf bifurcation"
            end

            if lab == "hh"
                lab = "Hopf-Hopf"
            end

            if lab in labels
                lab = ""
            end 

            if make_dimensional == true
                x_value_nondimensional = branch.specialpoint[i].param
                x_value = convert_to_dimensional(x_value_nondimensional, param)

                y_value_nondimensional = getfield(branch.specialpoint[i].printsol, var)
                y_value = convert_to_dimensional(y_value_nondimensional, var)
            else
                x_value = branch.specialpoint[i].param
                y_value = getfield(branch.specialpoint[i].printsol, var)
            end

            if plot_legend == true
                scatter!([x_value], [y_value], color = BifurcationKit.get_color(branch.specialpoint[i].type), m = get_marker(branch.specialpoint[i].type), label = lab, merge = true)
            else
                scatter!([x_value], [y_value], color = BifurcationKit.get_color(branch.specialpoint[i].type), m = get_marker(branch.specialpoint[i].type), label = "")
            end                

            push!(labels, lab)
        end
    end
end 

function plot_branch(branch, param::Symbol, var::Symbol; make_dimensional = true, plot_legend = true, label = "", xlabel = "", ylabel = "", title = "", xlimits = nothing, ylimits = nothing, dpi = 200, colstable = :dodgerblue, colunstable = :grey65, lwstable = 1, lwunstable = 1)
    col = [stb ? colstable : colunstable for stb in branch.stable]
    lw = [stb ? lwstable : lwunstable for stb in branch.stable]

    # converting back to dimensional form
    if make_dimensional == true
        x_values_nondimensional = branch.param
        x_values = convert_to_dimensional(x_values_nondimensional, param)

        y_values_nondimensional = getproperty(branch, var)
        y_values = convert_to_dimensional(y_values_nondimensional, var)
    else
        x_values = branch.param
        y_values = getproperty(branch, var)
    end

    p = plot(x_values, y_values, color = col, linewidth = lw, label = label, dpi = dpi, grid = false)
    title!(title)
    xlabel!(xlabel)
    ylabel!(ylabel)
     
    # adjusting x- and y-limits
    if xlimits !== nothing
        xlims!(xlimits[1], xlimits[2]) 
    end

    if ylimits !== nothing
        ylims!(ylimits[1], ylimits[2])
    end

    # plot the special points
    plot_special_points(branch, param, var; make_dimensional, plot_legend)

    return p
end

function plot_branch!(branch, param::Symbol, var::Symbol; make_dimensional = true, plot_legend = true, label = "", ls = :solid, xlabel = "", ylabel = "", title = "", xlimits = nothing, ylimits = nothing, dpi = 200, colstable = :dodgerblue, colunstable = :grey65, lwstable = 1, lwunstable = 1)
    col = [stb ? colstable : colunstable for stb in branch.stable]
    lw = [stb ? lwstable : lwunstable for stb in branch.stable]

    # converting back to dimensional form
    if make_dimensional == true
        x_values_nondimensional = branch.param
        x_values = convert_to_dimensional(x_values_nondimensional, param)

        y_values_nondimensional = getproperty(branch, var)
        y_values = convert_to_dimensional(y_values_nondimensional, var)
    else
        x_values = branch.param
        y_values = getproperty(branch, var)
    end

    p = plot!(x_values, y_values, color = col, linewidth = lw, linestyle = ls, label = label, dpi = dpi, grid = false)
    title!(title)
    xlabel!(xlabel)
    ylabel!(ylabel)
     
    # adjusting x- and y-limits
    if xlimits !== nothing
        xlims!(xlimits[1], xlimits[2]) 
    end

    if ylimits !== nothing
        ylims!(ylimits[1], ylimits[2])
    end

    # plot the special points
    plot_special_points(branch, param, var; make_dimensional, plot_legend)

    return p
end

end


####### backup 
# function convert_to_dimensional(values, quantity)
#     if quantity == :S2
#         converted_values = BM.dimensional_S(values)
#     elseif quantity == :μ4
#         converted_values = BM.dimensional_F(values)
#     elseif quantity == :T2
#         converted_values = BM.dimensional_T(values)
#     elseif quantity == :M
#         converted_values = BM.dimensional_M(values)
#     elseif quantity == :Δσ13 | quantity == :Δσ21 | quantity == :Δσ43 
#         converted_values = BM.dimensional_σ(values)
#     end

#     return converted_values
# end


# function plot_special_points(branch, var)
#     labels = []
#     for i in 1:length(branch.specialpoint) 
#         if branch.specialpoint[i].type != :endpoint # && branch.specialpoint[i].type != :hopf
#             lab = string(branch.specialpoint[i].type)
#             if lab in labels
#                 lab = ""
#             end 
#             value = getfield(branch.specialpoint[i].printsol, var)
#             scatter!([branch.specialpoint[i].param], [value], color = BifurcationKit.get_color(branch.specialpoint[i].type), label = lab, merge = true)
#             push!(labels, lab)
#         end
#     end
# end 

# function plot_branch(branch; var = :M, xlabel = "", ylabel = "", title = "", xlimits, ylimits, dpi = 200, colstable = :dodgerblue, colunstable = :grey65, lwstable = 1, lwunstable = 1)
#     col = [stb ? colstable : colunstable for stb in branch.stable]
#     lw = [stb ? lwstable : lwunstable for stb in branch.stable]

#     array = getproperty(branch, var)
#     p = plot(branch.param, array, color = col, linewidth = lw, label = "", title = title, xlabel = xlabel, ylabel = ylabel, dpi = dpi, grid = false)
#     xlims!(xlimits[1], xlimits[2])
#     ylims!(ylimits[1], ylimits[2])

#     labels = []
#     for i in 1:length(branch.specialpoint) 
#         if branch.specialpoint[i].type != :endpoint
#             lab = string(branch.specialpoint[i].type)
#             if lab in labels
#                 lab = ""
#             end 
#             value = getfield(branch.specialpoint[i].printsol, var)
#             scatter!([branch.specialpoint[i].param], [value], color = BifurcationKit.get_color(branch.specialpoint[i].type), label = lab, merge = true)
#             push!(labels, lab)
#         end
#     end
#     return p
# end