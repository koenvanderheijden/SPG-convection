module NoiseAuxiliaryFunctions

export calculate_ratio, integration, plot_integration_PDF

include("BornModel.jl")
import .BornModel as BM

using Plots, Plots.PlotMeasures
using DifferentialEquations
using Parameters
using Statistics
using KernelDensity

function calculate_ratio(M_yearly, threshold)
    # Split the array into two based on the condition (above or below threshold)
    M_above_threshold = [x for x in M_yearly if x > threshold]
    # M_below_threshold = [x for x in M_yearly if x <= threshold]

    # Count the lengths of the arrays
    n_above_threshold = length(M_above_threshold)
    n_total = length(M_yearly)

    # Calculate the ratio
    ratio = n_above_threshold / n_total
    return ratio
end

function count_below_above_threshold(M_yearly, threshold = 22)
    below_counts = Int[]
    above_counts = Int[]
    current_count = 1
    current_state = M_yearly[1] < threshold ? :below : :above

    for i in 2:length(M_yearly)
        if current_state == :below
            if M_yearly[i] < threshold
                current_count += 1
            else
                push!(below_counts, current_count)
                current_count = 1  # Reset count for new segment
                current_state = :above
            end
        else
            if M_yearly[i] >= threshold
                current_count += 1
            else
                push!(above_counts, current_count)
                current_count = 1  # Reset count for new segment
                current_state = :below
            end
        end
    end

    # Push the last segment count
    if current_state == :below
        push!(below_counts, current_count)
    else
        push!(above_counts, current_count)
    end

    return below_counts, above_counts
end


function integration(pars, noisevec, noise_seed, n_years, end_of_spinup = 10, threshold = 22, x0 = [BM.nondimensional_T(6), BM.nondimensional_S(35.5), BM.nondimensional_T(5), BM.nondimensional_S(35.5), 0, 0, 0])
    # x0 = [BM.nondimensional_T(6), BM.nondimensional_S(35.5), BM.nondimensional_T(5), BM.nondimensional_S(35.5), 0, 0, 0]
    timesteps_per_year = 3650
    dt = (1/timesteps_per_year)
    tspan = (0., n_years)

    sdeprob = SDEProblem(BM.stochastic_nondimensional_born!, BM.noisefunction(noisevec), x0, tspan, pars)
    sol = solve(sdeprob, EM(), dt = dt, seed = noise_seed)

    t = sol.t .- end_of_spinup

    offset = 2 #length(sol) % timesteps_per_year
    output_sol, output_t = sol[end_of_spinup*timesteps_per_year + offset:end], t[end_of_spinup*timesteps_per_year + offset:end]
    output_M = BM.dimensional_M(BM.calculate_M(output_sol, pars))

    M_yearly = zeros(Int(tspan[2]) - end_of_spinup)

    if length(output_M) > (n_years - end_of_spinup) * timesteps_per_year - 10
        for i in 1:length(M_yearly)
            M_yearly[i] = mean(output_M[(i-1) * timesteps_per_year + 1:(i*timesteps_per_year)])
        end
    end

    ratio = calculate_ratio(M_yearly, threshold)
    below, above = count_below_above_threshold(M_yearly, threshold)

    output = (M_yearly = M_yearly, ratio = ratio, below = below, above = above)
    # output = (sol = output_sol, t = output_t, M = output_M, M_yearly = M_yearly, ratio = ratio)
    return output
end

function plot_integration_PDF(output, title0 = "", title1 = "", title2 = "", title3 = "", xlimits12 = nothing, xlimits3 = nothing, ylimits12 = nothing, ylimits3 = nothing)
    # @unpack sol, t, M, M_yearly = output
    @unpack M_yearly, below, above = output

    p0 = plot(legend=false, framestyle=:none, grid=false, yticks = []) 
    annotate!(0, 0.5, text(title0, :black, 30, :left))

    xaxis = 4500:5000
    yfill = zeros(length(xaxis))
    fill!(yfill, 22)

    p1 = plot(xaxis, M_yearly[length(M_yearly) - 500 : end], label = "", grid = false, color = "grey65", 
        framestyle = :axes) # sol.t[timeslice:end] .- 10,
    plot!(xaxis, zeros(length(xaxis)), fillrange = yfill, color = "grey65", alpha = 0.15, label = "")
    xlims!(xaxis[1], xaxis[end])

    if ylimits12 !== nothing
        ylims!(ylimits12[1], ylimits12[2])
    end
    xlabel!("time [yr]")
    ylabel!("M [Sv]")
    title!(title1)

    k_M = kde(M_yearly)
    p2 = plot(k_M.density, k_M.x, lw = 2, grid = false, label = "", color = "dodgerblue", xticks = true)
    if xlimits12 !== nothing
        xlims!(xlimits12[1], xlimits12[2])
    end
    if ylimits12 !== nothing
        ylims!(ylimits12[1], ylimits12[2])
    end
    ylabel!("M [Sv]")
    xlabel!("density [a.u.]")
    title!(title2)
    
    p3 = plot(grid = false, yticks = true)
    if below != Int64[]
        k_below = kde(below) 
        plot!(k_below.x, k_below.density, lw = 2, color = "dodgerblue", label = "")
    end
    # if above != Int64[]
    #     k_above = kde(above)
    #     plot!(k_above.x, k_above.density, lw = 2, color = "dodgerblue", label = "")
    # end
    if xlimits3 !== nothing
        xlims!(xlimits3[1], xlimits3[2])
    end
    if ylimits3 !== nothing
        ylims!(ylimits3[1], ylimits3[2])
    end
    xlabel!("residence time [yr]")
    ylabel!("density [a.u.]")
    title!(title3)

    layout = @layout [a{0.075w} b{0.525w} c{0.2w} d{0.2w}]

    p = plot(p0, p1, p2, p3, dpi = 300, layout = layout,
        bottom_margin = 25px, left_margin = 20px, top_margin = 10px)

    return p
end

function plot_integration_PDF_noresidence(output, title0, title1 = "", title2 = "", xlimits12 = nothing, xlimits3 = nothing, ylimits12 = nothing, ylimits3 = nothing)
    # @unpack sol, t, M, M_yearly = output
    @unpack M_yearly, below, above = output

    p0 = plot(legend=false, framestyle=:none, grid=false, yticks = []) 
    annotate!(0, 0.5, text(title0, :black, 30, :left))

    xaxis = 4500:5000
    yfill = zeros(length(xaxis))
    fill!(yfill, 22)

    p1 = plot(xaxis, M_yearly[length(M_yearly) - 500 : end], label = "", grid = false, color = "grey65") # sol.t[timeslice:end] .- 10,
    if ylimits12 !== nothing
        ylims!(ylimits12[1], ylimits12[2])
    end
    plot!(xaxis, zeros(length(xaxis)), fillrange = yfill, color = "grey65", alpha = 0.15, label = "")
    xlims!(xaxis[1], xaxis[end])

    xlabel!("time [yr]")
    ylabel!("M [Sv]")
    title!(title1)

    k_M = kde(M_yearly)
    p2 = plot(k_M.density, k_M.x, lw = 2, grid = false, label = "", color = "dodgerblue", xticks = true)
    if xlimits12 !== nothing
        xlims!(xlimits12[1], xlimits12[2])
    end
    if ylimits12 !== nothing
        ylims!(ylimits12[1], ylimits12[2])
    end
    ylabel!("M [Sv]")
    xlabel!("density [a.u.]")
    title!(title2)
    
    layout = @layout [a{0.05w} b{0.75w} c{0.2w}]

    p = plot(p0, p1, p2, dpi = 300, layout = layout,
        bottom_margin = 25px, left_margin = 20px, top_margin = 10px)

    return p
end

end