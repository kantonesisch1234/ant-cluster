#!/usr/ds/julia/bin/julia

using Plots, LaTeXStrings
pyplot()
include("save_data.jl")
include("scintillation.jl")

function get_run_jld2_files(group_directory)
    items = [item for item in walkdir(group_directory)]
    jld2_files = []
    for item in items
        cwd = item[1]
        files = [joinpath(cwd, file) for file in item[3] if file != potential_file_name && splitext(file)[end] == ".jld2"]
        jld2_files = vcat(jld2_files, files)
    end
    return jld2_files
end

function save_averaged_jld2(jld2_files)
    file_numbers = size(jld2_files)[1]
    println("For verification, file_numbers = $(file_numbers)")
    total_hist = zeros(size(read_data_from_jld2(jld2_files[1]))...)
    for file in jld2_files
        arr = read_data_from_jld2(file)
        total_hist = total_hist + arr
    end
    averaged_arr = total_hist / file_numbers
    save_data_to_jld2(averaged_arr, joinpath(new_group_directory, "$(group_name)_avg.jld2"))
end

function save_averaged_squared_jld2(jld2_files)
    file_numbers = size(jld2_files)[1]
    println("For verification, file_numbers = $(file_numbers)")
    total_hist = zeros(size(read_data_from_jld2(jld2_files[1]))...)
    for file in jld2_files
        arr = read_data_from_jld2(file)
        total_hist = total_hist + arr.^2
    end
    averaged_arr = total_hist / file_numbers
    save_data_to_jld2(averaged_arr, joinpath(new_group_directory, "$(group_name)_avg_sq.jld2"))
end


function plot_hist()
    hist_file_avg = joinpath(new_group_directory, "$(group_name)_avg.jld2")
    hist_avg = read_data_from_jld2(hist_file_avg)
    # hist_file_avg_sq = joinpath(new_group_directory, "$(group_name)_avg_sq.jld2")
    # hist_max = maximum(hist)
    # the_plot = heatmap(x_hist, y_hist, hist, dpi=300, legends=false, clim=(0,hist_max * max_clim_ratio))
    the_plot = heatmap(x_hist, y_hist, log10.(hist_avg), dpi=300)
    savefig(the_plot, joinpath(data_directory, "plots","averaged_hist_clim_log_$(group_name)_.png"))
end

function save_scintillation()
    hist_file_avg = joinpath(new_group_directory, "$(group_name)_avg.jld2")
    hist_file_avg_sq = joinpath(new_group_directory, "$(group_name)_avg_sq.jld2")
    hist_avg = read_data_from_jld2(hist_file_avg)
    hist_avg_sq = read_data_from_jld2(hist_file_avg_sq)
    avg_density = hist_to_density(hist_avg)
    avg_sq_density = hist_to_density(hist_avg_sq)
    r_domain = 0.1:0.1:domain
    scintillation_arr = [scintillation_index(avg_density, avg_sq_density, r) for r in r_domain]
    save_data_to_jld2(scintillation_arr, joinpath(new_group_directory, "$(group_name)_scintillation.jld2"))
end

function plot_scintillation()
    r_domain = 0.1:0.1:domain
    scintillation_arr = read_data_from_jld2(joinpath(new_group_directory, "$(group_name)_scintillation.jld2"))
    scintillation_plot = plot(r_domain, scintillation_arr, xlabel=L"r", ylabel=L"S(r)", legends=false, dpi=300)
    savefig(scintillation_plot, joinpath(plot_directory, "$(group_name)_scintillation_index.png"))
end


group_name = ARGS[1]
new_group_name = ARGS[2]   # the group name for saving the averaged jld2
potential_file_name = "random_potential.jld2"  # We have to make sure that we don't wrongly use the potential for averaging
data_directory = "./data"
group_directory = joinpath(data_directory, "data_files", group_name)
plot_directory = joinpath(data_directory, "plots", group_name)
new_group_directory = joinpath(data_directory, "data_files", new_group_name)
mkpath(new_group_directory)


# For plotting hist and S(r), change if needed
domain = 40
freq_plot_intervals = 500

x_hist, y_hist = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=freq_plot_intervals, domain)

jld2_files = get_run_jld2_files(group_directory)
println("JLD2 files:")
for file in jld2_files
    println(file)
end
# save_averaged_jld2(jld2_files)
# save_averaged_squared_jld2(jld2_files)
plot_hist()
save_scintillation()
# plot_scintillation()