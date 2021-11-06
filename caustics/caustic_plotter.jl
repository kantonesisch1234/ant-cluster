#!/usr/ds/julia/bin/julia -t 1
include("save_data.jl")
println("save_data.jl loaded")
using Plots
pyplot()

# This code serves to plot only one run instead of an entire group, so we need to specify also the run_name
# This code serves to plot the structure of caustics

function create_dir(dir)
    if !isdir(dir)
        mkpath(dir)
    end
end

function plot_caustics_to_2dmap(freq_array, lastcoor_array, title, filename)
    filepath = joinpath(plot_path, join([filename, ".png"]))
    x_freq, y_freq = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=freq_plot_intervals, domain)
    clim = (0,maximum(freq_array))
    freq_plot = heatmap(x_freq, y_freq, freq_array', dpi=600, clim=clim, title=title, legends=false, xlims=(-10,10), ylims=(-10,10))
    lastcoor_x, lastcoor_y = lastcoor_array[1,:], lastcoor_array[2,:]
    scatter!(lastcoor_x, lastcoor_y, markersize=5, markercolor=:green)
    savefig(freq_plot, filepath)
end

function plot_caustics_2d(lastcoor_array, title, filename)
    filepath = joinpath(plot_path, join([filename, ".png"]))
    lastcoor_x, lastcoor_y = lastcoor_array[1,:], lastcoor_array[2,:]
    caustics_plot = scatter(lastcoor_x, lastcoor_y, markersize=5, legends=false)
    savefig(caustics_plot, filepath)
end

directory = "./data"
group_name = ARGS[1]
run_name = ARGS[2]
group_path = joinpath(directory, "data_files", group_name)
run_path = joinpath(group_path, "$(run_name).jld2")
lastcoor_path = joinpath(group_path, "$(run_name)_lastcoor.jld2")
plot_path = joinpath(directory, "plots", group_name)

# json_files = [joinpath(group_path, file) for file in readdir(group_path) if occursin("json",file)]

create_dir(plot_path)

# print out the simulation parameters; they should be the same for the same "group" so we only print them once by choosing the first parameter file
# read the file storing simulation parameters

parameter_filepath = joinpath(group_path, "$(run_name)_parameters.json")
simulation_args = read_data_from_json(parameter_filepath)

# store the parameters
const gaussian_dim = tryparse(Int64,simulation_args[1])
const n = tryparse(Int64, simulation_args[2])
const domain = tryparse(Float64,simulation_args[3])
const time_max = tryparse(Float64, simulation_args[4])
const init_y_min = tryparse(Float64, simulation_args[5])
const init_y_max = tryparse(Float64, simulation_args[6])
const time_interval = tryparse(Float64, simulation_args[7])
const potential_center = tryparse(Float64, simulation_args[8])
const potential_sd_x = tryparse(Float64, simulation_args[9])
const potential_sd_y = tryparse(Float64, simulation_args[10])
const A = tryparse(Float64,simulation_args[11])
const D2 = tryparse(Float64,simulation_args[12])
const threshold = tryparse(Float64,simulation_args[13])
const freq_plot_intervals = tryparse(Int64,simulation_args[14])

println("Parsed the arguments:")
println("group_name: $(group_name)")
if gaussian_dim == "1"
    println("Gaussian lens without y-grad")
    println("potential center = $(potential_center), potential sd = $(potential_sd)")
elseif gaussian_dim == "2"
    println("Gaussian lens with y-grad, potential center = ($(potential_center),0),
        potential sd = ($(potential_sd),0)")
end
println("Number of ants, n = $(n)")
println("Domain of simulation: [-$(domain), $(domain)] x [-$(domain), $(domain)]")
println("Initial y-coors range: [$(init_y_min), $(init_y_max)]")
println("Simulation time: $(time_max)")
println("Time interval: $(time_interval)")
println("Simulation parameters: A=$A, D2=$(D2), threshold=$(threshold)")
println("Frequency plot dimensions: $(freq_plot_intervals) x $(freq_plot_intervals)")
println("\n")

# hist_arr = []

# for jld2_file in jld2_files
#     # read the run_name
#     run_name = get_filename_wo_extension(jld2_file)
#     println("Reading data from $(run_name).jld2 ...")
#     hist = read_data_from_jld2(jld2_file)
#     push!(hist_arr, hist)
# end
# total_hist = sum(hist_arr)

println("Reading data from $(run_name).jld2 ...")
hist = read_data_from_jld2(run_path)
lastcoor = read_data_from_jld2(lastcoor_path)

# plot_caustics_to_2dmap(hist, lastcoor, "D2=$(D2)", "2D_density_D2=$(D2)")
plot_caustics_2d(lastcoor, "D2=$(D2)", "caustics_plot")