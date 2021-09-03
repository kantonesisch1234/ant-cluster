include("ODE_solver_reduced.jl")
println("ODE_solver_reduced.jl loaded")
include("save_data.jl")
println("save_data.jl loaded")
using Plots

function get_parameter_filepath(jld2_file)
    parameter_filename = "$(splitext(basename(jld2_file))[1])_parameters.json"
    return joinpath(group_path, parameter_filename)
end

function create_dir(dir)
    if !isdir(dir)
        # mkdir(dir)
        mkpath(dir)
    end
end

function plot_2D_density(freq_array, title, filename)
    filepath = joinpath(plot_path, join([filename, ".png"]))
    x_freq, y_freq = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=freq_plot_intervals, domain)
    clim = (0,maximum(freq_array))
    freq_plot = heatmap(x_freq, y_freq, freq_array', dpi=600, clim=clim, title=title)
    savefig(freq_plot, filepath)
end

function plot_1D_density(freq_array, dim, title, filename, focus_idx; xlims=nothing, ylims=nothing)
    filepath = joinpath(plot_path, join([filename, ".png"]))
    x_freq, y_freq = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=freq_plot_intervals, domain)
    
    if dim=="y"
        density_1D = freq_array[focus_idx[1],:]
        denplot = plot(x_freq,density_1D, dpi=600, title=title, legend=false,
        xlabel = "y", ylabel = "rho(y)", xlims=xlims, ylims=ylims)
        savefig(denplot, filepath)
    elseif dim=="x"
        density_1D = freq_array[:,focus_idx[2]]
        denplot = plot(x_freq,density_1D, dpi=600 , xlabel = "x", ylabel = "rho(x)",
            title=title, legend=false, xlims=xlims, ylims=ylims)
        savefig(denplot, filepath)
    end
end

# we have to locate the focus point using no-noise trajectories, so we will have to create the group_name "no_noise" as the standard


# cluster.jl: ARGS = [1D/2D, n, domain, time_max, time_interval, potential_center, potential_sd, A, D2, threshold, freq_plot_intervals, group_name, run_name]
# cluster_analysis.jl: ARGS = group_name

directory = raw".\data"
group_name = ARGS[1]
group_path = joinpath(directory, "data_files", group_name)
plot_path = joinpath(directory, "plots", group_name)

# filepath for focus_idx
focus_idx_file = joinpath(directory, "data_files", "standard", "focus_idx.jld2")
const focus_idx = read_data_from_jld2(focus_idx_file)

# println([joinpath(group_path, file) for file in readdir(group_path)])
jld2_files = [joinpath(group_path, file) for file in readdir(group_path) if occursin("jld2",file)]

create_dir(plot_path)

# print out the simulation parameters; they should be the same for the same "group" so we only print them once by choosing the first parameter file
# read the file storing simulation parameters
first_jld2_file = jld2_files[1]
parameter_filepath = get_parameter_filepath(first_jld2_file)
simulation_args = read_data_from_json(parameter_filepath)

# store the parameters
const gaussian_dim = tryparse(Int64,simulation_args[1])
const n = tryparse(Int64, simulation_args[2])
const domain = tryparse(Float64,simulation_args[3])
const time_max = tryparse(Float64, simulation_args[4])
const time_interval = tryparse(Float64, simulation_args[5])
const potential_center = tryparse(Float64, simulation_args[6])
const potential_sd = tryparse(Float64, simulation_args[7])
const A = tryparse(Float64,simulation_args[8])
const D2 = tryparse(Float64,simulation_args[9])
const threshold = tryparse(Float64,simulation_args[10])
const freq_plot_intervals = tryparse(Int64,simulation_args[11])
# const group_name = simulation_args[12]

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
println("Simulation time: $(time_max)")
println("Time interval: $(time_interval)")
println("Simulation parameters: A=$A, D2=$(D2), threshold=$(threshold)")
println("Frequency plot dimensions: $(freq_plot_intervals) x $(freq_plot_intervals)")
println("\n")

hist_arr = []

for jld2_file in jld2_files
    # read the run_name
    run_name = get_filename_wo_extension(jld2_file)
    println("Reading data from $(run_name).jld2 ...")
    hist = read_data_from_jld2(jld2_file)
    push!(hist_arr, hist)
end
total_hist = sum(hist_arr)
average_hist = total_hist/length(jld2_files)
plot_1D_density(average_hist, "x", "D2=$(D2)", "1D_density_x_D2=$(D2)", focus_idx)
plot_1D_density(average_hist, "y", "D2=$(D2)", "1D_density_y_D2=$(D2)", focus_idx)
plot_2D_density(average_hist, "D2=$(D2)", "2D_density_D2=$(D2)")