#!/usr/ds/bin/julia

# This code is for the implementation for ants in correlated random potential

include("random_potential_methods.jl")
println("random_potential_methods.jl loaded")
include("ODE_solver_reduced.jl")
println("ODE_solver_reduced.jl loaded")
include("save_data.jl")
println("save_data.jl loaded")

using Plots
pyplot()

# We will need run_name_1 and run_name_2 
# Same group_name means it is to be averaged for the calculation of scintillation index; same run_name_1 is for same field

const n = tryparse(Int64,ARGS[1])
const time_max = tryparse(Float64, ARGS[2])
const time_interval = tryparse(Float64, ARGS[3])
const A = tryparse(Float64,ARGS[4])
const D2 = tryparse(Float64,ARGS[5])
const threshold = tryparse(Float64,ARGS[6])
const freq_plot_intervals = tryparse(Int64,ARGS[7])
const group_name = ARGS[8]
const run_name_1 = ARGS[9]
const run_name_2 = ARGS[10]

data_directory = "./data"
group_directory = joinpath(data_directory, "data_files", group_name)
run_1_directory = joinpath(group_directory, run_name_1)
plot_directory = joinpath(data_directory, "plots", group_name)

potential_jld2 = joinpath(run_1_directory, "random_potential.jld2")

mkpath(data_directory)
mkpath(group_directory)
mkpath(plot_directory)
mkpath(run_1_directory)

potential_args_file = joinpath(run_1_directory, "potential_args.json")
potential_args = read_data_from_json(potential_args_file)

const domain = tryparse(Float64, potential_args[1])
const rho = tryparse(Float64, potential_args[2])
const fourier_interval = tryparse(Int64, potential_args[3])
const seed = tryparse(Int64, potential_args[4])

println("Ants in random potential")
println("group_name: $(group_name)")
println("Number of ants, n = $(n)")
println("Domain of simulation: [-$(domain), $(domain)] x [-$(domain), $(domain)]")
println("Correlation length of random potential : $(rho)")
println("Number of intervals for Fourier transformation in generation of correlated random potential: $(fourier_interval)")
println("Random seed for correlated potential generation: $(seed)")
println("Simulation time: $(time_max)")
println("Time interval: $(time_interval)")
println("Simulation parameters: A=$A, D2=$(D2), threshold=$(threshold)")
println("Frequency plot dimensions: $(freq_plot_intervals) x $(freq_plot_intervals)")

const sim = "fixed_speed"

const x_min= -1*domain
const x_max = domain
const y_min = -domain
const y_max =domain
const interval_x = fourier_interval
const interval_y = fourier_interval
x, y = range(-domain, length=fourier_interval, domain), range(-domain, length=fourier_interval, domain)
x_hist, y_hist = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=freq_plot_intervals, domain)

#const r = 0
const β = 1
const pot_type = "pheromone"

function simulate(n)
    time_steps, ant_coor_x, ant_coor_y, _ = ode_solver(n, speed, init_pos = pos, init_direction = direction, random_direction = false, eq=pot_type, domain=domain, time_interval=time_interval, tspan=tspan)

    hist = get_freq_array(time_steps, ant_coor_x, ant_coor_y)
    save_data_to_jld2(time_steps, joinpath(run_1_directory, "$(run_name_2)_time_steps.jld2"))
    save_data_to_jld2(ant_coor_x, joinpath(run_1_directory, "$(run_name_2)_ant_coor_x.jld2"))
    save_data_to_jld2(ant_coor_y, joinpath(run_1_directory, "$(run_name_2)_ant_coor_y.jld2"))
    save_data_to_jld2(hist, joinpath(run_1_directory, "$(run_name_2).jld2"))
    save_data_to_json(ARGS, joinpath(run_1_directory, "$(run_name_2)_parameters.json"))
end

function plot_trajectories()
    arr_t_file = joinpath(run_1_directory, "$(run_name_2)_time_steps.jld2")
    arr_x_file = joinpath(run_1_directory, "$(run_name_2)_ant_coor_x.jld2")
    arr_y_file = joinpath(run_1_directory, "$(run_name_2)_ant_coor_y.jld2")
    arr_t = read_data_from_jld2(arr_t_file)
    arr_x = read_data_from_jld2(arr_x_file)
    arr_y = read_data_from_jld2(arr_y_file)
    potential_arr = [v(i, j) for j in y, i in x]
    the_plot = heatmap(x, y, potential_arr, dpi=300, legends=false)
    # the_plot = plot(dpi=300, legends=false)
    for i in 1:n
        total_time = arr_t[i]
        plot!(arr_x[i,1:total_time], arr_y[i,1:total_time], color=:white, linewidth=0.5)
    end
    savefig(the_plot, joinpath(plot_directory, "trajectory_$(run_name_1)_$(run_name_2).png"))
end

function plot_hist()
    hist_file = joinpath(run_1_directory, "$(run_name_2).jld2")
    hist = read_data_from_jld2(hist_file)
    hist_max = maximum(hist)
    the_plot = heatmap(x_hist, y_hist, hist, dpi=300, legends=false, clim=(0,hist_max/5))
    savefig(the_plot, joinpath(plot_directory, "hist_$(run_name_1)_$(run_name_2).png"))
end

potential = read_data_from_jld2(potential_jld2)
v = interpolate_potential(potential)

parallel = false

# potential = [v(i, j) for j in y, i in x]
# potential_plot = contour(x, y, potential, fill=true, dpi=300)

# parallel
if parallel
    pos = [ones(n).*(-domain), range(-domain, length=n, domain)]
    direction = 0.0
else
    pos = (0.0, 0.0)
    direction = rand(n) * 2 * pi
    # direction = range(0,length=n,2*pi*(1-1/n))
end

const speed = 1.
const time_min = 0.0
tspan = (time_min, time_max)
const alpha_max = nothing

simulate(n)
# plot_trajectories()
plot_hist()