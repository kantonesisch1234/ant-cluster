#!/usr/lmp/julia/bin/julia
include("ODE_solver_reduced_new.jl")
println("ODE_solver_reduced_new.jl loaded")
include("random_potential.jl")
println("random_potential.jl loaded")
include("save_data.jl")
println("save_data.jl loaded")

# ARGS = [1D/2D, n, domain, time_max, time_interval, potential_center, potential_sd, A, D2, threshold, freq_plot_intervals, group_name, run_name]

# 1D Gaussian: 1, 2D Gaussian: 2
gaussian_dim = ARGS[1]
n = tryparse(Int64,ARGS[2])
domain = tryparse(Float64,ARGS[3])
time_max = tryparse(Float64, ARGS[4])
time_interval = tryparse(Float64, ARGS[5])
potential_center = tryparse(Float64, ARGS[6])
potential_sd = tryparse(Float64, ARGS[7])
A = tryparse(Float64,ARGS[8])
D2 = tryparse(Float64,ARGS[9])
threshold = tryparse(Float64,ARGS[10])
freq_plot_intervals = tryparse(Int64,ARGS[11])
group_name = ARGS[12]
run_name = ARGS[13]
chunks = tryparse(Int64,ARGS[14])

println("Parsed the arguments:")
println("group_name: $(group_name)")
println("run_name: $(run_name)")
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
println("Divide the simulation into $(chunks) chunks.")

function create_dir(dir)
    if !isdir(dir)
        mkdir(dir)
    end
end

data_directory = "./data"
json_file_directory = joinpath(data_directory, "json_files")
group_directory = joinpath(json_file_directory, group_name)
create_dir(data_directory)
create_dir(json_file_directory)
create_dir(group_directory)

rho, fourier_interval, seed, sim = 1, 512, 1, "fixed_speed"
x_min, x_max, y_min, y_max, interval_x, interval_y = -1*domain, domain, -1*domain, domain, fourier_interval, fourier_interval
x, y = range(-domain, length=fourier_interval, domain), range(-domain, length=fourier_interval, domain)

r = 0
const β = 1
pot_type = "pheromone"

# Constructing the Gaussian potential and drawing the potential contour plots here
potential_amp = 1

function v(x,y)
    # Standard normal distribution (Corrected with 1/2 term!)
    function snd(center, sd)
        return exp(-(x-center[1])^2/(2*sd[1]^2)-(y-center[2])^2/(2*sd[2]^2))
    end
    
    function snd_no_y_grad(center, sd)
        return exp(-((x-center)^2/(2*sd^2)))
    end
    
    if gaussian_dim == "1"
        return potential_amp*snd_no_y_grad(potential_center, potential_sd)
    elseif gaussian_dim == "2"
        potential_center_tuple = (potential_center, 0)
        potential_sd_tuple = (potential_sd, domain)
        return potential_amp*snd(potential_center_tuple, potential_sd_tuple)
    end
end

# potential = [v(i, j) for j in y, i in x]
# potential_plot = contour(x, y, potential, fill=true, dpi=300)

pos = [ones(n).*(-domain), range(-domain/2, length=n, domain/2)]
direction = 0
speed = 1.
time_min = 0.0
tspan = (time_min, time_max)
alpha_max = nothing

hist = ode_solver(n, speed, chunks, init_pos = pos, init_direction = direction, random_direction = false, eq=pot_type)
save_data_to_json(ARGS, joinpath(group_directory, "$(run_name)_parameters.json"))
save_data_to_json(hist, joinpath(group_directory, "$(run_name).json"))
