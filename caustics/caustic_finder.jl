#!/usr/ds/julia/bin/julia -t 1

include("ODE_solver_caustics.jl")
println("ODE_solver_caustics.jl loaded")
include("save_data.jl")
println("save_data.jl loaded")

# ARGS = [1D/2D, n, domain, time_max, time_interval, potential_center, potential_sd, A, D2, threshold, 
# freq_plot_intervals, group_name, run_name, mode]
# The argument "mode" can be "n" (no noise) or "s" (simulation). "n" is for creating a standard no-noise file that helps find focus point. 

# 1D Gaussian: 1, 2D Gaussian: 2
const gaussian_dim = ARGS[1]
const n = tryparse(Int64,ARGS[2])
const domain = tryparse(Float64,ARGS[3])
const init_y_min = tryparse(Float64,ARGS[4])
const init_y_max = tryparse(Float64,ARGS[5])         # y_min and y_max specifies the range of initial y-coor of trajectories 
const time_max = tryparse(Float64, ARGS[6])
const time_interval = tryparse(Float64, ARGS[7])
const potential_center = tryparse(Float64, ARGS[8])
const potential_sd_x = tryparse(Float64, ARGS[9])
const potential_sd_y = tryparse(Float64, ARGS[10])
const A = tryparse(Float64,ARGS[11])
const D2 = tryparse(Float64,ARGS[12])
const threshold = tryparse(Float64,ARGS[13])
const freq_plot_intervals = tryparse(Int64,ARGS[14])
const group_name = ARGS[15]
const run_name = ARGS[16]
const mode = ARGS[17]

@assert mode == "n" || mode == "s"

println("Parsed the arguments:")
if mode == "n"
    println("Mode \"n\": This run is for locating the focus point.")
elseif mode == "s"
    println("Mode \"s\": This run is for actual simulation.")
end
println("group_name: $(group_name)")
println("run_name: $(run_name)")
if gaussian_dim == "1"
    println("Gaussian lens without y-grad")
    println("potential center = $(potential_center), potential sd = $(potential_sd_x)")
elseif gaussian_dim == "2"
    println("Gaussian lens with y-grad, potential center = ($(potential_center),0),
        potential sd = ($(potential_sd_x),$(potential_sd_y))")
end
println("Number of ants, n = $(n)")
println("Domain of simulation: [-$(domain), $(domain)] x [-$(domain), $(domain)]")
println("Range of initial y-coors: [$(init_y_min), $(init_y_max)]")
println("Simulation time: $(time_max)")
println("Time interval: $(time_interval)")
println("Simulation parameters: A=$A, D2=$(D2), threshold=$(threshold)")
println("Frequency plot dimensions: $(freq_plot_intervals) x $(freq_plot_intervals)")

# flush(stdout)

data_directory = "./data"
json_file_directory = joinpath(data_directory, "data_files")
group_directory = joinpath(json_file_directory, group_name)
mkpath(data_directory)
mkpath(json_file_directory)
mkpath(group_directory)

rho, fourier_interval, seed, sim = 1, 512, 1, "fixed_speed"
const x_min= -1*domain
const x_max = domain
const y_min = -domain
const y_max =domain
const interval_x = fourier_interval
const interval_y = fourier_interval
x, y = range(-domain, length=fourier_interval, domain), range(-domain, length=fourier_interval, domain)

#const r = 0
const β = 1
pot_type = "pheromone"

# Constructing the Gaussian potential and drawing the potential contour plots here
const potential_amp = 1

function v(x,y)
    # Standard normal distribution (Corrected with 1/2 term!)
    function snd(center, sd)
        return exp(-(x-center[1])^2/(2*sd[1]^2)-(y-center[2])^2/(2*sd[2]^2))
    end

    function snd_no_y_grad(center, sd)
        return exp(-((x-center)^2/(2*sd^2)))
    end

    if gaussian_dim == "1"
        return potential_amp*snd_no_y_grad(potential_center, potential_sd_x)
    elseif gaussian_dim == "2"
        potential_center_tuple = (potential_center, 0)
        potential_sd_tuple = (potential_sd_x, potential_sd_y)
        return potential_amp*snd(potential_center_tuple, potential_sd_tuple)
    end
end

# Needed when not using ForwardDiff
# vx(x,y) = -potential_amp*(x-potential_center)*v(x,y)/potential_sd_x^2

# function vy(x,y)
#     if gaussian_dim == "1"
#         return 0
#     elseif gaussian_dim == "2"
#         return -potential_amp*y*v(x,y)/potential_sd_y^2
#     end
# end

function simulate(n)
    time_steps, ant_coor_x, ant_coor_y, ant_direction, ant_cross_prod = ode_solver(n, speed, init_pos = pos, init_direction = direction, random_direction = false, eq=pot_type)
    hist = get_freq_array(time_steps, ant_coor_x, ant_coor_y)
    last_coor_arr = get_last_coor(time_steps, ant_coor_x, ant_coor_y, ant_direction, ant_cross_prod)
    save_data_to_jld2(hist, joinpath(group_directory, "$(run_name).jld2"))
    save_data_to_jld2(last_coor_arr, joinpath(group_directory, "$(run_name)_lastcoor.jld2"))
    save_data_to_json(ARGS, joinpath(group_directory, "$(run_name)_parameters.json"))

    focus_idx = find_focus(hist)
    if mode == "n"
        save_data_to_jld2(focus_idx, joinpath(group_directory, "focus_idx.jld2"))
    end
end

# Function to find focus
function find_focus(hist)
    no_noise_density_along_x_axis = hist[:,trunc(Int, freq_plot_intervals/2)]
    focus_x = -domain + (findmax(no_noise_density_along_x_axis)[2]-1)*2*domain/freq_plot_intervals
    focus_coor = (focus_x,0)
    focus_idx = (trunc(Int,(domain+focus_coor[1])/(2*domain/freq_plot_intervals)), trunc(Int, freq_plot_intervals/2))
    return focus_idx
end

pos = [ones(n).*(-domain), range(init_y_min, length=n, init_y_max)]

const direction = 0.0
const speed = 1.
const time_min = 0.0
tspan = (time_min, time_max)
const alpha_max = nothing

simulate(n)
