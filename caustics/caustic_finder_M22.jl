#!/usr/ds/julia/bin/julia -t 1

using Plots, LaTeXStrings
pyplot()

include("ODE_solver_caustics_M22.jl")
println("ODE_solver_caustics_M22.jl loaded")
include("save_data.jl")
println("save_data.jl loaded")

# ARGS = [1D/2D, n, domain, time_max, time_interval, potential_center, potential_sd, A, D2, threshold, 
# freq_plot_intervals, group_name, run_name, mode]
# The argument "mode" can be "n" (no noise) or "s" (simulation). "n" is for creating a standard no-noise file that helps find focus point. 

# 1D Gaussian: 1, 2D Gaussian: 2
const gaussian_dim = ARGS[1]
const domain = tryparse(Float64,ARGS[2])
const time_max = tryparse(Float64, ARGS[3])
const time_interval = tryparse(Float64, ARGS[4])
const potential_center = tryparse(Float64, ARGS[5])
const potential_sd_x = tryparse(Float64, ARGS[6])
const potential_sd_y = tryparse(Float64, ARGS[7])
const A = tryparse(Float64,ARGS[8])
const D2 = tryparse(Float64,ARGS[9])
const threshold = tryparse(Float64,ARGS[10])
const group_name = ARGS[11]
const run_name = ARGS[12]

println("Parsed the arguments:")
println("group_name: $(group_name)")
println("run_name: $(run_name)")
if gaussian_dim == "1"
    println("Gaussian lens without y-grad")
    println("potential center = $(potential_center), potential sd = $(potential_sd_x)")
elseif gaussian_dim == "2"
    println("Gaussian lens with y-grad, potential center = ($(potential_center),0),
        potential sd = ($(potential_sd_x),$(potential_sd_y))")
end
println("Domain of simulation: [-$(domain), $(domain)] x [-$(domain), $(domain)]")
println("Simulation time: $(time_max)")
println("Time interval: $(time_interval)")
println("Simulation parameters: A=$A, D2=$(D2), threshold=$(threshold)")

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

function simulate()
    ant_coor_x, ant_M22 = ode_solver(speed)
    the_plot = plot(ant_coor_x, ant_M22, dpi=600, legends=false, xlabel=L"x", ylabel=L"M_{22}")
    savefig(the_plot, joinpath(group_directory, "$(run_name).png"))
end

const speed = 1.
const time_min = 0.0
tspan = (time_min, time_max)
simulate()
