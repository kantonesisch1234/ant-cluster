#!/usr/ds/bin/julia

# Run this file to create the correlated random potential fields

include("random_potential_methods.jl")
include("save_data.jl")

const domain = tryparse(Float64,ARGS[1])
const rho = tryparse(Float64,ARGS[2])
const fourier_interval = tryparse(Int64,ARGS[3])
const seed = tryparse(Int64,ARGS[4])
const group_name = ARGS[5]
const run_name_1 = ARGS[6]

data_directory = "./data"
group_directory = joinpath(data_directory, "data_files", group_name)
run_name_1_directory = joinpath(group_directory, run_name_1)

potential_arr = generate_potential_2()
save_data_to_jld2(potential_arr, joinpath(run_name_1_directory, "random_potential.jld2"))
# rm(joinpath(run_name_1_directory, "random_potential.txt"))
save_args()