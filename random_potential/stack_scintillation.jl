#!/usr/ds/julia/bin/julia

using Plots, LaTeXStrings
pyplot()
include("save_data.jl")

domain = 40

new_group_directory = joinpath(".", "data", "data_files", "random_potential_averaged")
plot_directory = joinpath(".", "data", "plots", "random_potential_averaged")

filenames = [joinpath(new_group_directory,"D2=$(D2)_scintillation.jld2") for D2 in ["0.0","0.004", "0.008","0.012","0.016","0.02","0.024"]]
labels = ["0.000","0.004","0.008","0.012","0.016","0.020","0.024"]

r_domain = 0.1:0.1:domain
scintillation_plot = plot(xlabel=L"r", ylabel=L"S(r)", legends=:outerright, dpi=300)

for (i,file) in enumerate(filenames)
    scintillation_arr = read_data_from_jld2(file)
    plot!(r_domain, scintillation_arr, label=L"D_2=%$(labels[i])")
end

savefig(scintillation_plot, "scintillation_index.png")


