using Plots; pyplot()
include("save_data.jl")

freq_plot_intervals = 500
domain = 40
x_hist, y_hist = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=freq_plot_intervals, domain)

plot_directory = "."

data_files_directory = "./data/data_files"
D2_arr = ["0.0","0.004","0.008","0.012","0.016","0.02","0.024"]
folders = [joinpath(data_files_directory, "D2=$(D2)", "seed=0") for D2 in D2_arr]

for (i, folder) in enumerate(folders)
    println(folder)
    files = [joinpath(folder, "run$(run).jld2") for run in 0:99]
    arr_sum = zeros(freq_plot_intervals,freq_plot_intervals)
    for file in files
        arr = read_data_from_jld2(file)
        arr_sum = arr_sum .+ arr
    end
    the_plot = heatmap(x_hist, y_hist, log10.(arr_sum), dpi=300)
    savefig(the_plot, joinpath(plot_directory, "single_hist_$(D2_arr[i]).png"))
end
