# ant-cluster

Run simulation:
Arguments for `cluster.jl`: [1D/2D, n, domain, time_max, time_interval, potential_center, potential_sd, A, D2, threshold, freq_plot_intervals, group_name, run_name, mode]

1D/2D: potential with/without y-gradient (1 or 2)

group_name: folder saving the files

run_name: name of data files

mode: n (no noise, with an extra `focus_idx.jld2` file indicating position of focus point, should choose group_name of `standard`); s (normal simulation)

Example: run the command `julia cluster.jl 2 10000 40 100 0.01 -30 0.5 1 0 0.05 50 standard run1 n`

Example: run the command `julia cluster.jl 2 10000 40 100 0.01 -30 0.5 1 0.002 0.05 50 D2=0.002 run1 s`. 

The Julia files needed are `ODE_solver_reduced.jl`, `cluster.jl` and `save_data.jl`. 


Run plotting:

One argument for `cluster_analysis.jl`: group_name

Example: run the command `julia cluster_analysis.jl D2=0.002`. Only works if the folder `standard` exists with a file `focus_idx.jld2` inside it.

