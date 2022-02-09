This folder contains the code that simulates correlated random potential.
Steps of using the codes:
1. Run run_random_potential_methods.jl with suitable arguments to create files of correlated random potential
Alternatively, run queuing_random_potential.py for a cluster implementation
2. After the folders are created and each of them contains the file random_potential.jld2 and potential_args.json, we can run cluster_random.jl (or plot trajectories/hist by changing the last three lines)
Alternatively, run queuing_simulate_random_potential for a cluster implementation
3. Run combine_jld2.jl (two arguments, group_name for simulation and new_group_name for saving jld2)
Can plot hist if wanted, but may have to change the source code (also in scintillation.jl) for correct domain and freq_plot_intervals
Alternatively, run queuing_combine_jld2.py for a cluster implementation
4. To stack scintillation indices to a single plot, run stack_scintillation.jl (may need to change source code)

*** seed is not random seed anymore due to change in code, so for each time the field is generated without seeds