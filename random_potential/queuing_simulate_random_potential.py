#!/usr/lmp/anaconda3/bin/python

# ARGS: n, time_max, time_interval, A, D2, threshold, freq_plot_intervals, group_name, run_name_1, run_name_2

import os
from math import ceil

n = 10000
time_max = 100
time_interval = 0.01
A = 1
# D2 = 0.004
threshold = 0.05
freq_plot_intervals = 500

for D2 in [0.008, 0.012, 0.016, 0.020, 0.024]:

    group_name = "D2={}".format(D2)

    run_number_per_seed = 100

    seed_list = range(20)
    run_name_1_list = ["seed={}".format(seed) for seed in seed_list]
    run_name_2_list = ["run{}".format(i) for i in range(run_number_per_seed)]

    QSUB="qsub -b y -cwd -e jobs -o jobs -q rostam.q"

    command1 = " ".join(["time ./cluster_random.jl"])

    for run_name_1 in run_name_1_list:
        for run_name_2 in run_name_2_list:
            whole_command = " ".join([QSUB, command1, str(n), str(time_max), str(time_interval), str(A), str(D2), str(threshold), 
            str(freq_plot_intervals), group_name, run_name_1, run_name_2])
            print(whole_command)
            # os.system(whole_command)