#!/usr/lmp/anaconda3/bin/python

# ARGS: domain, rho, fourier_interval, seed, group_name, run_name_1

import os
from math import ceil

domain = 40
rho = 1
fourier_interval = 512

for D2 in [0.000, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024]:
    group_name = "D2={}".format(D2)

    run_number = 20

    QSUB="qsub -b y -cwd -e jobs -o jobs -q rostam.q"

    command1 = " ".join(["time ./run_random_potential_methods.jl"])

    seed_list = [seed for seed in range(run_number)]
    run_name_1_list = ["seed={}".format(seed) for seed in seed_list]

    for i, seed in enumerate(seed_list):
        run_name_1 = run_name_1_list[i]
        whole_command = " ".join([QSUB, command1, str(domain), str(rho), str(fourier_interval), str(seed), group_name, run_name_1])
        print(whole_command)
        os.system(whole_command)