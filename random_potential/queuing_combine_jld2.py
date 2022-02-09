#!/usr/lmp/anaconda3/bin/python

# ARGS: group_name, new_group_name
import os
from math import ceil

new_group_name = "random_potential_averaged"

for D2 in [0., 0.004, 0.008, 0.012, 0.016, 0.020, 0.024]:
    group_name = "D2={}".format(D2)

    QSUB="qsub -b y -cwd -e jobs -o jobs -q rostam.q"

    command1 = " ".join(["time ./combine_jld2.jl"])

    whole_command = " ".join([QSUB, command1, group_name, new_group_name])
    print(whole_command)
    os.system(whole_command)