#!/usr/lmp/anaconda3/bin/python

import numpy as np
import sys
import os
import json

domain, rho, fourier_interval, seed, group_name, run_name_1 = float(sys.argv[1]), float(sys.argv[2]), int(float(sys.argv[3])), int(float(sys.argv[4])), sys.argv[5], sys.argv[6]

V0 = 1




np.random.seed(seed)

x_min, x_max, y_min, y_max, interval_x, interval_y = -1*domain, domain, -1*domain, domain, fourier_interval, fourier_interval

x = np.linspace(x_min, x_max, interval_x)
y = np.linspace(y_min, y_max, interval_y)

X, Y = np.meshgrid(x, y)

# Autocorrelation function
def c(x, y, V0, rho):
    return V0**2 * np.exp(-(x**2+y**2)/rho**2)

# Only apply for even interval_x and interval_y, zero paddings and symmetry about the origin needed for real output
matrix1 = np.random.rand(X.shape[0]//2-1, X.shape[1]//2-1)
matrix2 = np.random.rand(X.shape[0]//2-1, X.shape[1]//2-1)
phi1 = np.c_[np.zeros((X.shape[0]//2-1, 1)), matrix1, np.zeros((X.shape[0]//2-1, 1)), matrix2]
matrix3 = -np.flipud(np.fliplr(matrix2))
matrix4 = -np.flipud(np.fliplr(matrix1))
phi2 = np.c_[np.zeros((X.shape[0]//2-1, 1)), matrix3, np.zeros((X.shape[0]//2-1, 1)), matrix4]
phi = np.r_[np.zeros((1, X.shape[1])), phi1, np.zeros((1, X.shape[1])), phi2] 

# Discrete autocorrelation function in terms of numpy arrays
C = c(X, Y, V0, rho)

phase = np.exp(2j*np.pi*phi)

# Generated V(x, y) in discrete numpy array form
V = np.fft.ifft2(np.power(np.fft.fft2(C), 0.5)*phase)

data_directory = "./data"
group_directory = os.path.join(data_directory, "data_files", group_name)
run_name_1_directory = os.path.join(group_directory, run_name_1)
plot_directory = os.path.join(data_directory, "plots", group_name)
potential_txt_path = os.path.join(run_name_1_directory, "random_potential.txt")

# if not os.path.exists(group_path):
#     os.makedirs(group_path)
if not os.path.exists(run_name_1_directory):
    os.makedirs(run_name_1_directory)
if not os.path.exists(plot_directory):
    os.makedirs(plot_directory)

# Save the generated random potential array into a txt file
np.savetxt(potential_txt_path, V.real)

# Save the parameters
# with open(args_json_path, 'w') as f:
#     json.dump(sys.argv[1:], f)