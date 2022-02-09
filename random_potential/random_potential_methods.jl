#!/usr/ds/bin/julia

using DelimitedFiles, Interpolations, FileIO, FFTW
include("save_data.jl")
include("randpot.jl")

# cluster_python_path = "/usr/ds/julia/bin/julia"

function generate_potential()
    lc = rho/(2*domain)
    ϵ=1.0
    Nx,Ny = fourier_interval, fourier_interval
    corr=GaussCorr(n=[Nx,Ny],s=ϵ,l=lc)
    pot=RandPot2D{GaussCorr}(corr)
    pot_arr = pot.u
    pot_arr = pot_arr .- minimum(pot_arr)
    return pot_arr
end

function generate_potential_2()
    # This function generates the random potential based on convolution between uniform white noise and Gaussian function
    N = fourier_interval
    gauss(x,y,lc) = exp(-(x*x+y*y)/(2*lc^2))
    non_negative = false
    u = zeros(N,N)
    lc = rho/(2*domain)
    x_range, y_range = range(-0.5, length=N, 0.5), range(-0.5, length=N, 0.5)
    while !non_negative
        f = [gauss(x,y,lc) for x in x_range, y in y_range]
        p = rand(N,N)
        u = real.(ifft(fft(f).*fft(p)))/N
        non_negative = (minimum(u) >= 0)
    end
    return u
end

# function load_potential()
    
#     potential_filename = joinpath(run_name_1_directory, "random_potential.txt")
    
#     potential = readdlm(potential_filename, Float64)
#     potential = potential .- minimum(potential)
    
#     return potential
    
# end

function interpolate_potential(potential)

    x = range(x_min, length=fourier_interval, stop=x_max)
    y = range(y_min, length=fourier_interval, stop=y_max)

    v_unscaled = interpolate(potential, inter["cubic"])
    v_scaled = Interpolations.scale(v_unscaled, x, y)
    v = extrapolate(v_scaled, Periodic())
    
    return v
end

function save_args()
    save_data_to_json(ARGS, joinpath(run_name_1_directory, "potential_args.json"))
end

inter = Dict("constant" => BSpline(Constant()), 
"linear" => BSpline(Linear()), 
"quadratic" => BSpline(Quadratic(Line(OnCell()))),
"cubic" => BSpline(Cubic(Line(OnCell())))
)

