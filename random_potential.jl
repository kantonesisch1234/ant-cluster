using DelimitedFiles, Interpolations, Plots, FileIO, FileIO
gr()

function generate_potential()
    command = `python generate_random_potential.py $domain $rho $fourier_interval $seed $sim` 
    run(command)
end

function load_potential()
    
    potential_filename = join([directory, "random_potential.txt"])
    
    potential = readdlm(potential_filename, Float64)
    potential = potential .- minimum(potential)
    
    return potential
    
end

function interpolate_potential(potential)

    x = range(x_min, length=fourier_interval, stop=x_max)
    y = range(y_min, length=fourier_interval, stop=y_max)

    v_unscaled = interpolate(potential, inter["cubic"])
    v_scaled = Interpolations.scale(v_unscaled, x, y)
    v = extrapolate(v_scaled, Periodic())
    
    return v
end

inter = Dict("constant" => BSpline(Constant()), 
"linear" => BSpline(Linear()), 
"quadratic" => BSpline(Quadratic(Line(OnCell()))),
"cubic" => BSpline(Cubic(Line(OnCell())))
)
    
