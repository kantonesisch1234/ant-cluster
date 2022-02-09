#!/usr/ds/julia/bin/julia

using QuadGK, Interpolations
include("save_data.jl")

# Change whenever needed
domain = 40
freq_plot_intervals = 500

x_min, x_max = -domain, domain
y_min, y_max = -domain, domain
x = range(x_min, length=freq_plot_intervals, stop=x_max)
y = range(y_min, length=freq_plot_intervals, stop=y_max)

function hist_to_density(hist)
    density_unscaled = interpolate(hist, inter["constant"])
    density_scaled = Interpolations.scale(density_unscaled, x, y)

    return density_scaled
end

function density_mean(density, r)
    integral, _ = quadgk(theta -> density(r*cos(theta), r*sin(theta)), 0, 2*pi, rtol=1e-8)
    integral = integral/(2*pi)
    return integral
end

function scintillation_index(avg_density, avg_sq_density, r)
    return density_mean(avg_sq_density,r) * density_mean(avg_density,r)^(-2) - 1
end

inter = Dict("constant" => BSpline(Constant()), 
"linear" => BSpline(Linear()), 
"quadratic" => BSpline(Quadratic(Line(OnCell()))),
"cubic" => BSpline(Cubic(Line(OnCell())))
)