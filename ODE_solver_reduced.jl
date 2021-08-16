# using DifferentialEquations, PyCall, DelimitedFiles, ProgressMeter, Distributions, Dates
using DifferentialEquations, DelimitedFiles, Distributions, Dates
r=0

function deterministic(du, u, p, t, A=A, D2=D2, threshold=threshold)
    
    # p = initial speed of ant    
    alpha_temp = A*(-vx(u[1],u[2])*sin(u[3])+vy(u[1],u[2])*cos(u[3]))/(v(u[1],u[2])+threshold)
    if alpha_max == nothing
        alpha = alpha_temp
    else
        if abs(alpha_temp) >= alpha_max
            alpha = sign(alpha_temp)*alpha_max
        else
            alpha = alpha_temp
        end
    end
    
    du[1] = p*cos(u[3])
    du[2] = p*sin(u[3])  
    du[3] = alpha
    
end

function stochastic(du, u, p, t, A=A, D2=D2, threshold=threshold)
    
    du[1] = 0
    du[2] = 0
    du[3] = D2

end

function newtonian(du, u, p, t)
    
    pot = gradv(u[1], u[2])
    
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -A*pot[1]
    du[4] = -A*pot[2]
    du[5] = atan(u[4]/u[3])
    
end

# Condition to terminate the ODE if either the ant or the pheromone detection is out of the defined domain to prevent error
function condition(u, t, integrator, domain=domain)
    abs(u[1]) >= domain || abs(u[2]) >= domain
end

affect!(integrator) = terminate!(integrator)

# Initialize the speed randomly (Maxwell-Boltzmann distribution)

function maxwell_speeds_initialize(mps)
    stats = pyimport("scipy.stats")
    maxwell = stats.maxwell
    velocity = maxwell.rvs(size=1000, scale=speed_max*2^(-0.5))
    histogram(velocity, leg=false)
end

"""
Ants start moving at the origin in the same direction but different initial speeds. 
Total number of ants = n
"""

# init_pos: Tuple - gives the same initial position for all ants; Arrays - [[x1,x2, ..., xn],[y1,y2,...,yn]]
# init_direction: nothing - if random: random; if not random: uniform; Arrays - [theta1, theta2, ..., thetan]; Real number - all same direction
# eq = "newtonian" or "pheromone"

function ode_solver(n, speed; random_direction=true, init_pos=(0,0), init_direction=nothing, sim="fixed_speed", save_plot=false, eq="pheromone", domain=domain, time_interval=time_interval, tspan=tspan)
    
    # speed means speed of an ant in fixed_speed case and mps in Maxwell distribution case
    # direction takes in an angle, only applies for n=1
    
    cb = DiscreteCallback(condition,affect!)  # Terminate the ODE if the ant is out of domain

    # Simply define new ranges and cut out the boundary at which the ant will "detect" pheromone "out of defined domain"
    x_range = range(x_min, length=fourier_interval, stop=x_max)
    y_range = range(y_min, length=fourier_interval, stop=y_max)

    ant_coor_x = [[] for i=1:n]
    ant_coor_y = [[] for i=1:n]
#     ant_v_x = [[] for i=1:n]
#     ant_v_y = [[] for i=1:n]
    ant_direction = [[] for i=1:n]

    if isa(init_pos,Tuple{Real, Real})
        init_pos = [ones(n).*init_pos[1], ones(n).*init_pos[2]]
    end
    
    if init_direction==nothing
        if random_direction
            init_direction = rand(Uniform(-π,π),n)
        else
            init_direction = collect(0:n-1).*2π/n
        end
    elseif isa(init_direction, Real)
        init_direction = ones(n).*init_direction
    end
        

   # Plot the contour plot of pheromone density field again, the ant trajectories will be plotted above the contour plot.


    if sim == "fixed_speed"
        v0 = speed
        # p_ants = contour(x_range, y_range, potential, fill=true,
        # title=join([n, " ants on pheromone field, D2 = ", D2]), dpi=300)
    else
        mps = speed
        v0 = maxwell.rvs(size=n, scale=mps*2^(-0.5))
        # p_ants = contour(x_range, y_range, potential, fill=true,
        # title=join([n, " ants on pheromone field , mps = ", mps, ", D2 = ", D2]), dpi=300)

    end

    # Each theta is actually an ant, or may be considered as some dimensionless angle 
    for theta in 1:n
        if sim=="fixed_speed"
            p = v0
        else
            p = v0[theta]
        end

#         u0 = [init_pos[1][theta];init_pos[2][theta];p*cos(init_direction[theta]);p*sin(init_direction[theta]);init_direction[theta]]
        u0 = [init_pos[1][theta];init_pos[2][theta];init_direction[theta]]
        if eq == "pheromone"
            prob_sde = SDEProblem(deterministic,stochastic,u0,tspan,p,callback=cb)
        elseif eq == "newtonian"
            prob_sde = SDEProblem(newtonian,stochastic,u0,tspan,p,callback=cb)
        end
        
        sol = solve(prob_sde, saveat=time_interval)


        """ 
        To save the coordinates of all time in saveat is a little bit tricky here. 
        Each time when a callback is called, the element in sol.t repeats itself so to make the visualization of time evolution right,
        the repeated time moments have to be gotten rid of. time_cache serves for this purpose. 
        """
        time_cache = -1

        for (bin, time) in pairs(sol.t)
            if time != time_cache
                push!(ant_coor_x[theta], sol.u[bin][1])
                push!(ant_coor_y[theta], sol.u[bin][2])
#                 push!(ant_v_x[theta], sol.u[bin][3])
#                 push!(ant_v_y[theta], sol.u[bin][4])
                push!(ant_direction[theta], sol.u[bin][3])
            time_cache = time

            end
        end
        
#         if save_plot
#             plot!(sol,vars=(1,2), legend=false, color=:white, linewidth = 1, linestyle = :dot, xlims=(-domain+r, domain-r), ylims=(-domain+r, domain-r))
#         end

    end
    
    if save_plot
        if sim=="fixed_speed"
            savefig(p_ants, join([directory, "Ants_plot_n=", n, ".png"]))
        elseif sim=="maxwell"
            savefig(p_ants, join([directory, "ants_", n, "_", mps, "_", D2, ".png"]))
        end
    end
    
    return (ant_coor_x, ant_coor_y, ant_direction)
end

function get_freq_array(ant_coor_x, ant_coor_y)
    hist = zeros(freq_plot_intervals, freq_plot_intervals)
    function coor_to_histpos(x_coor, y_coor)
        unit_length = domain*2/freq_plot_intervals
        histpos_x = trunc(Int, (x_coor+domain)/unit_length+1)
        histpos_y = trunc(Int, (y_coor+domain)/unit_length+1)
        if histpos_x > freq_plot_intervals
            histpos_x = freq_plot_intervals
        end
        if histpos_y > freq_plot_intervals
            histpos_y = freq_plot_intervals
        end
        if histpos_x <= 0
            histpos_x = 1
        end
        if histpos_y <= 0
            histpos_y = 1
        end
        return histpos_x, histpos_y
    end
    for (i1, theta) in enumerate(ant_coor_x)
        for (i2, time_pt) in enumerate(theta)
            histpos = coor_to_histpos(ant_coor_x[i1][i2], ant_coor_y[i1][i2])
            hist[histpos[1],histpos[2]] = hist[histpos[1],histpos[2]] + 1
        end
    end
    return hist
end

function produce_freq_plot(ant_coor_x, ant_coor_y)
    the_array = get_freq_array(ant_coor_x, ant_coor_y)
    x_hist, y_hist = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=no_freq_intervals, domain)
    traj_plot = heatmap(x_hist, y_hist, the_array', dpi=300, title=join(["Frequency plot (t_interval = ", time_interval, ", x,y_interval = ", domain/freq_plot_intervals]))
    savefig(traj_plot, join([directory, time_interval, "_", traj_plot_intervals, ".png"]))
end

# coor chooses whether to plot vy against y or v_theta against theta. Fill in either "y" or "theta".
function plot_phase_space(t;typ="plot",coor="y")
    t_index = trunc(Int, t/time_interval+1)
    if coor=="y"
        yt = []
        vyt = []
        for i in 1:n
            push!(yt, ant_coor_y[i][t_index])
            push!(vyt, ant_v_y[i][t_index])
        end
        if typ=="plot"
            plot(yt, vyt, legend=false, dpi=300, title=join(["Time: ", round(t; digits=1)]), xlabel="y", ylabel="vy")
        elseif typ=="scatter"
            scatter(yt, vyt, legend=false, dpi=300, title=join(["Time: ", round(t; digits=1)]), xlabel="y", ylabel="vy", markersize=1)
        end
    elseif coor=="theta"
#         cart_to_polar_pos(vx,vy) = atan(vy, vx)
        cart_to_polar_pos(x,y) = atan(y, x)
        function cart_to_polar_vel(x,y,vx,vy)
            return (x*vy-y*vx)*(x^2+y^2)^(-1)
#             pot = gradv(x,y)
#             return A/speed*(-pot[1]*vy+pot[2]*vx)/v(x,y)
        end
        theta_t = []
        v_theta_t = []
        
        for i in 1:n
            push!(theta_t, cart_to_polar_pos(ant_coor_x[i][t_index],ant_coor_y[i][t_index]))
#             push!(theta_t, cart_to_polar_pos(ant_v_x[i][t_index],ant_v_y[i][t_index]))
            push!(v_theta_t, cart_to_polar_vel(ant_coor_x[i][t_index],ant_coor_y[i][t_index],ant_v_x[i][t_index],ant_v_y[i][t_index]))
        end
        if typ=="plot"
            plot(theta_t, v_theta_t, legend=false, dpi=300, title=join(["Time: ", round(t; digits=1)]), xlabel="theta", ylabel="d(theta)/dt")
        elseif typ=="scatter"
            scatter(theta_t, v_theta_t, legend=false, dpi=300, title=join(["Time: ", round(t; digits=1)]), xlabel="theta", ylabel="d(theta)/dt", markersize=0.1)
        end
        
    end
        
end

"""
# To make animation, some tricks are used in the following. Not the main point here, I think.
function make_ants_gif(ant_coor_x, ant_coor_y, ant_direction;fps=10, bg=true, trail=true, arrows=true)
    anim = @animate for (bin, t) in pairs(time_min:time_interval:time_max)

        ant_x = [bin<length(ant_coor_x[theta]) ? ant_coor_x[theta][bin] : ant_coor_x[theta][end] for theta=1:n]  
        ant_y = [bin<length(ant_coor_y[theta]) ? ant_coor_y[theta][bin] : ant_coor_y[theta][end] for theta=1:n]
        ant_angle = [bin<length(ant_direction[theta]) ? ant_direction[theta][bin] : ant_direction[theta][end] for theta=1:n]
        
        if bg
            color=:white
            x_range, y_range = range(x_min, length=fourier_interval, stop=x_max), range(y_min, length=fourier_interval, stop=y_max)
            contour(x_range, y_range, v, fill=true, 
                title=join(["n = ", n, ", threshold = ", threshold, ", Time = ", round((bin-1)*time_interval; digits=1)]),dpi=300)
            scatter!(ant_x, ant_y, lims=(-domain+r, domain-r), color=color, markersize=1, leg=false)
            if arrows
                quiver!(ant_x, ant_y, quiver=(cos.(ant_angle), sin.(ant_angle)), color=color)
            end
        else
            color=:black
            scatter(ant_x, ant_y, lims=(-domain+r, domain-r), color=color, markersize=1, leg=false,
            title=join(["n = ", n, ", D2 = ", D2, ", Time = ", round((bin-1)*time_interval; digits=1)]),dpi=300)
            if arrows
                quiver!(ant_x, ant_y, quiver=(cos.(ant_angle), sin.(ant_angle)), color=color)
            end
        end
        
        if trail
            @showprogress 1 for theta in 1:n
                if bin<length(ant_coor_x[theta]) && bin<length(ant_coor_y[theta])
                    plot!(ant_coor_x[theta][1:bin], ant_coor_y[theta][1:bin], color=color, linestyle=:dot, linewidth=0.5)
                else
                    plot!(ant_coor_x[theta][1:end], ant_coor_y[theta][1:end], color=color, linestyle=:dot, linewidth=0.5)
                end
            end
        end
    end
    if trail
        gif(anim, join([directory, "anim_", n, "_", D2, "_", threshold, ".gif"]), fps = fps)
    else
        gif(anim, join([directory, "anim_", n, "_", D2, "_", threshold, "_no_trails", ".gif"]), fps = fps)
    end
end

function screenshot_freq(t_screenshot)
    function ant_coor_trim(ant_coor)
        
    end
    
    
end


function screenshot(t_screenshot; bg=true, trail=true, arrows=true, show_plot=false, xlims=(-domain+r,domain-r), ylims=(-domain+r,domain-r),extra_info="")
    bin_screenshot = round(Int64, (t_screenshot - time_min)/time_interval + 1)

    ant_x = [bin_screenshot<length(ant_coor_x[theta]) ? ant_coor_x[theta][bin_screenshot] : ant_coor_x[theta][end] for theta=1:n]           
    ant_y = [bin_screenshot<length(ant_coor_y[theta]) ? ant_coor_y[theta][bin_screenshot] : ant_coor_y[theta][end] for theta=1:n]
    if arrows
        ant_angle = [bin_screenshot<length(ant_direction[theta]) ? ant_direction[theta][bin_screenshot] : ant_direction[theta][end] for theta=1:n]
    end
    
    if bg
        color=:white
        x_range, y_range = range(x_min, length=fourier_interval, stop=x_max), range(y_min, length=fourier_interval, stop=y_max)
        screenshot=contour(x_range, y_range, v, fill=true, 
            title=join(["n = ", n, ", threshold = ", threshold, ", Time = ", round((bin_screenshot-1)*time_interval; digits=1), extra_info]), dpi=300)
        scatter!(ant_x, ant_y, xlims=xlims, ylims=ylims,color=color, markersize=1, leg=false)
        if arrows
            quiver!(ant_x, ant_y, quiver=(cos.(ant_angle), sin.(ant_angle)), color=color)
        end
    else
        color=:black
        screenshot = scatter(ant_x, ant_y, xlims=xlims, ylims=ylims, color=:black, markersize=1, leg=false,
            title=join(["n = ", n, ", threshold = ", threshold, ", Time = ", round((bin_screenshot-1)*time_interval; digits=1), extra_info]), dpi=300) 
        if arrows
            quiver!(ant_x, ant_y, quiver=(cos.(ant_angle), sin.(ant_angle)), color=color)
        end
    end
    if trail
        @showprogress 1 for theta in 1:n
            if bin_screenshot<length(ant_coor_x[theta]) && bin_screenshot<length(ant_coor_y[theta])
                plot!(ant_coor_x[theta][1:bin_screenshot], ant_coor_y[theta][1:bin_screenshot], color=color, linestyle=:dot, linewidth=0.5)
            else
                plot!(ant_coor_x[theta][1:end], ant_coor_y[theta][1:end], color=color, linestyle=:dot, linewidth=0.5)
            end
        end
    end
    
    if show_plot
        display(screenshot)
    end
    extra_info_str = join(split(extra_info, "/"))
    savefig(screenshot, join([directory, "screenshot_", t_screenshot, "_", n, "_", threshold, extra_info_str, ".png"]))
end
"""

function create_parameters_file(;dir=directory)
    open(join([dir, "/parameters.txt"]), "w") do f
	write(f, join(["Time of simulation: ", Dates.now(), "\n\n"]))
        write(f, join(["n = ", n, "\n"]))
        write(f, join(["A = ", A, "\n"]))
        write(f, join(["D2 = ", D2, "\n"]))
        write(f, join(["domain = ", domain, "\n"]))
        write(f, join(["rho = ", rho, "\n"]))
        write(f, join(["fourier_interval = ", fourier_interval, "\n"]))
        write(f, join(["seed = ", seed, "\n"]))
        write(f, join(["tspan = ", tspan, "\n"]))
        write(f, join(["A = ", A, "\n"]))
        write(f, join(["time_interval = ", time_interval, "\n"]))
        write(f, join(["freq_plot_intervals = ", freq_plot_intervals, "\n"]))
        write(f, join(["bg_potential = ", bg_potential, "\n"]))
        write(f, join(["pot_type = ", pot_type, "\n"]))
        try
            write(f, join(["potential_center = ", potential_center, "\n"]))
            write(f, join(["potential_sd = ", potential_sd, "\n"]))
            write(f, join(["potential_amp = ", potential_amp, "\n"]))
            write(f, join(["alpha_max = ", alpha_max, "\n"]))
            write(f, join(["threshold = ", threshold, "\n"]))
            write(f, join(["pos = ", pos, "\n"]))
            write(f, join(["direction = ", direction, "\n"]))
        catch
            
            
        end

    end
end
