# using DifferentialEquations, PyCall, DelimitedFiles, ProgressMeter, Distributions, Dates
using DifferentialEquations, DelimitedFiles, Distributions
using ForwardDiff
# r=0

gradv(x,y) = ForwardDiff.gradient(z -> v(z[1], z[2]), [x, y])
hessianv(x,y) = ForwardDiff.hessian(z -> v(z[1],z[2]), [x, y])
caustic_points = []

function deterministic(du, u, p, t, A=A, D2=D2, threshold=threshold)
    
    # p = initial speed of ant    
    # pot = gradient of potential field (in 2-elements array [rho_x, rho_y])
    # vpot = potential field itself (as a floating point number)
    # (Notations may look confusing but it is the way it is)

    # u = [x,y,phi,M22,dotM22]

    pot = gradv(u[1], u[2])
    vpot = v(u[1], u[2])

    alpha_temp = A*(-pot[1]*sin(u[3])+pot[2]*cos(u[3]))/(v(u[1],u[2])+threshold)

    # When not using ForwardDiff:
    # alpha_temp = A*(-vx(u[1],u[2])*sin(u[3])+vy(u[1],u[2])*cos(u[3]))/(v(u[1],u[2])+threshold)

    if alpha_max == nothing
        alpha = alpha_temp
    else
        if abs(alpha_temp) >= alpha_max
            alpha = sign(alpha_temp)*alpha_max
        else
            alpha = alpha_temp
        end
    end

    # Adding three extra variables: M22, M32 and dotM22, where M22 and M32 are entries of monodromy matrix and dotM22 its time derivative
    # We have to define K31,K32,K33 before that, where dM/dt = KM, and to do that, we will need the second-order partial derivatives of potential

    rho_x, rho_y = pot[1], pot[2]
    hessian_matrix = hessianv(u[1],u[2])
    rho_xx = hessian_matrix[1,1]
    rho_xy = hessian_matrix[1,2]
    rho_yy = hessian_matrix[2,2]

    # Define matrix K
    K11, K12, K21, K22 = 0, 0, 0, 0
    K13 = -sin(u[3])
    K23 = cos(u[3])
    K31 = A*(-rho_xx*sin(u[3])+rho_xy*cos(u[3]))/(threshold+vpot) - A*rho_x*(-rho_x*sin(u[3]+rho_y*cos(u[3])))/(threshold+vpot)^2
    K32 = A*(rho_yy*cos(u[3])-rho_xy*sin(u[3]))/(threshold+vpot) - A*rho_y*(-rho_x*sin(u[3])+rho_y*cos(u[3]))/(threshold+vpot)^2
    K33 = A*(-rho_y*sin(u[3])-rho_x*cos(u[3]))/(threshold+vpot)

    K_matrix = [K11 K12 K13 ; K21 K22 K23 ; K31 K32 K33]
    M_matrix = [u[4] u[5] u[6] ; u[7] u[8] u[9] ; u[10] u[11] u[12]]
    KM_matrix = K_matrix*M_matrix

    # u[4] = M11, u[5] = M12, u[6] = M13, ..., u[11] = M32, u[12] = M33
    
    du[1] = p*cos(u[3])
    du[2] = p*sin(u[3])
    du[3] = alpha
    du[4] = KM_matrix[1,1]
    du[5] = KM_matrix[1,2]
    du[6] = KM_matrix[1,3]
    du[7] = KM_matrix[2,1]
    du[8] = KM_matrix[2,2]
    du[9] = KM_matrix[2,3]
    du[10] = KM_matrix[3,1]
    du[11] = KM_matrix[3,2]
    du[12] = KM_matrix[3,3]
end

function stochastic(du, u, p, t, A=A, D2=D2, threshold=threshold)
    
    du[1] = 0
    du[2] = 0
    du[3] = D2
    du[4],du[5],du[6],du[7],du[8],du[9],du[10],du[11],du[12]=0,0,0,0,0,0,0,0,0

end

# Condition to terminate the ODE if either the ant or the pheromone detection is out of the defined domain to prevent error
function hit_wall_condition(u, t, integrator, domain=domain)
    abs(u[1]) >= domain || abs(u[2]) >= domain
end


# General caustic condition, continuous callback
function caustic_condition(u, t, integrator)
    cos(u[3])*u[8] - sin(u[3])*u[5]
end

# Stop the integrator when an ant is out of domain
affect!(integrator) = terminate!(integrator)


# Initialize the speed randomly (Maxwell-Boltzmann distribution)

# function maxwell_speeds_initialize(mps)
#     stats = pyimport("scipy.stats")
#     maxwell = stats.maxwell
#     velocity = maxwell.rvs(size=1000, scale=speed_max*2^(-0.5))
#     histogram(velocity, leg=false)
# end

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
    
    hit_wall_cb = DiscreteCallback(hit_wall_condition,affect!)  # Terminate the ODE if the ant is out of domain
    caustic_cb = ContinuousCallback(caustic_condition,affect!)  # Terminate the ODE when hitting caustic condition
    # cb = CallbackSet(hit_wall_cb, caustic_cb)
    cb = caustic_cb

    # Simply define new ranges and cut out the boundary at which the ant will "detect" pheromone "out of defined domain"
    x_range = range(x_min, length=fourier_interval, stop=x_max)
    y_range = range(y_min, length=fourier_interval, stop=y_max)


    m = trunc(Int64,ceil(time_max/time_interval*1.2)+1)
    ant_coor_x = zeros(n,m)
    ant_coor_y = zeros(n,m)
    ant_direction = zeros(n,m)
    ant_cross_prod = zeros(n,m)
    
    time_steps = [0 for i in 1:n]

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
    else
        mps = speed
        v0 = maxwell.rvs(size=n, scale=mps*2^(-0.5))

    end

    # Each theta is actually an ant, or may be considered as some dimensionless angle (this is for starting at origin, not parallel case)
    for theta in 1:n
        if sim=="fixed_speed"
            p = v0
        else
            p = v0[theta]
        end

        # Initial condtions

        u0 = [init_pos[1][theta];init_pos[2][theta];init_direction[theta];1;0;0;0;1;0;0;0;1]

        if eq == "pheromone"
            prob_sde = SDEProblem(deterministic,stochastic,u0,tspan,p,callback=cb)
        elseif eq == "newtonian"
            prob_sde = SDEProblem(newtonian,stochastic,u0,tspan,p,callback=cb)
        end
        
        sol = solve(prob_sde, saveat=time_interval)

        for (bin, time) in enumerate(sol.t)
            ant_coor_x[theta, bin] = sol.u[bin][1]
            ant_coor_y[theta, bin] = sol.u[bin][2]
            ant_direction[theta, bin] = sol.u[bin][3]
            ant_cross_prod[theta, bin] = cos(sol.u[bin][3])*sol.u[bin][8] - sin(sol.u[bin][3])*sol.u[bin][5]
        end

        time_steps[theta] = length(sol.t)
    end
    
    if save_plot
        if sim=="fixed_speed"
            savefig(p_ants, join([directory, "Ants_plot_n=", n, ".png"]))
        elseif sim=="maxwell"
            savefig(p_ants, join([directory, "ants_", n, "_", mps, "_", D2, ".png"]))
        end
    end
    
    return (time_steps, ant_coor_x, ant_coor_y, ant_direction, ant_cross_prod)
end

function get_freq_array(time_steps, ant_coor_x, ant_coor_y)
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

    for (theta, time_step) in enumerate(time_steps)
        coor_x_arr = ant_coor_x[theta, 1:time_step]
        coor_y_arr = ant_coor_y[theta, 1:time_step]
        for (coor_x, coor_y) in zip(coor_x_arr, coor_y_arr)
            histpos = coor_to_histpos(coor_x, coor_y)
            hist[histpos[1],histpos[2]] = hist[histpos[1],histpos[2]] + 1
        end
    end
    return hist
end

function produce_freq_plot(time_steps, ant_coor_x, ant_coor_y)
    the_array = get_freq_array(time_steps, ant_coor_x, ant_coor_y)
    x_hist, y_hist = range(-domain, length=freq_plot_intervals, domain), range(-domain, length=no_freq_intervals, domain)
    traj_plot = heatmap(x_hist, y_hist, the_array', dpi=300, title=join(["Frequency plot (t_interval = ", time_interval, ", x,y_interval = ", domain/freq_plot_intervals]))
    savefig(traj_plot, join([directory, time_interval, "_", traj_plot_intervals, ".png"]))
end

function get_last_coor(time_steps, ant_coor_x, ant_coor_y, ant_direction, ant_cross_prod)
    last_coor_arr = zeros(4,n)
    for theta in 1:n
        last_time_step = time_steps[theta]
        last_coor_arr[1,theta] = ant_coor_x[theta, last_time_step]
        last_coor_arr[2,theta] = ant_coor_y[theta, last_time_step]
        last_coor_arr[3,theta] = ant_direction[theta, last_time_step]
        last_coor_arr[4,theta] = ant_cross_prod[theta, last_time_step]
    end
    return last_coor_arr
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
