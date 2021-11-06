using DifferentialEquations, ForwardDiff

gradv(x,y) = ForwardDiff.gradient(z -> v(z[1], z[2]), [x, y])
hessianv(x,y) = ForwardDiff.hessian(z -> v(z[1],z[2]), [x, y])

function deterministic(du, u, p, t, A=A, D2=D2, threshold=threshold)
    
    # p = initial speed of ant    
    # pot = gradient of potential field (in 2-elements array [rho_x, rho_y])
    # vpot = potential field itself (as a floating point number)
    # (Notations may look confusing but it is the way it is)

    x = u[1]
    y = 0

    pot = gradv(x,y)
    vpot = v(x,y)

    # Adding three extra variables: M22, M32 and dotM22, where M22 and M32 are entries of monodromy matrix and dotM22 its time derivative
    # We have to define K31,K32,K33 before that, where dM/dt = KM, and to do that, we will need the second-order partial derivatives of potential

    rho_x, rho_y = pot[1], pot[2]
    hessian_matrix = hessianv(x,y)
    rho_xx = hessian_matrix[1,1]
    rho_xy = hessian_matrix[1,2]
    rho_yy = hessian_matrix[2,2]

    # K31 = A*rho_xy/(threshold+vpot) - A*rho_x*rho_y/(threshold+vpot)^2
    K32 = A*rho_yy/(threshold+vpot) - A*rho_y^2/(threshold+vpot)^2
    K33 = -A*rho_x/(threshold+vpot)

    # u[2] = M22, u[3] = M32
    
    du[1] = p
    du[2] = u[3]
    du[3] = K32*u[2] + K33*u[3]
end

function stochastic(du, u, p, t, A=A, D2=D2, threshold=threshold)
    
    du[1] = 0
    du[2] = 0
    du[3] = 0

end

"""
Ants start moving at the origin in the same direction but different initial speeds. 
Total number of ants = n
"""

# init_pos: Tuple - gives the same initial position for all ants; Arrays - [[x1,x2, ..., xn],[y1,y2,...,yn]]
# init_direction: nothing - if random: random; if not random: uniform; Arrays - [theta1, theta2, ..., thetan]; Real number - all same direction
# eq = "newtonian" or "pheromone"

function ode_solver(speed; domain=domain, time_interval=time_interval, tspan=tspan)
    
    # speed means speed of an ant in fixed_speed case and mps in Maxwell distribution case
    # direction takes in an angle, only applies for n=1

    # Simply define new ranges and cut out the boundary at which the ant will "detect" pheromone "out of defined domain"
    x_range = range(x_min, length=fourier_interval, stop=x_max)
    y_range = range(y_min, length=fourier_interval, stop=y_max)

    m = trunc(Int64,ceil(time_max/time_interval*1.2)+1)
    ant_coor_x = []
    ant_M2 = []
    
    v0 = speed
 
    # Each theta is actually an ant, or may be considered as some dimensionless angle (this is for starting at origin, not parallel case)
    p = v0

    # Initial condtions
    u0 = [-domain;1;0]
    prob_sde = SDEProblem(deterministic,stochastic,u0,tspan,p)

    sol = solve(prob_sde, saveat=time_interval)

    for (bin, time) in enumerate(sol.t)
        push!(ant_coor_x, sol.u[bin][1])
        push!(ant_M2, sol.u[bin][2])
    end
    
    return (ant_coor_x, ant_M2)
end