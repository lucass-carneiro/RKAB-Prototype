using HDF5
using LinearAlgebra
using SummationByPartsOperators

function sw_phi(A, kx, t, x)
    omega = sqrt(kx * kx)
    return A * cos(2 * pi * (kx * x - omega * t))
end

function sw_pi(A, kx, t, x)
    omega = sqrt(kx * kx)
    return 2 * A * omega * pi * sin(2 * pi * (kx * x - omega * t))
end

function sw_dx(A, kx, t, x)
    omega = sqrt(kx * kx)
    return -2 * A * kx * pi * sin(2 * pi * (kx * x - omega * t))
end

function compute_rhs!(D, y, dy)
    Pi = y[2]
    Dx = y[3]
    
    Phi_rhs = dy[1]
    Pi_rhs = dy[2]
    Dx_rhs = dy[3]
    
    copy!(Phi_rhs, Pi)
    mul!(Pi_rhs, D, Dx)
    mul!(Dx_rhs, D, Pi)
end

function euler_step!(dt, y, dy)
    Phi = y[1]
    Pi = y[2]
    Dx = y[3]
    
    Phi_rhs = dy[1]
    Pi_rhs = dy[2]
    Dx_rhs = dy[3]

    @. Phi = Phi + dt * Phi_rhs
    @. Pi = Pi + dt * Pi_rhs
    @. Dx = Dx + dt * Dx_rhs
end

function evolve()
    # IO
    h5_file = h5open("output.h5", "w")
    state_group = create_group(h5_file, "state")
    grid_group = create_group(h5_file, "grid")
    
    # Domain
    x0 = -1.0
    xf = 1.0
    num_pts = 51
    dx = (xf - x0) / (num_pts - 1)

    # Derivative operator
    D = derivative_operator(
        DienerDorbandSchnetterTiglio2007(),
        derivative_order=1,
        accuracy_order=2,
        xmin=x0,
        xmax=xf,
        N=num_pts
    )

    # Write Domain
    write(grid_group, "x_coords", collect(grid(D)))

    # Time steps
    cfl = 0.25
    dt = cfl * dx
    last_iter = 100

    # State
    Phi = zeros(Float64, num_pts)
    Pi = zeros(Float64, num_pts)
    Dx = zeros(Float64, num_pts)
    y = [Phi, Pi, Dx]

    # RHS
    Phi_rhs = zeros(Float64, num_pts)
    Pi_rhs = zeros(Float64, num_pts)
    Dx_rhs = zeros(Float64, num_pts)
    dy = [Phi_rhs, Pi_rhs, Dx_rhs]
    
    # Standing wave params
    A = 1.0
    kx = 0.25
    
    # Init State
    for i in 0:(num_pts - 1)
        x = x0 + i * dx
        Phi[i + 1] = sw_phi(A, kx, 0.0, x)
        Pi[i + 1] = sw_pi(A, kx, 0.0, x)
        Dx[i + 1] = sw_dx(A, kx, 0.0, x)
    end

    # Write Initial state
    write(state_group, "Phi_$(lpad(0, 4, "0"))", Phi)

    # Iterate
    for i in 1:last_iter
        t = i * dt
        println("Iteration $i, t = $t")

        # Compute RHS
        compute_rhs!(D, y, dy)

        # Time Step
        euler_step!(dt, y, dy)

        # Enforce BCs
        Phi[1] = sw_phi(A, kx, 0.0, x0)
        Phi[end] = sw_phi(A, kx, 0.0, xf)

        Pi[1] = sw_pi(A, kx, 0.0, x0)
        Pi[end] = sw_pi(A, kx, 0.0, xf)
        
        Dx[1] = sw_dx(A, kx, 0.0, x0)
        Dx[end] = sw_dx(A, kx, 0.0, xf)

        # Write state
        write(state_group, "Phi_$(lpad(i, 4, "0"))", Phi)
    end

    close(h5_file)
end

function plot()
    h5_file = h5open("output.h5", "r")
    x_coords = read(h5_file["grid/x_coords"])
    Phi = read(h5_file["state/Phi_0001"])
    println(Phi)
    println(x_coords)
    close(h5_file)
end

function main()
    evolve()
    #plot()
end

main()