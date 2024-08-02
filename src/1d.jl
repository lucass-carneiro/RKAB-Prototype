using HDF5
using LinearAlgebra
using SummationByPartsOperators

function sw_Phi(A, kx, t, x)
    omega = sqrt(kx * kx)
    return A*cos(2*omega*pi*t)*sin(2*kx*pi*x)
end

function sw_Pi(A, kx, t, x)
    omega = sqrt(kx * kx)
    return -2*A*omega*pi*sin(2*omega*pi*t)*sin(2*kx*pi*x)
end

function sw_Dx(A, kx, t, x)
    omega = sqrt(kx * kx)
    return 2*A*kx*pi*cos(2*omega*pi*t)*cos(2*kx*pi*x)
end

struct GridFuncs
    Phi::Array{Float64}
    Pi::Array{Float64}
    Dx::Array{Float64}

    function GridFuncs(num_points)
        new(
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points)
        )
    end
end

struct Substeps
    k0_Phi::Array{Float64}
    k0_Pi::Array{Float64}
    k0_Dx::Array{Float64}
    
    k1_Phi::Array{Float64}
    k1_Pi::Array{Float64}
    k1_Dx::Array{Float64}
    
    k2_Phi::Array{Float64}
    k2_Pi::Array{Float64}
    k2_Dx::Array{Float64}
    
    k3_Phi::Array{Float64}
    k3_Pi::Array{Float64}
    k3_Dx::Array{Float64}

    function Substeps(num_points)
        new(
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points)
        )
    end
end

function compute_Phi_rhs!(Pi, Phi_rhs)
    copy!(Phi_rhs, Pi)
end

function compute_Pi_rhs!(D, Dx, Pi_rhs)
    mul!(Pi_rhs, D, Dx)
end

function compute_Dx_rhs!(D, Pi, Dx_rhs)
    mul!(Dx_rhs, D, Pi)
end

function compute_rhs!(D, y::GridFuncs, dy::GridFuncs)
    compute_Phi_rhs!(y.Pi, dy.Phi)
    compute_Pi_rhs!(D, y.Dx, dy.Pi)
    compute_Dx_rhs!(D, y.Pi, dy.Dx)
end

function compute_rhs!(D, Phi, Pi, Dx, dy::GridFuncs)
    compute_Phi_rhs!(Phi, dy.Phi)
    compute_Pi_rhs!(D, Dx, dy.Pi)
    compute_Dx_rhs!(D, Pi, dy.Dx)
end

function compute_k0!(D, ks::Substeps, yp::GridFuncs, dy::GridFuncs)
    # Compute the RHS on the previous time step and store it in dy
    compute_rhs!(D, yp, dy)

    # Copy dy to k0
    copy!(ks.k0_Phi, dy.Phi)
    copy!(ks.k0_Pi, dy.Pi)
    copy!(ks.k0_Dx, dy.Dx)
end

function compute_k1!(h, cs, D, ks::Substeps, y::GridFuncs, dy::GridFuncs)
    c0 = cs[1]

    # Step 1: Create the argument of the rhs call by storing it into k1
    @. ks.k1_Phi = y.Phi + h * c0 * ks.k0_Phi
    @. ks.k1_Pi = y.Pi + h * c0 * ks.k0_Pi
    @. ks.k1_Dx = y.Dx+ h * c0 * ks.k0_Dx

    # Step 2: Compute the RHS using the values in ks.k1_[...] as state
    compute_rhs!(D, ks.k1_Phi, ks.k1_Pi, ks.k1_Dx, dy)

    # Step 3: Copy the results in dy back to ks.k1_[...]
    copy!(ks.k1_Phi, dy.Phi)
    copy!(ks.k1_Pi, dy.Pi)
    copy!(ks.k1_Dx, dy.Dx)
end

function compute_k2!(h, cs, D, ks::Substeps, y::GridFuncs, dy::GridFuncs)
    c1 = cs[2]
    c2 = cs[3]

    # Step 1: Create the argument of the rhs call by storing it into k2
    @. ks.k2_Phi = y.Phi + h * (c1 * ks.k1_Phi + c2 * ks.k0_Phi)
    @. ks.k2_Pi = y.Pi + h * (c1 * ks.k1_Pi + c2 * ks.k0_Pi)
    @. ks.k2_Dx = y.Dx + h * (c1 * ks.k1_Dx + c2 * ks.k0_Dx)

    # Step 2: Compute the RHS using the values in ks.k2_[...] as state
    compute_rhs!(D, ks.k2_Phi, ks.k2_Pi, ks.k2_Dx, dy)

    # Step 3: Copy the results in dy back to ks.k2_[...]
    copy!(ks.k2_Phi, dy.Phi)
    copy!(ks.k2_Pi, dy.Pi)
    copy!(ks.k2_Dx, dy.Dx)
end

function compute_k3!(h, cs, D, ks::Substeps, y::GridFuncs, dy::GridFuncs)
    c3 = cs[4]
    c4 = cs[5]
    c5 = cs[6]

    # Step 1: Create the argument of the rhs call by storing it into k3
    @. ks.k3_Phi = y.Phi + h * (c3 * ks.k2_Phi + c4 * ks.k1_Phi + c5 * ks.k0_Phi)
    @. ks.k3_Pi = y.Pi + h * (c3 * ks.k2_Pi + c4 * ks.k1_Pi + c5 * ks.k0_Pi)
    @. ks.k3_Dx = y.Dx + h * (c3 * ks.k2_Dx + c4 * ks.k1_Dx + c5 * ks.k0_Dx)

    # Step 2: Compute the RHS using the values in ks.k3_[...] as state
    compute_rhs!(D, ks.k3_Phi, ks.k3_Pi, ks.k3_Dx, dy)

    # Step 3: Copy the results in dy back to ks.k3_[...]
    copy!(ks.k3_Phi, dy.Phi)
    copy!(ks.k3_Pi, dy.Pi)
    copy!(ks.k3_Dx, dy.Dx)
end

function rkab_step!(h, cs, D, ks::Substeps, yp::GridFuncs, y::GridFuncs, dy::GridFuncs)
    c6 = cs[7]
    c7 = cs[8]
    c8 = cs[9]
    c9 = cs[10]
    
    # Step 1: Compute ks
    compute_k0!(D, ks, yp, dy)
    compute_k1!(h, cs, D, ks, y, dy)
    compute_k2!(h, cs, D, ks, y, dy)
    compute_k3!(h, cs, D, ks, y, dy)

    # Step 2: Store the current state vector as the previous state vector
    copy!(yp.Phi, y.Phi)
    copy!(yp.Pi, y.Pi)
    copy!(yp.Dx, y.Dx)

    # Step 3: Apply the evolution formula
    @. y.Phi =  y.Phi + h * (c6 * ks.k0_Phi + c7 * ks.k1_Phi + c8 * ks.k2_Phi + c9 * ks.k3_Phi)
    @. y.Pi =  y.Pi + h * (c6 * ks.k0_Pi + c7 * ks.k1_Pi + c8 * ks.k2_Pi + c9 * ks.k3_Pi)
    @. y.Dx =  y.Dx + h * (c6 * ks.k0_Dx + c7 * ks.k1_Dx + c8 * ks.k2_Dx + c9 * ks.k3_Dx)
end

function euler_step!(h, D, y::GridFuncs, dy::GridFuncs)
    # Step 1: Compute RHS
    compute_rhs!(D, y::GridFuncs, dy::GridFuncs)

    # Step 2: Apply the evolution formula
    @. y.Phi = y.Phi + h * dy.Phi
    @. y.Pi = y.Pi + h * dy.Pi
    @. y.Dx = y.Dx + h * dy.Dx
end

function apply_dirichlet_bcs!(A, kx, t, x0, xf, y::GridFuncs)
    y.Phi[1] = sw_Phi(A, kx, t, x0)
    y.Pi[1] = sw_Pi(A, kx, t, x0)
    y.Dx[1] = sw_Dx(A, kx, t, x0)

    y.Phi[end] = sw_Phi(A, kx, t, xf)
    y.Pi[end] = sw_Pi(A, kx, t, xf)
    y.Dx[end] = sw_Dx(A, kx, t, xf)
end

function write_state(group, i, y::GridFuncs)
    write(group, "Phi_$(lpad(i, 4, "0"))", y.Phi)
    write(group, "Pi_$(lpad(i, 4, "0"))", y.Pi)
    write(group, "Dx_$(lpad(i, 4, "0"))", y.Dx)
end

function write_rhs(group, i, dy::GridFuncs)
    write(group, "Phi_rhs_$(lpad(i, 4, "0"))", dy.Phi)
    write(group, "Pi_rhs_$(lpad(i, 4, "0"))", dy.Pi)
    write(group, "Dx_rhs_$(lpad(i, 4, "0"))", dy.Dx)
end

function evolve()
    # Domain
    x0 = -1.0
    xf = 1.0
    num_pts = 81
    dx = (xf - x0) / (num_pts - 1)

    # Time steps
    cfl = 0.25
    dt = cfl * dx
    last_iter = 160

    # Standing wave params
    A = 1.0
    kx = 1.0

    # RKAB cs
    c0 = 0
    c1 = 0.3736646857963324
    c2 = 0.03127973625120939
    c3 = 0.6231453199819736
    c4 = -1.322215937221026
    c5 = -0.25177800477107587
    c6 = -3.4674940572428117
    c7 = -1.1964762610986783
    c8 = 1.7835202250496323
    c9 = 1 - c6 - c7 - c8
    
    cs = [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9]

    # IO
    h5_file = h5open("output.h5", "w")
    state_group = create_group(h5_file, "state")
    rhs_group = create_group(h5_file, "rhs")
    grid_group = create_group(h5_file, "grid")

    # Derivative operator
    D = derivative_operator(
        DienerDorbandSchnetterTiglio2007(),
        derivative_order=1,
        accuracy_order=4,
        xmin=x0,
        xmax=xf,
        N=num_pts
    )

    # GridFuncs
    y = GridFuncs(num_pts)
    
    # Previous state
    yp = GridFuncs(num_pts)
    
    # RHS
    dy = GridFuncs(num_pts)

    # RKAB Ks (k0 .. k3 for each variable)
    ks = Substeps(num_pts)
    
    # Init state and previous state
    for i in 0:(num_pts - 1)
        x = x0 + i * dx
        y.Phi[i + 1] = sw_Phi(A, kx, 0.0, x)
        y.Pi[i + 1] = sw_Pi(A, kx, 0.0, x)
        y.Dx[i + 1] = sw_Dx(A, kx, 0.0, x)

        yp.Phi[i + 1] = sw_Phi(A, kx, -dt, x)
        yp.Pi[i + 1] = sw_Pi(A, kx, -dt, x)
        yp.Dx[i + 1] = sw_Dx(A, kx, -dt, x)
    end

    # Write grid and initial state
    write(grid_group, "x_coords", collect(grid(D)))
    write_state(state_group, 0, y)
    write_rhs(rhs_group, 0, dy)

    # Iterate
    for i in 1:last_iter
        t = i * dt
        println("Iteration $i, t = $t")
        println("  Stepping")
        #rkab_step!(dt, cs, D, ks, yp, y, dy)
        euler_step!(dt, D, y, dy)

        println("  Applying BCs")
        apply_dirichlet_bcs!(A, kx, t, x0, xf, y)

        println("  Saving")
        write_state(state_group, i, y)
        write_rhs(rhs_group, i, dy)
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