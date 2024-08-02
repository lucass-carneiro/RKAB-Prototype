using RKAB

using HDF5
using LinearAlgebra
using SummationByPartsOperators

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
    c0 = 0.0
    c1 = 0.3736646857963324
    c2 = 0.03127973625120939
    c3 = -0.14797683066152537
    c4 = 0.33238257148754524
    c5 = -0.0010981891892632696
    c6 = -0.0547559191353386
    c7 = 2.754535159970365
    c8 = 3.414713672966062
    c9 = 1 - c6 - c7 - c8
    
    cs = [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9]

    # IO
    h5_file = h5open("1d_output.h5", "w")
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
    attributes(h5_file)["last_iter"]  = last_iter
    attributes(h5_file)["cfl"]  = cfl
    attributes(h5_file)["dt"]  = dt

    attributes(h5_file)["A"]  = A
    attributes(h5_file)["kx"]  = kx
    
    write(grid_group, "x_coords", collect(grid(D)))
    write_state(state_group, 0, y)
    write_rhs(rhs_group, 0, dy)

    # Iterate
    for i in 1:last_iter
        t = i * dt
        println("Iteration $i, t = $t")
        println("  Stepping")
        rkab_step!(dt, cs, D, ks, yp, y, dy)
        #euler_step!(dt, D, y, dy)

        println("  Applying BCs")
        apply_dirichlet_bcs!(A, kx, t, x0, xf, y)

        println("  Saving")
        write_state(state_group, i, y)
        write_rhs(rhs_group, i, dy)
     end

    close(h5_file)
end

evolve()