using RKAB

using HDF5
using LinearAlgebra
using SummationByPartsOperators

function evolve()
    if size(ARGS, 1) == 1
        config_file = ARGS[1]
    else
        config_file = "param_1d.yaml"
    end
    
    @info "Using parameter file $config_file"

    config_data = Params2D(config_file)

    # Domain
    r0 = config_data.domain_start
    rf = config_data.domain_end
    num_pts = config_data.domain_points
    dr = (rf - r0) / (num_pts - 1)
    
    # Time steps
    final_time = config_data.time_final
    cfl = config_data.time_cfl
    dt = cfl * dr
    last_iter = convert(Int, ceil(final_time / dt))

    @info "Spatial step = $dr"
    @info "Time step = $dt"
    @info "Total iterations: $last_iter"

    # Standing wave params
    A = config_data.standing_wave_A
    kx = config_data.standing_wave_kx
    ky = config_data.standing_wave_ky

    # RKAB cs    
    cs = config_data.RKAB_coeffs

    # IO
    h5_file = h5open(config_data.output_file, "w")
    state_group = create_group(h5_file, "state")
    rhs_group = create_group(h5_file, "rhs")
    grid_group = create_group(h5_file, "grid")

    # Derivative operator
    D = derivative_operator(
        DienerDorbandSchnetterTiglio2007(),
        derivative_order=1,
        accuracy_order=4,
        xmin=r0,
        xmax=rf,
        N=num_pts
    )

    # GridFuncs
    y = GridFuncs2D(num_pts)
    
    # Previous state
    yp = GridFuncs2D(num_pts)
    
    # RHS
    dy = GridFuncs2D(num_pts)

    # RKAB Ks (k0 .. k3 for each variable)
    ks = Substeps2D(num_pts)

    # Scratch space for derivative values
    d = Derivatives2D(num_pts)

    # Init state and previous state
    for i in 0:(num_pts - 1)
        for j in 0:(num_pts - 1)
            X = r0 + i * dr
            Y = r0 + j * dr
            
            y.Phi[i + 1, j + 1] = sw_Phi(A, kx, ky, 0.0, X, Y)
            y.Pi[i + 1, j + 1] = sw_Pi(A, kx, ky, 0.0, X, Y)
            y.Dx[i + 1, j + 1] = sw_Dx(A, kx, ky, 0.0, X, Y)
            y.Dy[i + 1, j + 1] = sw_Dy(A, kx, ky, 0.0, X, Y)

            dy.Phi[i + 1, j + 1] = sw_Phi(A, kx, ky, -dt, X, Y)
            dy.Pi[i + 1, j + 1] = sw_Pi(A, kx, ky, -dt, X, Y)
            dy.Dx[i + 1, j + 1] = sw_Dx(A, kx, ky, -dt, X, Y)
            dy.Dy[i + 1, j + 1] = sw_Dy(A, kx, ky, -dt, X, Y)
        end
    end

    # Write grid and initial state
    attributes(h5_file)["last_iter"]  = last_iter
    attributes(h5_file)["cfl"]  = cfl
    attributes(h5_file)["dt"]  = dt

    attributes(h5_file)["A"]  = A
    attributes(h5_file)["kx"]  = kx
    attributes(h5_file)["ky"]  = ky
    
    write(grid_group, "x_coords", collect(grid(D)))
    write(grid_group, "y_coords", collect(grid(D)))
    write_state(state_group, 0, y)
    write_rhs(rhs_group, 0, dy)

    for i in 1:last_iter
        t = i * dt
        
        @info "Iteration $i, t = $t"
        @info "  Stepping"
        rkab_step!(dt, cs, D, ks, d, yp, y, dy)
        #euler_step!(dt, D, d, y, dy)

        @info "  Applying BCs"
        apply_dirichlet_bcs!(A, kx, ky, t, r0, dr, num_pts, y)

        @info "  Saving"
        write_state(state_group, i, y)
        write_rhs(rhs_group, i, dy)
    end

    close(h5_file)
end

evolve()