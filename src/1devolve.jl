using RKAB

using HDF5
using LinearAlgebra
using SummationByPartsOperators

function main()
    if size(ARGS, 1) == 1
        config_file = ARGS[1]
    else
        config_file = "param_1d.yaml"
    end
    
    @info "Using parameter file $config_file"

    config_data = Params1D(config_file)

    # Domain
    x0 = config_data.domain_start
    xf = config_data.domain_end
    num_pts = config_data.domain_points
    dx = (xf - x0) / (num_pts - 1)
    
    # Time steps
    final_time = config_data.time_final
    cfl = config_data.time_cfl
    dt = cfl * dx
    last_iter = floor(Int, final_time / dt) + 1

    @info "Spatial step = $dx"
    @info "Time step = $dt"
    @info "Total iterations: $last_iter"

    # Standing wave params
    A = config_data.standing_wave_A
    kx = config_data.standing_wave_kx

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
        
        @info "Iteration $i, t = $t"
        @info "  Stepping"
        rkab_step!(dt, cs, D, ks, yp, y, dy)

        @info "  Applying BCs"
        apply_dirichlet_bcs!(A, kx, t, x0, xf, y)

        @info "  Saving"
        write_state(state_group, i, y)
        write_rhs(rhs_group, i, dy)
     end

    close(h5_file)
end

main()