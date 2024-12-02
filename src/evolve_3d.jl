using RKAB

using HDF5
using LinearAlgebra
using SummationByPartsOperators
using Random

function evolve_3d(config_file)
    @info "Using parameter file $config_file"

    config_data = Params3D(config_file)

    # Domain
    r0 = config_data.domain_start
    rf = config_data.domain_end
    
    config_num_pts = config_data.domain_points
    
    if config_data.boundary_type == "periodic"
        num_pts = config_num_pts
        dr = (rf - r0) / num_pts
    elseif config_data.boundary_type == "dirichlet"
        num_pts = config_num_pts + 1
        dr = (rf - r0) / config_num_pts
    else
         @error "Unrecognized boundary type $(config_data.derivatives_type)"
    end
    
    # Time steps
    final_time = config_data.time_final
    cfl = config_data.time_cfl
    dt = cfl * dr
    last_iter = floor(Int, final_time / dt) + 1

    @info "Spatial step = $dr"
    @info "Time step = $dt"
    @info "Total iterations: $last_iter"
    @info "Save every $(config_data.output_every) iterations"

    # IO
    h5_file = h5open(config_data.output_file, "w")
    state_group = create_group(h5_file, "state")
    rhs_group = create_group(h5_file, "rhs")
    grid_group = create_group(h5_file, "grid")

    # Derivative operator
    if config_data.boundary_type == "periodic"
        D = periodic_derivative_operator(
            derivative_order=1,
            accuracy_order=config_data.derivative_order,
            xmin=r0,
            xmax=rf,
            N=num_pts
        )
    else
        D = derivative_operator(
            DienerDorbandSchnetterTiglio2007(),
            derivative_order=1,
            accuracy_order=config_data.derivative_order,
            xmin=r0,
            xmax=rf,
            N=num_pts
        )
    end
    
    @info "Derivative operators used: $D"
    
    # Grids
    grids = grid(D)
    grids_array = collect(grids)
    @info "Grids used: $grids"

    # GridFuncs
    y = GridFuncs3D(num_pts)
    
    # Previous state
    yp = GridFuncs3D(num_pts)
    
    # RHS
    dy = GridFuncs3D(num_pts)

    # RKAB cs    
    cs = config_data.RKAB_coeffs

    # RKAB Ks (k0 .. k3 for each variable)
    ks = Substeps3D(num_pts)

    # Scratch space for derivative values
    d = Derivatives3D(num_pts)

    # Init state and previous state
    @info "Initializing state vector and previous state vector"
    if config_data.id_type == "standing"
        A = config_data.standing_wave_A
        kx = config_data.standing_wave_kx
        ky = config_data.standing_wave_ky
        kz = config_data.standing_wave_kz

        attributes(h5_file)["id_type"] = "standing"
        attributes(h5_file)["A"] = A
        attributes(h5_file)["kx"] = kx
        attributes(h5_file)["ky"] = ky
        attributes(h5_file)["kz"] = kz
        
        for i in 0:(num_pts - 1)
            for j in 0:(num_pts - 1)
                for k in 0:(num_pts - 1)
                    X = r0 + i * dr
                    Y = r0 + j * dr
                    Z = r0 + k * dr
                    
                    y.Phi[i + 1, j + 1, k + 1] = sw_Phi(A, kx, ky, kz, 0.0, X, Y, Z)
                    y.Pi[i + 1, j + 1, k + 1] = sw_Pi(A, kx, ky, kz, 0.0, X, Y, Z)
                    y.Dx[i + 1, j + 1, k + 1] = sw_Dx(A, kx, ky, kz, 0.0, X, Y, Z)
                    y.Dy[i + 1, j + 1, k + 1] = sw_Dy(A, kx, ky, kz, 0.0, X, Y, Z)
                    y.Dz[i + 1, j + 1, k + 1] = sw_Dz(A, kx, ky, kz, 0.0, X, Y, Z)

                    yp.Phi[i + 1, j + 1, k + 1] = sw_Phi(A, kx, ky, kz, -dt, X, Y, Z)
                    yp.Pi[i + 1, j + 1, k + 1] = sw_Pi(A, kx, ky, kz, -dt, X, Y, Z)
                    yp.Dx[i + 1, j + 1, k + 1] = sw_Dx(A, kx, ky, kz, -dt, X, Y, Z)
                    yp.Dy[i + 1, j + 1, k + 1] = sw_Dy(A, kx, ky, kz, -dt, X, Y, Z)
                    yp.Dz[i + 1, j + 1, k + 1] = sw_Dz(A, kx, ky, kz, -dt, X, Y, Z)
                end
            end
        end
    elseif config_data.id_type == "noise"
        seed = config_data.noise_seed
        range = config_data.noise_range

        attributes(h5_file)["id_type"] = "noise"
        attributes(h5_file)["seed"] = seed
        attributes(h5_file)["range"] = range
        
        for i in 0:(num_pts - 1)
            for j in 0:(num_pts - 1)
                for k in 0:(num_pts - 1)
                    X = r0 + i * dr
                    Y = r0 + j * dr
                    Z = r0 + k * dr

                    noise_val = range * sin(seed * X / dr) * sin(seed * Y / dr) * sin(seed * Z / dr)

                    y.Phi[i + 1, j + 1, k + 1] = noise_val
                    y.Pi[i + 1, j + 1, k + 1] = noise_val
                    y.Dx[i + 1, j + 1, k + 1] = noise_val
                    y.Dy[i + 1, j + 1, k + 1] = noise_val
                    y.Dz[i + 1, j + 1, k + 1] = noise_val

                    yp.Phi[i + 1, j + 1, k + 1] = noise_val
                    yp.Pi[i + 1, j + 1, k + 1] = noise_val
                    yp.Dx[i + 1, j + 1, k + 1] = noise_val
                    yp.Dy[i + 1, j + 1, k + 1] = noise_val
                    yp.Dz[i + 1, j + 1, k + 1] = noise_val
                end
            end
        end
    else
        @error "Unrecognized initial data $(config_data.id_type)"
    end

    # Write grid and initial state
    attributes(h5_file)["final_time"] = final_time
    attributes(h5_file)["cfl"] = cfl
    attributes(h5_file)["last_iter"] = last_iter
    attributes(h5_file)["dt"] = dt
    attributes(h5_file)["time_stepping_method"] = config_data.time_method
    
    attributes(h5_file)["derivative_order"] = config_data.derivative_order
    attributes(h5_file)["boundary_type"] = config_data.boundary_type
    attributes(h5_file)["dr"] = dr
    
    write(grid_group, "x_coords", grids_array)
    write(grid_group, "y_coords", grids_array)
    write(grid_group, "z_coords", grids_array)
    
    write_state(state_group, 0, y)
    if config_data.output_rhs
        write_rhs(rhs_group, 0, dy)
    end

    out_every_counter = 1

    for i in 1:last_iter
        t = i * dt
        
        # Step
        @info "Iteration $i, t = $t"
        @info "  Stepping"
        if config_data.time_method == "RKAB"
            #rkab_step!(dt, cs, D, ks, d, yp, y, dy)
            rkab2_step!(dt, cs, D, ks, d, yp, y, dy)
        elseif config_data.time_method == "Euler"
            euler_step!(dt, D, d, y, dy)
        else
            @error "Unrecognized time stepping method $(config_data.time_method)"
        end
        
        # BCs    
        if config_data.boundary_type == "dirichlet"
            @info "  Applying BCs"
            apply_dirichlet_bcs!(A, kx, ky, kz, t, r0, dr, num_pts, y)
        end

        # Save
        if out_every_counter == config_data.output_every
            @info "  Saving"
            write_state(state_group, i, y)
            
            if config_data.output_rhs
                write_rhs(rhs_group, i, dy)
            end
            
            out_every_counter = 1
        else
            out_every_counter += 1
        end
    end

    close(h5_file)
end