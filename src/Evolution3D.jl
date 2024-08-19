using LinearAlgebra
using SummationByPartsOperators

function du_dx_3d!(D, u, dudx)
    for j in axes(u, 2)
        for k in axes(u, 3)
            mul!(view(dudx, :, j, k), D, view(u, :, j, k))
        end
    end
end

function du_dy_3d!(D, u, dudy)
    for i in axes(u, 1)
        for k in axes(u, 3)
            mul!(view(dudy, i, :, k), D, view(u, i, :, k))
        end
    end
end

function du_dz_3d!(D, u, dudz)
    for i in axes(u, 1)
        for j in axes(u, 2)
            mul!(view(dudz, i, j, :), D, view(u, i, j, :))
        end
    end
end

function compute_Pi_rhs!(D::DerivativeOperator, d::Derivatives3D, Dx, Dy, Dz, Pi_rhs)
    du_dx_3d!(D, Dx, d.dDx_dx)
    du_dy_3d!(D, Dy, d.dDy_dy)
    du_dz_3d!(D, Dz, d.dDz_dz)
    @. Pi_rhs = d.dDx_dx + d.dDy_dy + d.dDz_dz
end

function compute_rhs!(D::DerivativeOperator, d::Derivatives3D, y::GridFuncs3D, dy::GridFuncs3D)
    compute_Phi_rhs!(y.Pi, dy.Phi)
    compute_Pi_rhs!(D, d, y.Dx, y.Dy, y.Dz, dy.Pi)
    du_dx_3d!(D, y.Pi, dy.Dx)
    du_dy_3d!(D, y.Pi, dy.Dy)
    du_dz_3d!(D, y.Pi, dy.Dz)
end

function compute_rhs!(D::DerivativeOperator, d::Derivatives3D, Pi, Dx, Dy, Dz, dy::GridFuncs3D)
    compute_Phi_rhs!(Pi, dy.Phi)
    compute_Pi_rhs!(D, d, Dx, Dy, Dz, dy.Pi)
    du_dx_3d!(D, Pi, Dx)
    du_dy_3d!(D, Pi, Dy)
    du_dz_3d!(D, Pi, Dz)
end

function compute_k0!(D::DerivativeOperator, ks::Substeps3D, d::Derivatives3D, yp::GridFuncs3D, dy::GridFuncs3D)
    # Compute the RHS on the previous time step and store it in dy
    compute_rhs!(D, d, yp, dy)

    # Copy dy to k0
    copy!(ks.k0_Phi, dy.Phi)
    copy!(ks.k0_Pi, dy.Pi)
    copy!(ks.k0_Dx, dy.Dx)
    copy!(ks.k0_Dy, dy.Dy)
    copy!(ks.k0_Dz, dy.Dz)
end

function compute_k1!(h, cs, D::DerivativeOperator, ks::Substeps3D, d::Derivatives3D, y::GridFuncs3D, dy::GridFuncs3D)
    c0 = cs[1]

    # Step 1: Create the argument of the rhs call by storing it into k1
    @. ks.k1_Phi = y.Phi + h * c0 * ks.k0_Phi
    @. ks.k1_Pi = y.Pi + h * c0 * ks.k0_Pi
    @. ks.k1_Dx = y.Dx + h * c0 * ks.k0_Dx
    @. ks.k1_Dy = y.Dy + h * c0 * ks.k0_Dy
    @. ks.k1_Dz = y.Dz + h * c0 * ks.k0_Dz

    # Step 2: Compute the RHS using the values in ks.k1_[...] as state
    compute_rhs!(
        D,
        d,
        ks.k1_Pi,
        ks.k1_Dx,
        ks.k1_Dy,
        ks.k1_Dz,
        dy
    )

    # Step 3: Copy the results in dy back to ks.k1_[...]
    copy!(ks.k1_Phi, dy.Phi)
    copy!(ks.k1_Pi, dy.Pi)
    copy!(ks.k1_Dx, dy.Dx)
    copy!(ks.k1_Dy, dy.Dy)
    copy!(ks.k1_Dz, dy.Dz)
end

function compute_k2!(h, cs, D::DerivativeOperator, ks::Substeps3D, d::Derivatives3D, y::GridFuncs3D, dy::GridFuncs3D)
    c1 = cs[2]
    c2 = cs[3]

    # Step 1: Create the argument of the rhs call by storing it into k2
    @. ks.k2_Phi = y.Phi + h * (c1 * ks.k1_Phi + c2 * ks.k0_Phi)
    @. ks.k2_Pi = y.Pi + h * (c1 * ks.k1_Pi + c2 * ks.k0_Pi)
    @. ks.k2_Dx =  y.Dx + h * (c1 * ks.k1_Dx + c2 * ks.k0_Dx)
    @. ks.k2_Dy =  y.Dy + h * (c1 * ks.k1_Dy + c2 * ks.k0_Dy)
    @. ks.k2_Dz =  y.Dz + h * (c1 * ks.k1_Dz + c2 * ks.k0_Dz)

    # Step 2: Compute the RHS using the values in ks.k2_[...] as state
    compute_rhs!(
        D,
        d,
        ks.k2_Pi,
        ks.k2_Dx,
        ks.k2_Dy,
        ks.k2_Dz,
        dy
    )

    # Step 3: Copy the results in dy back to ks.k2_[...]
    copy!(ks.k2_Phi, dy.Phi)
    copy!(ks.k2_Pi, dy.Pi)
    copy!(ks.k2_Dx, dy.Dx)
    copy!(ks.k2_Dy, dy.Dy)
    copy!(ks.k2_Dz, dy.Dz)
end

function compute_k3!(h, cs, D::DerivativeOperator, ks::Substeps3D, d::Derivatives3D, y::GridFuncs3D, dy::GridFuncs3D)
    c3 = cs[4]
    c4 = cs[5]
    c5 = cs[6]

    # Step 1: Create the argument of the rhs call by storing it into k3
    @. ks.k3_Phi = y.Phi + h * (c3 * ks.k2_Phi + c4 * ks.k1_Phi + c5 * ks.k0_Phi)
    @. ks.k3_Pi = y.Pi + h * (c3 * ks.k2_Pi + c4 * ks.k1_Pi + c5 * ks.k0_Pi)
    @. ks.k3_Dx = y.Dx + h * (c3 * ks.k2_Dx + c4 * ks.k1_Dx + c5 * ks.k0_Dx)
    @. ks.k3_Dy = y.Dy + h * (c3 * ks.k2_Dy + c4 * ks.k1_Dy + c5 * ks.k0_Dy)
    @. ks.k3_Dz = y.Dz + h * (c3 * ks.k2_Dz + c4 * ks.k1_Dz + c5 * ks.k0_Dz)

    # Step 2: Compute the RHS using the values in ks.k3_[...] as state
    compute_rhs!(
        D,
        d,
        ks.k3_Pi,
        ks.k3_Dx,
        ks.k3_Dy,
        ks.k3_Dz,
        dy
    )

    # Step 3: Copy the results in dy back to ks.k2_[...]
    copy!(ks.k3_Phi, dy.Phi)
    copy!(ks.k3_Pi, dy.Pi)
    copy!(ks.k3_Dx, dy.Dx)
    copy!(ks.k3_Dy, dy.Dy)
    copy!(ks.k3_Dz, dy.Dz)
end

function rkab_step!(h, cs, D::DerivativeOperator, ks::Substeps3D, d::Derivatives3D, yp::GridFuncs3D, y::GridFuncs3D, dy::GridFuncs3D)
    c6 = cs[7]
    c7 = cs[8]
    c8 = cs[9]
    c9 = cs[10]
    
    # Step 1: Compute ks
    compute_k0!(D, ks, d, yp, dy)
    compute_k1!(h, cs, D, ks, d, y, dy)
    compute_k2!(h, cs, D, ks, d, y, dy)
    compute_k3!(h, cs, D, ks, d, y, dy)

    # Step 2: Store the current state vector as the previous state vector
    copy!(yp.Phi, y.Phi)
    copy!(yp.Pi, y.Pi)
    copy!(yp.Dx, y.Dx)
    copy!(yp.Dy, y.Dy)
    copy!(yp.Dz, y.Dz)

    # Step 3: Apply the evolution formula
    @. y.Phi =  y.Phi + h * (c6 * ks.k0_Phi + c7 * ks.k1_Phi + c8 * ks.k2_Phi + c9 * ks.k3_Phi)
    @. y.Pi =  y.Pi + h * (c6 * ks.k0_Pi + c7 * ks.k1_Pi + c8 * ks.k2_Pi + c9 * ks.k3_Pi)
    @. y.Dx =  y.Dx + h * (c6 * ks.k0_Dx + c7 * ks.k1_Dx + c8 * ks.k2_Dx + c9 * ks.k3_Dx)
    @. y.Dy =  y.Dy + h * (c6 * ks.k0_Dy + c7 * ks.k1_Dy + c8 * ks.k2_Dy + c9 * ks.k3_Dy)
    @. y.Dz =  y.Dz + h * (c6 * ks.k0_Dz + c7 * ks.k1_Dz + c8 * ks.k2_Dz + c9 * ks.k3_Dz)
end

function euler_step!(h, D::DerivativeOperator, d::Derivatives3D, y::GridFuncs3D, dy::GridFuncs3D)
    # Step 1: Compute RHS
    compute_rhs!(D, d, y, dy)

    # Step 2: Apply the evolution formula
    @. y.Phi = y.Phi + h * dy.Phi
    @. y.Pi = y.Pi + h * dy.Pi
    @. y.Dx = y.Dx + h * dy.Dx
    @. y.Dy = y.Dy + h * dy.Dy
    @. y.Dz = y.Dz + h * dy.Dz
end

function apply_dirichlet_bcs!(A, kx, ky, kz, t, r0, dr, num_pts, y::GridFuncs3D)
    for i in 0:(num_pts - 1)
        for j in 0:(num_pts - 1)
            for k in 0:(num_pts - 1)
                x_bnd = i == 0 || i == (num_pts - 1)
                y_bnd = j == 0 || j == (num_pts - 1)
                z_bnd = k == 0 || k == (num_pts - 1)
                
                if x_bnd || y_bnd || z_bnd
                    X = r0 + i * dr
                    Y = r0 + j * dr
                    Z = r0 + k * dr

                    y.Phi[i + 1, j + 1, k + 1] = sw_Phi(A, kx, ky, kz, t, X, Y, Z)
                    y.Pi[i + 1, j + 1, k + 1] = sw_Pi(A, kx, ky, kz, t, X, Y, Z)
                    y.Dx[i + 1, j + 1, k + 1] = sw_Dx(A, kx, ky, kz, t, X, Y, Z)
                    y.Dy[i + 1, j + 1, k + 1] = sw_Dy(A, kx, ky, kz, t, X, Y, Z)
                    y.Dz[i + 1, j + 1, k + 1] = sw_Dz(A, kx, ky, kz, t, X, Y, Z)
                end
            end
        end
    end
end