using LinearAlgebra
using SummationByPartsOperators

function du_dx!(D, u, dudx)
    for j in axes(u, 2)
        mul!(view(dudx, :, j), D, view(u, :, j))
    end
end

function du_dy!(D, u, dudy)
    for i in axes(u, 1)
        mul!(view(dudy, i, :), D, view(u, i, :))
    end
end

function compute_Pi_rhs!(D::DerivativeOperator, d::Derivatives2D, Dx, Dy, Pi_rhs)
    du_dx!(D, Dx, d.dDx_dx)
    du_dy!(D, Dy, d.dDy_dy)
    @. Pi_rhs = d.dDx_dx + d.dDy_dy
end

function compute_rhs!(D::DerivativeOperator, d::Derivatives2D, y::GridFuncs2D, dy::GridFuncs2D)
    compute_Phi_rhs!(y.Pi, dy.Phi)
    compute_Pi_rhs!(D, d, y.Dx, y.Dy, dy.Pi)
    du_dx!(D, y.Pi, dy.Dx)
    du_dy!(D, y.Pi, dy.Dy)
end

function compute_rhs!(D::DerivativeOperator, d::Derivatives2D, Pi, Dx, Dy, dy::GridFuncs2D)
    compute_Phi_rhs!(Pi, dy.Phi)
    compute_Pi_rhs!(D, d, Dx, Dy, dy.Pi)
    du_dx!(D, Pi, dy.Dx)
    du_dy!(D, Pi, dy.Dy)
end

function compute_k0!(D::DerivativeOperator, ks::Substeps2D, d::Derivatives2D, yp::GridFuncs2D, dy::GridFuncs2D)
    # Compute the RHS on the previous time step and store it in dy
    compute_rhs!(D, d, yp, dy)

    # Copy dy to k0
    copy!(ks.k0_Phi, dy.Phi)
    copy!(ks.k0_Pi, dy.Pi)
    copy!(ks.k0_Dx, dy.Dx)
    copy!(ks.k0_Dy, dy.Dy)
end

function compute_k1!(h, cs, D::DerivativeOperator, ks::Substeps2D, d::Derivatives2D, y::GridFuncs2D, dy::GridFuncs2D)
    c0 = cs[1]

    # Step 1: Create the argument of the rhs call by storing it into k1
    @. ks.k1_Phi = y.Phi + h * c0 * ks.k0_Phi
    @. ks.k1_Pi = y.Pi + h * c0 * ks.k0_Pi
    @. ks.k1_Dx = y.Dx + h * c0 * ks.k0_Dx
    @. ks.k1_Dy = y.Dy + h * c0 * ks.k0_Dy

    # Step 2: Compute the RHS using the values in ks.k1_[...] as state
    compute_rhs!(
        D,
        d,
        ks.k1_Pi,
        ks.k1_Dx,
        ks.k1_Dy,
        dy
    )

    # Step 3: Copy the results in dy back to ks.k1_[...]
    copy!(ks.k1_Phi, dy.Phi)
    copy!(ks.k1_Pi, dy.Pi)
    copy!(ks.k1_Dx, dy.Dx)
    copy!(ks.k1_Dx, dy.Dy)
end

function compute_k2!(h, cs, D::DerivativeOperator, ks::Substeps2D, d::Derivatives2D, y::GridFuncs2D, dy::GridFuncs2D)
    c1 = cs[2]
    c2 = cs[3]

    # Step 1: Create the argument of the rhs call by storing it into k2
    @. ks.k2_Phi = y.Phi + h * (c1 * ks.k1_Phi + c2 * ks.k0_Phi)
    @. ks.k2_Pi = y.Pi + h * (c1 * ks.k1_Pi + c2 * ks.k0_Pi)
    @. ks.k2_Dx =  y.Dx + h * (c1 * ks.k1_Dx + c2 * ks.k0_Dx)
    @. ks.k2_Dy =  y.Dy + h * (c1 * ks.k1_Dy + c2 * ks.k0_Dy)

    # Step 2: Compute the RHS using the values in ks.k2_[...] as state
    compute_rhs!(
        D,
        d,
        ks.k2_Pi,
        ks.k2_Dx,
        ks.k2_Dy,
        dy
    )

    # Step 3: Copy the results in dy back to ks.k2_[...]
    copy!(ks.k2_Phi, dy.Phi)
    copy!(ks.k2_Pi, dy.Pi)
    copy!(ks.k2_Dx, dy.Dx)
    copy!(ks.k2_Dx, dy.Dy)
end

function compute_k3!(h, cs, D::DerivativeOperator, ks::Substeps2D, d::Derivatives2D, y::GridFuncs2D, dy::GridFuncs2D)
    c3 = cs[4]
    c4 = cs[5]
    c5 = cs[6]

    # Step 1: Create the argument of the rhs call by storing it into k3
    @. ks.k3_Phi = y.Phi + h * (c3 * ks.k2_Phi + c4 * ks.k1_Phi + c5 * ks.k0_Phi)
    @. ks.k3_Pi = y.Pi + h * (c3 * ks.k2_Pi + c4 * ks.k1_Pi + c5 * ks.k0_Pi)
    @. ks.k3_Dx = y.Dx + h * (c3 * ks.k2_Dx + c4 * ks.k1_Dx + c5 * ks.k0_Dx)
    @. ks.k3_Dy = y.Dy + h * (c3 * ks.k2_Dy + c4 * ks.k1_Dy + c5 * ks.k0_Dy)

    # Step 2: Compute the RHS using the values in ks.k3_[...] as state
    compute_rhs!(
        D,
        d,
        ks.k3_Pi,
        ks.k3_Dx,
        ks.k3_Dy,
        dy
    )

    # Step 3: Copy the results in dy back to ks.k2_[...]
    copy!(ks.k3_Phi, dy.Phi)
    copy!(ks.k3_Pi, dy.Pi)
    copy!(ks.k3_Dx, dy.Dx)
    copy!(ks.k3_Dx, dy.Dy)
end

function rkab_step!(h, cs, D::DerivativeOperator, ks::Substeps2D, d::Derivatives2D, yp::GridFuncs2D, y::GridFuncs2D, dy::GridFuncs2D)
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

    # Step 3: Apply the evolution formula
    @. y.Phi =  y.Phi + h * (c6 * ks.k0_Phi + c7 * ks.k1_Phi + c8 * ks.k2_Phi + c9 * ks.k3_Phi)
    @. y.Pi =  y.Pi + h * (c6 * ks.k0_Pi + c7 * ks.k1_Pi + c8 * ks.k2_Pi + c9 * ks.k3_Pi)
    @. y.Dx =  y.Dx + h * (c6 * ks.k0_Dx + c7 * ks.k1_Dx + c8 * ks.k2_Dx + c9 * ks.k3_Dx)
    @. y.Dy =  y.Dy + h * (c6 * ks.k0_Dy + c7 * ks.k1_Dy + c8 * ks.k2_Dy + c9 * ks.k3_Dy)
end

function euler_step!(h, D::DerivativeOperator, d::Derivatives2D, y::GridFuncs2D, dy::GridFuncs2D)
    # Step 1: Compute RHS
    compute_rhs!(D, d, y, dy)

    # Step 2: Apply the evolution formula
    @. y.Phi = y.Phi + h * dy.Phi
    @. y.Pi = y.Pi + h * dy.Pi
    @. y.Dx = y.Dx + h * dy.Dx
    @. y.Dy = y.Dy + h * dy.Dy
end

function apply_dirichlet_bcs!(A, kx, ky, t, r0, dr, num_pts, y::GridFuncs2D)
    for i in 0:(num_pts - 1)
        for j in 0:(num_pts - 1)
            x_bnd = i == 0 || i == (num_pts - 1)
            y_bnd = j == 0 || j == (num_pts - 1)
            
            if x_bnd || y_bnd
                X = r0 + i * dr
                Y = r0 + j * dr

                y.Phi[i + 1, j + 1] = sw_Phi(A, kx, ky, t, X, Y)
                y.Pi[i + 1, j + 1] = sw_Pi(A, kx, ky, t, X, Y)
                y.Dx[i + 1, j + 1] = sw_Dx(A, kx, ky, t, X, Y)
                y.Dy[i + 1, j + 1] = sw_Dy(A, kx, ky, t, X, Y)
            end
        end
    end
end