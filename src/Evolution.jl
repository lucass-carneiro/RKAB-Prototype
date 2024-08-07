using LinearAlgebra
using SummationByPartsOperators

function compute_Phi_rhs!(Pi, Phi_rhs)
    copy!(Phi_rhs, Pi)
end

function compute_Pi_rhs!(D::DerivativeOperator, Dx, Pi_rhs)
    mul!(Pi_rhs, D, Dx)
end

function compute_Dx_rhs!(D::DerivativeOperator, Pi, Dx_rhs)
    mul!(Dx_rhs, D, Pi)
end

function compute_rhs!(D::DerivativeOperator, y::GridFuncs, dy::GridFuncs)
    compute_Phi_rhs!(y.Pi, dy.Phi)
    compute_Pi_rhs!(D, y.Dx, dy.Pi)
    compute_Dx_rhs!(D, y.Pi, dy.Dx)
end

function compute_rhs!(D::DerivativeOperator, Pi, Dx, dy::GridFuncs)
    compute_Phi_rhs!(Pi, dy.Phi)
    compute_Pi_rhs!(D, Dx, dy.Pi)
    compute_Dx_rhs!(D, Pi, dy.Dx)
end

function compute_k0!(D::DerivativeOperator, ks::Substeps, yp::GridFuncs, dy::GridFuncs)
    # Compute the RHS on the previous time step and store it in dy
    compute_rhs!(D, yp, dy)

    # Copy dy to k0
    copy!(ks.k0_Phi, dy.Phi)
    copy!(ks.k0_Pi, dy.Pi)
    copy!(ks.k0_Dx, dy.Dx)
end

function compute_k1!(h, cs, D::DerivativeOperator, ks::Substeps, y::GridFuncs, dy::GridFuncs)
    c0 = cs[1]

    # Step 1: Create the argument of the rhs call by storing it into k1
    @. ks.k1_Phi = y.Phi + h * c0 * ks.k0_Phi
    @. ks.k1_Pi = y.Pi + h * c0 * ks.k0_Pi
    @. ks.k1_Dx = y.Dx + h * c0 * ks.k0_Dx

    # Step 2: Compute the RHS using the values in ks.k1_[...] as state
    compute_rhs!(D, ks.k1_Pi, ks.k1_Dx, dy)

    # Step 3: Copy the results in dy back to ks.k1_[...]
    copy!(ks.k1_Phi, dy.Phi)
    copy!(ks.k1_Pi, dy.Pi)
    copy!(ks.k1_Dx, dy.Dx)
end

function compute_k2!(h, cs, D::DerivativeOperator, ks::Substeps, y::GridFuncs, dy::GridFuncs)
    c1 = cs[2]
    c2 = cs[3]

    # Step 1: Create the argument of the rhs call by storing it into k2
    @. ks.k2_Phi = y.Phi + h * (c1 * ks.k1_Phi + c2 * ks.k0_Phi)
    @. ks.k2_Pi = y.Pi + h * (c1 * ks.k1_Pi + c2 * ks.k0_Pi)
    @. ks.k2_Dx = y.Dx + h * (c1 * ks.k1_Dx + c2 * ks.k0_Dx)

    # Step 2: Compute the RHS using the values in ks.k2_[...] as state
    compute_rhs!(D, ks.k2_Pi, ks.k2_Dx, dy)

    # Step 3: Copy the results in dy back to ks.k2_[...]
    copy!(ks.k2_Phi, dy.Phi)
    copy!(ks.k2_Pi, dy.Pi)
    copy!(ks.k2_Dx, dy.Dx)
end

function compute_k3!(h, cs, D::DerivativeOperator, ks::Substeps, y::GridFuncs, dy::GridFuncs)
    c3 = cs[4]
    c4 = cs[5]
    c5 = cs[6]

    # Step 1: Create the argument of the rhs call by storing it into k3
    @. ks.k3_Phi = y.Phi + h * (c3 * ks.k2_Phi + c4 * ks.k1_Phi + c5 * ks.k0_Phi)
    @. ks.k3_Pi = y.Pi + h * (c3 * ks.k2_Pi + c4 * ks.k1_Pi + c5 * ks.k0_Pi)
    @. ks.k3_Dx = y.Dx + h * (c3 * ks.k2_Dx + c4 * ks.k1_Dx + c5 * ks.k0_Dx)

    # Step 2: Compute the RHS using the values in ks.k3_[...] as state
    compute_rhs!(D, ks.k3_Pi, ks.k3_Dx, dy)

    # Step 3: Copy the results in dy back to ks.k3_[...]
    copy!(ks.k3_Phi, dy.Phi)
    copy!(ks.k3_Pi, dy.Pi)
    copy!(ks.k3_Dx, dy.Dx)
end

function rkab_step!(h, cs, D::DerivativeOperator, ks::Substeps, yp::GridFuncs, y::GridFuncs, dy::GridFuncs)
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

function euler_step!(h, D::DerivativeOperator, y::GridFuncs, dy::GridFuncs)
    # Step 1: Compute RHS
    compute_rhs!(D, y, dy)

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