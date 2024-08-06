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

function compute_Pi_rhs!(D::DerivativeOperator, d::Derivatives2D, Dx::Vector{Float64}, Dy::Vector{Float64}, Pi_rhs::Vector{Float64})
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

function compute_rhs!(D::DerivativeOperator, d::Derivatives2D, Pi::Vector{Float64}, Dx::Vector{Float64}, Dy::Vector{Float64}, dy::GridFuncs2D)
    compute_Phi_rhs!(Pi, dy.Phi)
    compute_Pi_rhs!(D, d, Dx, Dy, dy.Pi)
    du_dx!(D, Pi, dy.Dx)
    du_dy!(D, Pi, dy.Dy)
end