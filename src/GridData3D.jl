struct GridFuncs3D
    Phi::Array{Float64}
    Pi::Array{Float64}
    Dx::Array{Float64}
    Dy::Array{Float64}
    Dz::Array{Float64}

    function GridFuncs3D(n)
        num_points = (n, n, n)
        new(
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points)
        )
    end
end

struct Substeps3D
    k0_Phi::Array{Float64}
    k0_Pi::Array{Float64}
    k0_Dx::Array{Float64}
    k0_Dy::Array{Float64}
    k0_Dz::Array{Float64}
    
    k1_Phi::Array{Float64}
    k1_Pi::Array{Float64}
    k1_Dx::Array{Float64}
    k1_Dy::Array{Float64}
    k1_Dz::Array{Float64}
    
    k2_Phi::Array{Float64}
    k2_Pi::Array{Float64}
    k2_Dx::Array{Float64}
    k2_Dy::Array{Float64}
    k2_Dz::Array{Float64}
    
    k3_Phi::Array{Float64}
    k3_Pi::Array{Float64}
    k3_Dx::Array{Float64}
    k3_Dy::Array{Float64}
    k3_Dz::Array{Float64}

    function Substeps3D(n)
        num_points = (n, n, n)
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

struct Derivatives3D
    dDx_dx::Array{Float64}
    dDy_dy::Array{Float64}
    dDz_dz::Array{Float64}

    function Derivatives3D(n)
        num_points = (n, n, n)
        new(
            zeros(Float64, num_points),
            zeros(Float64, num_points),
            zeros(Float64, num_points)
        )
    end
end