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