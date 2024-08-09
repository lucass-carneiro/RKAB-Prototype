using HDF5
using YAML

struct Params1D
    domain_start::Float64
    domain_end::Float64
    domain_points::Int

    time_final::Float64
    time_cfl::Float64

    standing_wave_A::Float64
    standing_wave_kx::Float64

    SBP_accuracy_order::Int

    RKAB_coeffs::Array{Float64}

    output_file::String

    function Params1D(config_file)
        config_data = YAML.load_file(config_file)

        domain_start = config_data["1d_params"]["domain"]["start"]
        domain_end = config_data["1d_params"]["domain"]["end"]
        domain_points = config_data["1d_params"]["domain"]["points"]

        time_final = config_data["1d_params"]["time"]["final"]
        time_cfl = config_data["1d_params"]["time"]["cfl"]

        standing_wave_A = config_data["1d_params"]["standing-wave"]["A"]
        standing_wave_kx = config_data["1d_params"]["standing-wave"]["kx"]

        SBP_accuracy_order = config_data["1d_params"]["SBP-accuracy-order"]

        c0 = 0.0
        c1 = config_data["1d_params"]["RKAB-coeffs"]["c1"]
        c2 = config_data["1d_params"]["RKAB-coeffs"]["c2"]
        c3 = config_data["1d_params"]["RKAB-coeffs"]["c3"]
        c4 = config_data["1d_params"]["RKAB-coeffs"]["c4"]
        c5 = config_data["1d_params"]["RKAB-coeffs"]["c5"]
        c6 = config_data["1d_params"]["RKAB-coeffs"]["c6"]
        c7 = config_data["1d_params"]["RKAB-coeffs"]["c7"]
        c8 = config_data["1d_params"]["RKAB-coeffs"]["c8"]
        c9 = 1 - c6 - c7 - c8

        RKAB_coeffs = [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9]

        output_file = config_data["1d_params"]["output"]

        return new(
            domain_start,
            domain_end,
            domain_points,
            time_final,
            time_cfl,
            standing_wave_A,
            standing_wave_kx,
            SBP_accuracy_order,
            RKAB_coeffs,
            output_file
        )
    end
end

function write_state(group, i, y::GridFuncs)
    write(group, "Phi_$(lpad(i, 4, "0"))", y.Phi)
    write(group, "Pi_$(lpad(i, 4, "0"))", y.Pi)
    write(group, "Dx_$(lpad(i, 4, "0"))", y.Dx)
end

function write_rhs(group, i, dy::GridFuncs)
    write(group, "Phi_rhs_$(lpad(i, 4, "0"))", dy.Phi)
    write(group, "Pi_rhs_$(lpad(i, 4, "0"))", dy.Pi)
    write(group, "Dx_rhs_$(lpad(i, 4, "0"))", dy.Dx)
end

function write_state(group, i, y::GridFuncs2D)
    write(group, "Phi_$(lpad(i, 4, "0"))", y.Phi)
    write(group, "Pi_$(lpad(i, 4, "0"))", y.Pi)
    write(group, "Dx_$(lpad(i, 4, "0"))", y.Dx)
    write(group, "Dy_$(lpad(i, 4, "0"))", y.Dy)
end

function write_rhs(group, i, dy::GridFuncs2D)
    write(group, "Phi_rhs_$(lpad(i, 4, "0"))", dy.Phi)
    write(group, "Pi_rhs_$(lpad(i, 4, "0"))", dy.Pi)
    write(group, "Dx_rhs_$(lpad(i, 4, "0"))", dy.Dx)
    write(group, "Dy_rhs_$(lpad(i, 4, "0"))", dy.Dy)
end