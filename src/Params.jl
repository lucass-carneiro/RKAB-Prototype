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

struct Params2D
    domain_start::Float64
    domain_end::Float64
    domain_points::Int

    time_final::Float64
    time_cfl::Float64

    standing_wave_A::Float64
    standing_wave_kx::Float64
    standing_wave_ky::Float64

    SBP_accuracy_order::Int

    RKAB_coeffs::Array{Float64}

    output_file::String

    function Params2D(config_file)
        config_data = YAML.load_file(config_file)

        domain_start = config_data["2d_params"]["domain"]["start"]
        domain_end = config_data["2d_params"]["domain"]["end"]
        domain_points = config_data["2d_params"]["domain"]["points"]

        time_final = config_data["2d_params"]["time"]["final"]
        time_cfl = config_data["2d_params"]["time"]["cfl"]

        standing_wave_A = config_data["2d_params"]["standing-wave"]["A"]
        standing_wave_kx = config_data["2d_params"]["standing-wave"]["kx"]
        standing_wave_ky = config_data["2d_params"]["standing-wave"]["ky"]

        SBP_accuracy_order = config_data["2d_params"]["SBP-accuracy-order"]

        c0 = 0.0
        c1 = config_data["2d_params"]["RKAB-coeffs"]["c1"]
        c2 = config_data["2d_params"]["RKAB-coeffs"]["c2"]
        c3 = config_data["2d_params"]["RKAB-coeffs"]["c3"]
        c4 = config_data["2d_params"]["RKAB-coeffs"]["c4"]
        c5 = config_data["2d_params"]["RKAB-coeffs"]["c5"]
        c6 = config_data["2d_params"]["RKAB-coeffs"]["c6"]
        c7 = config_data["2d_params"]["RKAB-coeffs"]["c7"]
        c8 = config_data["2d_params"]["RKAB-coeffs"]["c8"]
        c9 = 1 - c6 - c7 - c8

        RKAB_coeffs = [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9]

        output_file = config_data["2d_params"]["output"]

        return new(
            domain_start,
            domain_end,
            domain_points,
            time_final,
            time_cfl,
            standing_wave_A,
            standing_wave_kx,
            standing_wave_ky,
            SBP_accuracy_order,
            RKAB_coeffs,
            output_file
        )
    end
end

struct Params3D
    domain_start::Float64
    domain_end::Float64
    domain_points::Int

    time_final::Float64
    time_cfl::Float64

    standing_wave_A::Float64
    standing_wave_kx::Float64
    standing_wave_ky::Float64
    standing_wave_kz::Float64

    SBP_accuracy_order::Int

    RKAB_coeffs::Array{Float64}

    output_file::String

    function Params3D(config_file)
        config_data = YAML.load_file(config_file)

        domain_start = config_data["3d_params"]["domain"]["start"]
        domain_end = config_data["3d_params"]["domain"]["end"]
        domain_points = config_data["3d_params"]["domain"]["points"]

        time_final = config_data["3d_params"]["time"]["final"]
        time_cfl = config_data["3d_params"]["time"]["cfl"]

        standing_wave_A = config_data["3d_params"]["standing-wave"]["A"]
        standing_wave_kx = config_data["3d_params"]["standing-wave"]["kx"]
        standing_wave_ky = config_data["3d_params"]["standing-wave"]["ky"]
        standing_wave_kz = config_data["3d_params"]["standing-wave"]["kz"]

        SBP_accuracy_order = config_data["3d_params"]["SBP-accuracy-order"]

        c0 = 0.0
        c1 = config_data["3d_params"]["RKAB-coeffs"]["c1"]
        c2 = config_data["3d_params"]["RKAB-coeffs"]["c2"]
        c3 = config_data["3d_params"]["RKAB-coeffs"]["c3"]
        c4 = config_data["3d_params"]["RKAB-coeffs"]["c4"]
        c5 = config_data["3d_params"]["RKAB-coeffs"]["c5"]
        c6 = config_data["3d_params"]["RKAB-coeffs"]["c6"]
        c7 = config_data["3d_params"]["RKAB-coeffs"]["c7"]
        c8 = config_data["3d_params"]["RKAB-coeffs"]["c8"]
        c9 = 1 - c6 - c7 - c8

        RKAB_coeffs = [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9]

        output_file = config_data["3d_params"]["output"]

        return new(
            domain_start,
            domain_end,
            domain_points,
            time_final,
            time_cfl,
            standing_wave_A,
            standing_wave_kx,
            standing_wave_ky,
            standing_wave_kz,
            SBP_accuracy_order,
            RKAB_coeffs,
            output_file
        )
    end
end