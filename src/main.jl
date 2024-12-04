using ArgParse
using Profile
using PProf

include("evolve_1d.jl")
include("evolve_2d.jl")
include("evolve_3d.jl")

function main(args)
    s = ArgParseSettings(
        description = "Runge-Kutta + Adams–Bashforth hybrid wave equation integrator",
        version = "Version 1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--cpu-profile"
            help = "Run the CPU profiler"
            nargs = 0 
        "--mem-profile"
            help = "Run the memory profiler"
            nargs = 0 
        "mode"
            help = "Mode of the evolution"
            required = true
            default = "3D"
        "param_file"
            help = "Parameter file for driving the evolution"
            required = true
    end

    parsed_args = parse_args(args, s)

    if parsed_args["cpu-profile"]
        @info "Profiling CPU execution"
        
        # Collect a profile
        Profile.clear() 

        # Export pprof profile and open interactive profiling web interface.
        if parsed_args["mode"] == "1D" || parsed_args["mode"] == "1d"
            @profile evolve_1d(parsed_args["param_file"])
        elseif parsed_args["mode"] == "2D" || parsed_args["mode"] == "2d"
            @profile evolve_2d(parsed_args["param_file"])
        elseif parsed_args["mode"] == "3D" || parsed_args["mode"] == "3d"
            @profile evolve_3d(parsed_args["param_file"])
        else
            @error "Unrecognized mode of operation $(parsed_args["mode"])"
            return 1
        end

        pprof()
    elseif parsed_args["mem-profile"]
        @info "Profiling memory usage"
    elseif parsed_args["cpu-profile"] && parsed_args["mem-profile"]
        @info "Profiling CPU execution and memory usage"
    else
        @info "Regular evolution"
        if parsed_args["mode"] == "1D" || parsed_args["mode"] == "1d"
            evolve_1d(parsed_args["param_file"])
        elseif parsed_args["mode"] == "2D" || parsed_args["mode"] == "2d"
            evolve_2d(parsed_args["param_file"])
        elseif parsed_args["mode"] == "3D" || parsed_args["mode"] == "3d"
            evolve_3d(parsed_args["param_file"])
        else
            @error "Unrecognized mode of operation $(parsed_args["mode"])"
            return 1
        end
    end

    return 0
end

return main(ARGS)