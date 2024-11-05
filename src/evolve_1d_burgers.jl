using RKAB

using HDF5
using LinearAlgebra
using SummationByPartsOperators
using Random

function evolve_1d_hydro(config_file)
    @info "Using parameter file $config_file"

    config_data = Params1DBurgers(config_file)

    # Domain
    r0 = config_data.domain_start
    rf = config_data.domain_end
    num_cells = config_data.domain_cells
    dr = (rf - r0) / num_cells

    # Time steps
    final_time = config_data.time_final
    cfl = config_data.time_cfl
    dt = cfl * dr
    last_iter = floor(Int, final_time / dt) + 1
    save_every = config_data.output_every

    @info "Spatial step = $dr"
    @info "Time step = $dt"
    @info "Total iterations: $last_iter"
    @info "Save every $(save_every) iterations"

    # IO
    h5_file = h5open(config_data.output_file, "w")
    #state_group = create_group(h5_file, "state")
    #rhs_group = create_group(h5_file, "rhs")
    #grid_group = create_group(h5_file, "grid")

    close(h5_file)
end