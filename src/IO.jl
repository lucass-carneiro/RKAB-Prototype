using HDF5

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