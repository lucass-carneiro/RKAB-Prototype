using HDF5

function write_state(group, i, y::GridFuncs)
    write(group, "Phi_$(lpad(i, 8, "0"))", y.Phi)
    write(group, "Pi_$(lpad(i, 8, "0"))", y.Pi)
    write(group, "Dx_$(lpad(i, 8, "0"))", y.Dx)
end

function write_rhs(group, i, dy::GridFuncs)
    write(group, "Phi_rhs_$(lpad(i, 8, "0"))", dy.Phi)
    write(group, "Pi_rhs_$(lpad(i, 8, "0"))", dy.Pi)
    write(group, "Dx_rhs_$(lpad(i, 8, "0"))", dy.Dx)
end

function write_state(group, i, y::GridFuncs2D)
    write(group, "Phi_$(lpad(i, 8, "0"))", y.Phi)
    write(group, "Pi_$(lpad(i, 8, "0"))", y.Pi)
    write(group, "Dx_$(lpad(i, 8, "0"))", y.Dx)
    write(group, "Dy_$(lpad(i, 8, "0"))", y.Dy)
end

function write_rhs(group, i, dy::GridFuncs2D)
    write(group, "Phi_rhs_$(lpad(i, 8, "0"))", dy.Phi)
    write(group, "Pi_rhs_$(lpad(i, 8, "0"))", dy.Pi)
    write(group, "Dx_rhs_$(lpad(i, 8, "0"))", dy.Dx)
    write(group, "Dy_rhs_$(lpad(i, 8, "0"))", dy.Dy)
end

function write_state(group, i, y::GridFuncs3D)
    write(group, "Phi_$(lpad(i, 8, "0"))", y.Phi)
    write(group, "Pi_$(lpad(i, 8, "0"))", y.Pi)
    write(group, "Dx_$(lpad(i, 8, "0"))", y.Dx)
    write(group, "Dy_$(lpad(i, 8, "0"))", y.Dy)
    write(group, "Dz_$(lpad(i, 8, "0"))", y.Dz)
end

function write_rhs(group, i, dy::GridFuncs3D)
    write(group, "Phi_rhs_$(lpad(i, 8, "0"))", dy.Phi)
    write(group, "Pi_rhs_$(lpad(i, 8, "0"))", dy.Pi)
    write(group, "Dx_rhs_$(lpad(i, 8, "0"))", dy.Dx)
    write(group, "Dy_rhs_$(lpad(i, 8, "0"))", dy.Dy)
    write(group, "Dz_rhs_$(lpad(i, 8, "0"))", dy.Dz)
end