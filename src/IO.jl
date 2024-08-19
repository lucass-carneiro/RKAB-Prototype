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

# The functions below use the "permutedims" construct to save row-major ordered arrays.
# This is required because Numpy (which is used for reading and plotting the data)
# uses row major ordering.

function write_state(group, i, y::GridFuncs3D)
    write(group, "Phi_$(lpad(i, 4, "0"))", permutedims(y.Phi, reverse(1:ndims(y.Phi))))
    write(group, "Pi_$(lpad(i, 4, "0"))", permutedims(y.Pi, reverse(1:ndims(y.Pi))))
    write(group, "Dx_$(lpad(i, 4, "0"))", permutedims(y.Dx, reverse(1:ndims(y.Dx))))
    write(group, "Dy_$(lpad(i, 4, "0"))", permutedims(y.Dy, reverse(1:ndims(y.Dy))))
    write(group, "Dz_$(lpad(i, 4, "0"))", permutedims(y.Dz, reverse(1:ndims(y.Dz))))
end

function write_rhs(group, i, dy::GridFuncs3D)
    write(group, "Phi_rhs_$(lpad(i, 4, "0"))", permutedims(dy.Phi, reverse(1:ndims(dy.Phi))))
    write(group, "Pi_rhs_$(lpad(i, 4, "0"))", permutedims(dy.Pi, reverse(1:ndims(dy.Pi))))
    write(group, "Dx_rhs_$(lpad(i, 4, "0"))", permutedims(dy.Dx, reverse(1:ndims(dy.Dx))))
    write(group, "Dy_rhs_$(lpad(i, 4, "0"))", permutedims(dy.Dy, reverse(1:ndims(dy.Dy))))
    write(group, "Dz_rhs_$(lpad(i, 4, "0"))", permutedims(dy.Dz, reverse(1:ndims(dy.Dz))))
end