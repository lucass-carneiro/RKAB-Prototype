
function sw_Phi(A, kx, ky, kz, t, x, y, z)
    omega = sqrt(kx * kx + ky * ky + kz * kz)
    return A*cos(2*omega*pi*t)*sin(2*kx*pi*x)*sin(2*ky*pi*y)*sin(2*kz*pi*z)
end

function sw_Pi(A, kx, ky, kz, t, x, y, z)
    omega = sqrt(kx * kx + ky * ky + kz * kz)
    return -2*A*omega*pi*sin(2*omega*pi*t)*sin(2*kx*pi*x)*sin(2*ky*pi*y)*sin(2*kz*pi*z)
end

function sw_Dx(A, kx, ky, kz, t, x, y, z)
    omega = sqrt(kx * kx + ky * ky + kz * kz)
    return 2*A*kx*pi*cos(2*omega*pi*t)*cos(2*kx*pi*x)*sin(2*ky*pi*y)*sin(2*kz*pi*z)
end

function sw_Dy(A, kx, ky, kz, t, x, y, z)
    omega = sqrt(kx * kx + ky * ky + kz * kz)
    return 2*A*ky*pi*cos(2*omega*pi*t)*sin(2*kx*pi*x)*cos(2*ky*pi*y)*sin(2*kz*pi*z)
end

function sw_Dz(A, kx, ky, kz, t, x, y, z)
    omega = sqrt(kx * kx + ky * ky + kz * kz)
    return 2*A*kz*pi*cos(2*omega*pi*t)*sin(2*kx*pi*x)*sin(2*ky*pi*y)*cos(2*kz*pi*z)
end