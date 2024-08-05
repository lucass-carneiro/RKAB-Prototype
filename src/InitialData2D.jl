
function sw_Phi(A, kx, ky, t, x, y)
    omega = sqrt(kx * kx + ky * ky)
    return A*cos(2*omega*pi*t)*sin(2*kx*pi*x)*sin(2*ky*pi*y)
end

function sw_Pi(A, kx, ky, t, x, y)
    omega = sqrt(kx * kx + ky * ky)
    return -2*A*omega*pi*sin(2*omega*pi*t)*sin(2*kx*pi*x)*sin(2*ky*pi*y)
end

function sw_Dx(A, kx, ky, t, x, y)
    omega = sqrt(kx * kx + ky * ky)
    return 2*A*kx*pi*cos(2*omega*pi*t)*cos(2*kx*pi*x)*sin(2*ky*pi*y)
end

function sw_Dy(A, kx, ky, t, x, y)
    omega = sqrt(kx * kx + ky * ky)
    return 2*A*ky*pi*cos(2*omega*pi*t)*sin(2*kx*pi*x)*cos(2*ky*pi*y)
end