function sw_Phi(A, kx, t, x)
    omega = sqrt(kx * kx)
    return A*cos(2*omega*pi*t)*sin(2*kx*pi*x)
end

function sw_Pi(A, kx, t, x)
    omega = sqrt(kx * kx)
    return -2*A*omega*pi*sin(2*omega*pi*t)*sin(2*kx*pi*x)
end

function sw_Dx(A, kx, t, x)
    omega = sqrt(kx * kx)
    return 2*A*kx*pi*cos(2*omega*pi*t)*cos(2*kx*pi*x)
end