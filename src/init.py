import parameters as par
import state_vector as stv

import numpy as np


def standing_wave_initialize(parameters: par.Parameters, state_vector: stv.StateVector):
    A = parameters.standing_wave_A

    kx = parameters.standing_wave_Kx
    ky = parameters.standing_wave_Ky
    kz = parameters.standing_wave_Kz

    omega = np.sqrt(kx * kx + ky * ky + kz * kz)

    def u(t, x, y, z):
        return A * np.cos(2 * np.pi * omega * t) * np.cos(2 * np.pi * kx * x) * \
            np.cos(2 * np.pi * ky * y) * np.cos(2 * np.pi * kz * z)

    def rho(t, x, y, z):
        return A * (-2 * np.pi * omega) * np.sin(2 * np.pi * omega * t) * \
            np.cos(2 * np.pi * kx * x) * np.cos(2 * np.pi * ky * y) * \
            np.cos(2 * np.pi * kz * z)

    for p in parameters.domain.all_points():
        i, j, k = p.idx()
        x, y, z = p.pos()

        # Current time level
        state_vector.u[0, i, j, k] = u(0, x, y, z)
        state_vector.rho[0, i, j, k] = rho(0, x, y, z)

        # Previous time level
        state_vector.u[1, i, j, k] = u(-parameters.domain.dt, x, y, z)
        state_vector.rho[1, i, j, k] = rho(-parameters.domain.dt, x, y, z)


def initialize(parameters: par.Parameters, state_vector: stv.StateVector):
    print("Initializing")

    if (parameters.id_type == "standing-wave"):
        standing_wave_initialize(parameters, state_vector)
