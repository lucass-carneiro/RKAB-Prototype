import parameters as par
import state_vector as stv

import numpy as np


def standing_wave_initialize(parameters: par.Parameters, state_vector: stv.StateVector):
    A = parameters.standing_wave_A

    kx = parameters.standing_wave_Kx
    ky = parameters.standing_wave_Ky
    kz = parameters.standing_wave_Kz

    omega = np.sqrt(kx * kx + ky * ky + kz * kz)

    for p in parameters.domain.all_points():
        t = 0
        i, j, k = p.idx()
        x, y, z = p.pos()

        u = A * np.cos(2 * np.pi * omega * t) * np.cos(2 * np.pi * kx * x) * \
            np.cos(2 * np.pi * ky * y) * np.cos(2 * np.pi * kz * z)

        rho = A * (-2 * np.pi * omega) * np.sin(2 * np.pi * omega * t) * \
            np.cos(2 * np.pi * kx * x) * np.cos(2 * np.pi * ky * y) * \
            np.cos(2 * np.pi * kz * z)

        state_vector.u[0, i, j, k] = u
        state_vector.rho[0, i, j, k] = rho


def initialize(parameters: par.Parameters, state_vector: stv.StateVector):
    print("Initializing")

    if (parameters.id_type == "standing-wave"):
        standing_wave_initialize(parameters, state_vector)
