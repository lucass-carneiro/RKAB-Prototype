import state_vector as stv
import right_hand_side as rhs
import domain as dom


def compute_rhs(domain: dom.Domain, state: stv.StateVector, rhs: rhs.RightHandSide):
    h = 1.0 / (domain.delta * domain.delta)

    for p in domain.interior_points():
        i, j, k = p.idx()

        d2udx2 = (state.u[0, i - 1, j, k] - 2.0 *
                  state.u[0, i, j, k] + state.u[0, i + 1, j, k]) * h

        d2udy2 = (state.u[0, i, j - 1, k] - 2.0 *
                  state.u[0, i, j, k] + state.u[0, i, j + 1, k]) * h

        d2udz2 = (state.u[0, i, j, k - 1] - 2.0 *
                  state.u[0, i, j, k] + state.u[0, i, j, k + 1]) * h

        rhs.u_rhs[i, j, k] = state.rho[0, i, j, k]
        rhs.rho_rhs[i, j, k] = d2udx2 + d2udy2 + d2udz2


def update(domain: dom.Domain, state: stv.StateVector, rhs: rhs.RightHandSide):
    for p in domain.interior_points():
        i, j, k = p.idx()

        # 1. Copy the current time level (0) to the previous time level (1)
        state.u[1, i, j, k] = state.u[0, i, j, k]
        state.rho[1, i, j, k] = state.rho[0, i, j, k]

        # 2. Update the current time level
        state.u[0, i, j, k] = state.u[0, i, j, k] + \
            domain.dt * rhs.u_rhs[i, j, k]
        state.rho[0, i, j, k] = state.rho[0, i, j, k] + \
            domain.dt * rhs.rho_rhs[i, j, k]


def apply_bcs(domain: dom.Domain, state: stv.StateVector):
    for p in domain.ghost_points():
        i, j, k = p.idx()

        state.u[0, i, j, k] = 0.0
        state.rho[0, i, j, k] = 0.0
