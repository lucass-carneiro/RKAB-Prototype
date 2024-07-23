import domain as dom
import state_vector as stv

import numpy as np


class RightHandSide:
    def __init__(self, domain: dom.Domain):
        shape = (domain.n_all, domain.n_all, domain.n_all)
        self.u_rhs = np.zeros(shape)
        self.rho_rhs = np.zeros(shape)
