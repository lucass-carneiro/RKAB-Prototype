import domain as dom

import numpy as np


class StateVector:
    def __init__(self, time_levels: int, domain: dom.Domain):
        shape = (time_levels, domain.n_all, domain.n_all, domain.n_all)
        self.u = np.zeros(shape)
        self.rho = np.zeros(shape)
