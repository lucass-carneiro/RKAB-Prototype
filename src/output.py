import state_vector as stv
import domain as dom

import adios2
import numpy as np


def write_domain_data(stream: adios2.Stream, domain: dom.Domain):
    stream.write_attribute("all_start", domain.all_start)
    stream.write_attribute("all_end", domain.all_end)

    stream.write_attribute("int_start", domain.int_start)
    stream.write_attribute("int_end", domain.int_end)

    stream.write_attribute("n_all", domain.n_all)
    stream.write_attribute("n_int", domain.n_int)

    stream.write_attribute("ghosts", domain.ghosts)

    stream.write_attribute("delta", domain.delta)


def write_time_step(stream: adios2.Stream, iteration: int, time: float, state_vector: stv.StateVector):
    print("Saving")

    stream.begin_step()

    stream.write("iteration", iteration)
    stream.write("time", time)

    shape = np.shape(state_vector.u[0])
    start = [0] * len(shape)
    count = shape

    stream.write("u", state_vector.u[0], shape, start, count)
    stream.write("rho", state_vector.rho[0], shape, start, count)

    stream.end_step()
