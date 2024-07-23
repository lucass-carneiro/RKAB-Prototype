"""RKAB-Prototype.

Usage:
  main.py evolve <parameter-file>
  main.py plot <parameter-file>
  main.py (-h | --help)
  main.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
"""
from docopt import docopt

import parameters as par
import state_vector as stv
import right_hand_side as rhs
import output as out
import plots as plt
import update as upd
import init

import adios2


def evolve(args):
    # Parse parameters
    parameters = par.Parameters(args["<parameter-file>"])

    # Allocate grid functions
    state_vector = stv.StateVector(2, parameters.domain)
    right_hand_side = rhs.RightHandSide(parameters.domain)

    # Create output stream
    output_stream = adios2.Stream("output.bp", "w")
    out.write_domain_data(output_stream, parameters.domain)

    # Initialize the state vector
    init.initialize(parameters, state_vector)
    out.write_time_step(output_stream, 0, 0, state_vector)

    # Time stepping loop
    dt = parameters.domain.delta * parameters.courant_factor

    for iteration in range(1, parameters.last_iter + 1):
        # Current time
        t = iteration * dt

        print(f"Computing iteration {iteration}, t = {t}")

        # Compute RHS
        upd.compute_rhs(parameters.domain, state_vector, right_hand_side)

        # Time step
        upd.update(parameters.domain, dt, state_vector, right_hand_side)

        # Apply BCs
        upd.apply_bcs(parameters.domain, state_vector)

        # Output
        out.write_time_step(output_stream, iteration, t, state_vector)

    # Close output stream
    output_stream.close()


def plot(args):
    parameters = par.Parameters(args["<parameter-file>"])

    input_stream = adios2.Stream("output.bp", "r")

    plt.plot_time_steps(parameters.domain, input_stream)

    input_stream.close()


def main(args):
    if (args["evolve"]):
        evolve(args)
    elif (args["plot"]):
        plot(args)


if __name__ == "__main__":
    args = docopt(__doc__, version="RKAB-Prototype 1.0")
    main(args)
