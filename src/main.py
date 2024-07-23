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
