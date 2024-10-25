"""RKAB Plotter

Usage:
  plot.py 1d [options] <data-file>
  plot.py 2d [options] <data-file>
  plot.py 3d [options] <data-file>
  plot.py conv-time-1d [options] <coarse-data> <medium-data> <fine-data>
  plot.py (-h | --help)
  plot.py --version

Options:
  -h --help               Show this screen.
  --version               Show version.
  --iterations=<number>   Plot a specific iteration or all [default: all]
  --font-size=<size>      The size of the font in plots [default: 18]
  --slice-value=<value>   The value of the slicing coordinate [default: 0.0]
  --autorange             Automatically select plot ranges. If set, ignores varmin and varmax.
  --varmin=<value>        Min. value for the plot variable [default: -1.0].
  --varmax=<value>        Max. value for the plot variable [default: 1.0].
"""
from docopt import docopt

from plot3d import plot_3d

import h5py
import matplotlib as mpl
import sys

import logging
logger = logging.getLogger(__name__)


def main(args):
    logging.basicConfig(
        format="[%(asctime)s] [PID: %(process)d] %(levelname)s: %(message)s",
        # filename="mpx.log",
        # filemode="w",
        stream=sys.stdout,
        level=logging.INFO,
    )

    font_size = int(args["--font-size"])

    mpl.use("agg")
    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["font.family"] = "Latin Modern Roman"
    mpl.rcParams["xtick.labelsize"] = font_size
    mpl.rcParams["ytick.labelsize"] = font_size

    if args["1d"]:
        data_file_name = args["<data-file>"]
        h5_file = h5py.File(data_file_name, "r")

        # plot_1d(args, font_size, h5_file)

        h5_file.close()
    elif args["2d"]:
        data_file_name = args["<data-file>"]
        h5_file = h5py.File(data_file_name, "r")

        # plot_2d(args, font_size, h5_file)

        h5_file.close()
    elif args["3d"]:
        data_file_name = args["<data-file>"]
        h5_file = h5py.File(data_file_name, "r")

        plot_3d(args, font_size, h5_file)

        h5_file.close()
    elif args["conv-time-1d"]:
        coarse_data = args["<coarse-data>"]
        medium_data = args["<medium-data>"]
        fine_data = args["<fine-data>"]

        h5_coarse = h5py.File(coarse_data, "r")
        h5_medium = h5py.File(medium_data, "r")
        h5_fine = h5py.File(fine_data, "r")

        # conv_time_1d(args, font_size, h5_coarse, h5_medium, h5_fine)

        h5_coarse.close()
        h5_medium.close()
        h5_fine.close()


if __name__ == '__main__':
    arguments = docopt(__doc__, version="RKAB Plot 1.0")
    main(arguments)
