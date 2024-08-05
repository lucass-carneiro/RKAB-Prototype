"""RKAB Plotter

Usage:
  plot.py 1d [options] <data-file>
  plot.py (-h | --help)
  plot.py --version

Options:
  -h --help               Show this screen.
  --version               Show version.
  --iterations=<number>   Plot a specific iteration or all [default: all]
  --font-size=<size>      The size of the font in plots [default: 18]
"""
from docopt import docopt

import numpy as np

import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt

import os
import sys

import logging
logger = logging.getLogger(__name__)


def Phi_1d(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return A*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)


def Pi_1d(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return -2*A*omega*np.pi*np.sin(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)


def Dx_1d(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return 2*A*kx*np.pi*np.cos(2*omega*np.pi*t)*np.cos(2*kx*np.pi*x)


def plot_gfs_1d(h5_file, prefix, gfs, iteration, iteration_string, dt, x, font_size, path):
    t = iteration * dt

    plt.close("all")

    for (name, color) in gfs:
        gf_path = f"{prefix}/{name}_{iteration_string}"
        plt.plot(x, h5_file[gf_path][:], color=color, label=name)

    plt.xlim(-1.0, 1.0)

    plt.xlabel("$x$", size=font_size)
    plt.ylabel("$f(x)$", size=font_size)

    plt.title(f"{prefix} at iteration {iteration}, $t = {t}$")
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            path,
            f"1D_{prefix}_it_{iteration_string}.pdf"
        )
    )


def plot_expected_1d(h5_file, prefix, gfs, iteration, iteration_string, dt, x, font_size, path):
    A = h5_file.attrs["A"]
    kx = h5_file.attrs["kx"]
    t = iteration * dt

    plt.close("all")

    for (name, color) in gfs:
        gf_path = f"{prefix}/{name}_{iteration_string}"
        expected = eval(f"{name}_1d")(A, kx, t, x)
        data = np.abs(h5_file[gf_path][:] - expected)
        plt.plot(x, data, color=color, label=f"|True{name} - {name}|")

    plt.xlim(-1.0, 1.0)

    plt.xlabel("$x$", size=font_size)
    plt.ylabel("$f(x)$", size=font_size)

    plt.title(f"{prefix} at iteration {iteration}, $t = {t}$")
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            path,
            f"1D_expected_{prefix}_it_{iteration_string}.pdf"
        )
    )


def plot_1d(args, font_size, gfs, rhs_gfs, h5_file):
    last_iter = h5_file.attrs["last_iter"]
    dt = h5_file.attrs["dt"]
    x = h5_file["grid/x_coords"][:]

    if args["--iterations"] != "all":
        iteration_range = range(
            int(args["--iterations"]),
            int(args["--iterations"]) + 1
        )
    else:
        iteration_range = range(
            0,
            last_iter + 1
        )

    plots_dir = "1d_plots"
    state_dir = os.path.join(plots_dir, "state")
    rhs_dir = os.path.join(plots_dir, "rhs")
    expected_dir = os.path.join(plots_dir, "expected")

    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    if not os.path.exists(state_dir):
        os.mkdir(state_dir)

    if not os.path.exists(rhs_dir):
        os.mkdir(rhs_dir)

    if not os.path.exists(expected_dir):
        os.mkdir(expected_dir)

    for iteration in iteration_range:
        iteration_string = f"{iteration:04}"

        logging.info(f"Plotting iteration {iteration_string}")

        plot_gfs_1d(
            h5_file,
            "state",
            gfs,
            iteration,
            iteration_string,
            dt,
            x,
            font_size,
            state_dir
        )

        plot_gfs_1d(
            h5_file,
            "rhs",
            rhs_gfs,
            iteration,
            iteration_string,
            dt,
            x,
            font_size,
            rhs_dir
        )

        plot_expected_1d(
            h5_file,
            "state",
            gfs,
            iteration,
            iteration_string,
            dt,
            x,
            font_size,
            expected_dir
        )


def main(args):
    logging.basicConfig(
        format="[%(asctime)s] [PID: %(process)d] %(levelname)s: %(message)s",
        # filename="mpx.log",
        # filemode="w",
        stream=sys.stdout,
        level=logging.INFO,
    )

    font_size = int(args["--font-size"])
    data_file_name = args["<data-file>"]

    h5_file = h5py.File(data_file_name, "r")

    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["font.family"] = "Latin Modern Roman"
    mpl.rcParams["xtick.labelsize"] = font_size
    mpl.rcParams["ytick.labelsize"] = font_size

    gfs = [
        ("Phi", "red"),
        ("Pi", "black"),
        ("Dx", "blue")
    ]

    rhs_gfs = [
        ("Phi_rhs", "red"),
        ("Pi_rhs", "black"),
        ("Dx_rhs", "blue")
    ]

    if args["1d"]:
        plot_1d(args, font_size, gfs, rhs_gfs, h5_file)

    h5_file.close()


if __name__ == '__main__':
    arguments = docopt(__doc__, version="RKAB Plot 1.0")
    main(arguments)
