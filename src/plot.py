"""RKAB Plotter

Usage:
  plot.py 1d [options] <data-file>
  plot.py 2d [options] <data-file>
  plot.py conv-time-1d [options] <coarse-data> <medium-data> <fine-data>
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

import concurrent.futures

import logging
logger = logging.getLogger(__name__)


def sw_Phi_1d(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return A*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)


def sw_Pi_1d(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return -2*A*omega*np.pi*np.sin(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)


def sw_Dx_1d(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return 2*A*kx*np.pi*np.cos(2*omega*np.pi*t)*np.cos(2*kx*np.pi*x)


def sw_Phi_2d(A, kx, ky, t, x, y):
    omega = np.sqrt(kx * kx + ky * ky)
    return A*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)*np.sin(2*ky*np.pi*y)


def sw_Pi_2d(A, kx, ky, t, x, y):
    omega = np.sqrt(kx * kx + ky * ky)
    return -2*A*omega*np.pi*np.sin(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)*np.sin(2*ky*np.pi*y)


def sw_Dx_2d(A, kx, ky, t, x, y):
    omega = np.sqrt(kx * kx + ky * ky)
    return 2*A*kx*np.pi*np.cos(2*omega*np.pi*t)*np.cos(2*kx*np.pi*x)*np.sin(2*ky*np.pi*y)


def sw_Dy_2d(A, kx, ky, t, x, y):
    omega = np.sqrt(kx * kx + ky * ky)
    return 2*A*ky*np.pi*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)*np.cos(2*ky*np.pi*y)


def plot_gfs_1d(x, data, name, font_size, prefix, iteration, iteration_string, t, path):
    plt.close("all")

    plt.plot(x, data, color="black")

    plt.xlim(-1.0, 1.0)

    plt.xlabel("$x$", size=font_size)
    plt.ylabel("$f(x)$", size=font_size)

    plt.title(f"{prefix}/{name} at iteration {iteration}, $t = {t}$")

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            path,
            f"1D_{prefix}_{name}_it_{iteration_string}.pdf"
        )
    )


def plot_expected_1d(x, data, name, A, kx, font_size, prefix, iteration, iteration_string, t, path):
    plt.close("all")

    expected = eval(f"{name}_1d")(A, kx, t, x)
    y = np.abs(data - expected)
    plt.plot(x, y, color="black")

    plt.xlim(-1.0, 1.0)
    plt.ylim(0.0, 8.0e-3)

    plt.xlabel("$x$", size=font_size)
    plt.ylabel(f"|True{name} - {name}|", size=font_size)

    plt.title(f"{prefix}/{name} error at iteration {iteration}, $t = {t}$")

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            path,
            f"1D_expected_{prefix}_{name}_it_{iteration_string}.pdf"
        )
    )


def plot_1d(args, font_size, h5_file):
    last_iter = h5_file.attrs["last_iter"]
    dt = h5_file.attrs["dt"]

    x = h5_file["grid/x_coords"][:]

    A = h5_file.attrs["A"]
    kx = h5_file.attrs["kx"]

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

    gfs = [
        "Phi",
        "Pi",
        "Dx",
    ]

    rhs_gfs = [
        "Phi_rhs",
        "Pi_rhs",
        "Dx_rhs",
    ]

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

    level_array = np.linspace(-1.0, 1.0, endpoint=True, num=101)
    err_level_array = np.linspace(0.0, 0.003, endpoint=True, num=101)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in iteration_range:
            logger.info(f"Dispatching plot of iteration {i}")

            iteration_string = f"{i:04}"
            t = i * dt

            # Plot State
            for gf in gfs:
                path = f"state/{gf}_{iteration_string}"
                data = h5_file[path][:]

                executor.submit(
                    plot_gfs_1d,
                    x,
                    data,
                    gf,
                    font_size,
                    "state",
                    i,
                    iteration_string,
                    t,
                    state_dir
                )

            # Plot RHS
            for gf in rhs_gfs:
                path = f"rhs/{gf}_{iteration_string}"
                data = h5_file[path][:]

                executor.submit(
                    plot_gfs_1d,
                    x,
                    data,
                    gf,
                    font_size,
                    "rhs",
                    i,
                    iteration_string,
                    t,
                    rhs_dir
                )

            # Plot error
            for gf in gfs:
                path = f"state/{gf}_{iteration_string}"
                data = h5_file[path][:]

                executor.submit(
                    plot_expected_1d,
                    x,
                    data,
                    gf,
                    A,
                    kx,
                    font_size,
                    "state",
                    i,
                    iteration_string,
                    t,
                    expected_dir
                )


def plot_gfs_2d(x, y, data, levels, font_size, prefix, name, iteration, iteration_string, t, path):
    plt.close("all")
    plt.contourf(
        x,
        y,
        data,
        levels=levels,
        cmap="RdBu"
    )

    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)

    plt.xlabel("$x$", size=font_size)
    plt.ylabel("$y$", size=font_size)

    plt.title(f"{prefix}/{name} at iteration {iteration}, $t = {t}$")

    cb = plt.colorbar()
    cb.ax.set_ylabel(name)

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            path,
            f"2D_{prefix}_{name}_it_{iteration_string}.png"
        )
    )


def plot_expected_2d(x, y, z, levels, font_size, prefix, name, iteration, iteration_string, t, path):
    plt.close("all")

    # We transpose the data because Julia uses column major ordering
    plt.contourf(
        x,
        y,
        z,
        levels=levels
    )

    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.0, 1.0)

    plt.xlabel("$x$", size=font_size)
    plt.ylabel("$y$", size=font_size)

    plt.title(f"{prefix}/{name} error at iteration {iteration}, $t = {t}$")

    cb = plt.colorbar()
    cb.ax.set_ylabel(name)

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            path,
            f"2D_{prefix}_{name}_it_{iteration_string}.png"
        )
    )


def plot_2d(args, font_size, h5_file):
    last_iter = h5_file.attrs["last_iter"]
    dt = h5_file.attrs["dt"]

    x = h5_file["grid/x_coords"][:]
    y = h5_file["grid/y_coords"][:]
    X, Y = np.meshgrid(x, y)

    A = h5_file.attrs["A"]
    kx = h5_file.attrs["kx"]
    ky = h5_file.attrs["kx"]

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

    gfs = [
        "Phi",
        "Pi",
        "Dx",
        "Dy",
    ]

    rhs_gfs = [
        "Phi_rhs",
        "Pi_rhs",
        "Dx_rhs",
        "Dy_rhs",
    ]

    plots_dir = "2d_plots"
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

    level_array = np.linspace(-1.0, 1.0, endpoint=True, num=101)
    err_level_array = np.linspace(0.0, 0.003, endpoint=True, num=101)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in iteration_range:
            logger.info(f"Dispatching plot of iteration {i}")

            iteration_string = f"{i:04}"
            t = i * dt

            # Plot State
            for gf in gfs:
                path = f"state/{gf}_{iteration_string}"
                data = h5_file[path][:]

                if gf == "Phi":
                    levels = level_array
                else:
                    levels = 100

                executor.submit(
                    plot_gfs_2d,
                    x,
                    y,
                    data,
                    levels,
                    font_size,
                    "state",
                    gf,
                    i,
                    iteration_string,
                    t,
                    state_dir
                )

            # Plot RHS
            for gf in rhs_gfs:
                path = f"rhs/{gf}_{iteration_string}"
                data = h5_file[path][:]

                executor.submit(
                    plot_gfs_2d,
                    x,
                    y,
                    data,
                    100,
                    font_size,
                    "rhs",
                    gf,
                    i,
                    iteration_string,
                    t,
                    rhs_dir
                )

            # Plot Error
            for gf in gfs:
                path = f"state/{gf}_{iteration_string}"
                data = h5_file[path][:]

                z = np.abs(
                    eval(f"sw_{gf}_2d")(A, kx, ky, t, X, Y) - data
                )

                executor.submit(
                    plot_expected_2d,
                    x,
                    y,
                    z,
                    err_level_array,
                    font_size,
                    "state",
                    gf,
                    i,
                    iteration_string,
                    t,
                    expected_dir
                )

        logger.info(f"Waiting for plot workers to finish")


def conv_time_1d(args, font_size, h5_coarse, h5_medium, h5_fine):
    A = h5_coarse.attrs["A"]
    kx = h5_coarse.attrs["kx"]

    x = h5_coarse["grid/x_coords"][:]

    last_iter_coarse = h5_coarse.attrs["last_iter"]
    last_iter_medium = h5_medium.attrs["last_iter"]
    last_iter_fine = h5_fine.attrs["last_iter"]

    dt_coarse = h5_coarse.attrs["dt"]
    dt_medium = h5_medium.attrs["dt"]
    dt_fine = h5_fine.attrs["dt"]

    coarse_times = np.linspace(
        0.0,
        last_iter_coarse * dt_coarse,
        endpoint=True,
        num=(last_iter_coarse + 1)
    )

    medium_times = np.linspace(
        0.0,
        last_iter_medium * dt_medium,
        endpoint=True,
        num=(last_iter_medium + 1)
    )

    fine_times = np.linspace(
        0.0,
        last_iter_fine * dt_fine,
        endpoint=True,
        num=(last_iter_fine + 1)
    )

    medium_coarse_mask = []
    for i in range(last_iter_medium + 1):
        m_time = i * dt_medium

        for c_time in coarse_times:
            if np.isclose(m_time, c_time):
                medium_coarse_mask.append(i)

    fine_coarse_mask = []
    for i in range(last_iter_fine + 1):
        f_time = i * dt_fine

        for c_time in coarse_times:
            if np.isclose(f_time, c_time):
                fine_coarse_mask.append(i)

    coarse_test_iteration_idx = 1
    test_time = coarse_test_iteration_idx * dt_coarse

    medium_test_iteration_idx = medium_coarse_mask[coarse_test_iteration_idx]
    fine_test_iteration_idx = fine_coarse_mask[coarse_test_iteration_idx]

    coarse_test_iteration_idx_string = f"{coarse_test_iteration_idx:04}"
    medium_test_iteration_idx_string = f"{medium_test_iteration_idx:04}"
    fine_test_iteration_idx_string = f"{fine_test_iteration_idx:04}"

    coarse_gf_string = f"state/Phi_{coarse_test_iteration_idx_string}"
    medium_gf_string = f"state/Phi_{medium_test_iteration_idx_string}"
    fine_gf_string = f"state/Phi_{fine_test_iteration_idx_string}"

    coarse_gf = h5_coarse[coarse_gf_string][:]
    medium_gf = h5_medium[medium_gf_string][:]
    fine_gf = h5_fine[fine_gf_string][:]

    exact_gf = sw_Phi_1d(A, kx, test_time, x)

    print(np.linalg.norm(exact_gf - coarse_gf))
    print(np.linalg.norm(exact_gf - medium_gf))
    print(np.linalg.norm(exact_gf - fine_gf))


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

        plot_1d(args, font_size, h5_file)

        h5_file.close()
    elif args["2d"]:
        data_file_name = args["<data-file>"]
        h5_file = h5py.File(data_file_name, "r")

        plot_2d(args, font_size, h5_file)

        h5_file.close()
    elif args["conv-time-1d"]:
        coarse_data = args["<coarse-data>"]
        medium_data = args["<medium-data>"]
        fine_data = args["<fine-data>"]

        h5_coarse = h5py.File(coarse_data, "r")
        h5_medium = h5py.File(medium_data, "r")
        h5_fine = h5py.File(fine_data, "r")

        conv_time_1d(args, font_size, h5_coarse, h5_medium, h5_fine)

        h5_coarse.close()
        h5_medium.close()
        h5_fine.close()


if __name__ == '__main__':
    arguments = docopt(__doc__, version="RKAB Plot 1.0")
    main(arguments)
