import numpy as np

import matplotlib.pyplot as plt

import os
import concurrent.futures

import logging
logger = logging.getLogger(__name__)


def sw_Phi_3d(A, kx, ky, kz, t, x, y, z):
    omega = np.sqrt(kx * kx + ky * ky + kz * kz)
    return A*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)*np.sin(2*ky*np.pi*y)*np.sin(2*kz*np.pi*z)


def sw_Pi_3d(A, kx, ky, kz, t, x, y, z):
    omega = np.sqrt(kx * kx + ky * ky + kz * kz)
    return -2*A*omega*np.pi*np.sin(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)*np.sin(2*ky*np.pi*y)*np.sin(2*kz*np.pi*z)


def sw_Dx_3d(A, kx, ky, kz, t, x, y, z):
    omega = np.sqrt(kx * kx + ky * ky + kz * kz)
    return 2*A*kx*np.pi*np.cos(2*omega*np.pi*t)*np.cos(2*kx*np.pi*x)*np.sin(2*ky*np.pi*y)*np.sin(2*kz*np.pi*z)


def sw_Dy_3d(A, kx, ky, kz, t, x, y, z):
    omega = np.sqrt(kx * kx + ky * ky + kz * kz)
    return 2*A*ky*np.pi*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)*np.cos(2*ky*np.pi*y)*np.sin(2*kz*np.pi*z)


def sw_Dz_3d(A, kx, ky, kz, t, x, y, z):
    omega = np.sqrt(kx * kx + ky * ky + kz * kz)
    return 2*A*kz*np.pi*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)*np.sin(2*ky*np.pi*y)*np.cos(2*kz*np.pi*z)


def plot_gfs_3d(x, y, data, levels, font_size, prefix, name, iteration, iteration_string, t, path):
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
            f"3D_{prefix}_{name}_it_{iteration_string}.png"
        )
    )


def plot_expected_3d(x, y, z, levels, font_size, prefix, name, iteration, iteration_string, t, path):
    plt.close("all")
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
            f"3D_{prefix}_{name}_it_{iteration_string}.png"
        )
    )


def plot_3d(args, font_size, h5_file):
    last_iter = h5_file.attrs["last_iter"]
    dt = h5_file.attrs["dt"]

    x = h5_file["grid/x_coords"][:]
    y = h5_file["grid/y_coords"][:]
    z = h5_file["grid/z_coords"][:]

    slice_value = float(args["--slice-value"])
    slice_idx = np.where(np.isclose(z, slice_value))

    slice_value = z[slice_idx]
    logger.info(f"Actual slice used: {slice_value}")

    if (np.size(slice_idx) == 0):
        logger.error(f"Slice value {slice_value} of coordinate z not found")
        return

    slice_idx = slice_idx[0][0]
    slice_value = z[slice_idx]

    logger.info(
        f"Slicing coordinate z at value {slice_value}, index {slice_idx}"
    )

    X, Y = np.meshgrid(x, y)

    A = h5_file.attrs["A"]
    kx = h5_file.attrs["kx"]
    ky = h5_file.attrs["ky"]
    kz = h5_file.attrs["kz"]

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
        "Dz",
    ]

    rhs_gfs = [
        "Phi_rhs",
        "Pi_rhs",
        "Dx_rhs",
        "Dy_rhs",
        "Dz_rhs",
    ]

    plots_dir = "3d_plots"
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
    err_level_array = np.linspace(0.0, 1.0, endpoint=True, num=101)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in iteration_range:
            logger.info(f"Dispatching plot of iteration {i}")

            iteration_string = f"{i:04}"
            t = i * dt

            # Plot State
            for gf in gfs:
                path = f"state/{gf}_{iteration_string}"
                data = h5_file[path][:, :, slice_idx]

                if gf == "Phi":
                    levels = level_array
                else:
                    levels = 100

                executor.submit(
                    plot_gfs_3d,
                    X,
                    Y,
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
                data = h5_file[path][:, :, slice_idx]

                executor.submit(
                    plot_gfs_3d,
                    X,
                    Y,
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
                data = h5_file[path][:, :, slice_idx]

                z = np.abs(
                    eval(f"sw_{gf}_3d")(
                        A,
                        kx,
                        ky,
                        kz,
                        t,
                        x,
                        y,
                        slice_value
                    ) - data
                )

                executor.submit(
                    plot_expected_3d,
                    X,
                    Y,
                    z,
                    100,
                    font_size,
                    "state",
                    gf,
                    i,
                    iteration_string,
                    t,
                    expected_dir
                )

        logger.info(f"Waiting for plot workers to finish")
