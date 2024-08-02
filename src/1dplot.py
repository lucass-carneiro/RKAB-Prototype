import numpy as np

import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt


def Phi(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return A*np.cos(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)


def Pi(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return -2*A*omega*np.pi*np.sin(2*omega*np.pi*t)*np.sin(2*kx*np.pi*x)


def Dx(A, kx, t, x):
    omega = np.sqrt(kx * kx)
    return 2*A*kx*np.pi*np.cos(2*omega*np.pi*t)*np.cos(2*kx*np.pi*x)


def plot_gfs(h5_file, prefix, gfs, iteration, iteration_string, dt, x, font_size):
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

    plt.savefig(f"1D_{prefix}_it_{iteration}.pdf")


def plot_expected(h5_file, prefix, gfs, iteration, iteration_string, dt, x, font_size):
    A = h5_file.attrs["A"]
    kx = h5_file.attrs["kx"]
    t = iteration * dt

    plt.close("all")

    for (name, color) in gfs:
        gf_path = f"{prefix}/{name}_{iteration_string}"
        expected = eval(name)(A, kx, t, x)
        data = h5_file[gf_path][:] - expected
        plt.plot(x, data, color=color, label=f"True{name} - {name}")

    plt.xlim(-1.0, 1.0)

    plt.xlabel("$x$", size=font_size)
    plt.ylabel("$f(x)$", size=font_size)

    plt.title(f"{prefix} at iteration {iteration}, $t = {t}$")
    plt.legend()

    plt.tight_layout()

    plt.savefig(f"1D_expected_{prefix}_it_{iteration}.pdf")


def main():
    iteration = 1
    data_file_name = "1d_output.h5"

    font_size = 18

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

    h5_file = h5py.File(data_file_name, "r")

    dt = h5_file.attrs["dt"]
    x = h5_file["grid/x_coords"][:]

    iteration_string = f"{iteration:04}"

    plot_gfs(h5_file, "state", gfs, iteration,
             iteration_string, dt, x, font_size)
    plot_gfs(h5_file, "rhs", rhs_gfs, iteration,
             iteration_string, dt, x, font_size)

    plot_expected(h5_file, "state", gfs, iteration,
                  iteration_string, dt, x, font_size)

    h5_file.close()


main()
