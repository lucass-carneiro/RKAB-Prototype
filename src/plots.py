import domain as dom

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import adios2


font_size = 18
line_thickness = 2.0
line_color = "black"

mpl.rcParams["mathtext.fontset"] = "cm"
mpl.rcParams["font.family"] = "Latin Modern Roman"
mpl.rcParams["xtick.labelsize"] = font_size
mpl.rcParams["ytick.labelsize"] = font_size


def plot_time_steps(domain: dom.Domain, stream: adios2.Stream):
    levels = np.linspace(-1.0, 1.0, 100)

    for _ in stream.steps():
        iter = stream.current_step()
        print(f"Plotting iteration {iter}")
        u = stream.read("u")

        x_data = []
        y_data = []
        z_data = []

        for p in domain.interior_points():
            i, j, k = p.idx()
            x, y, z = p.pos()

            if (np.isclose(z, 0.0)):
                x_data.append(x)
                y_data.append(y)
                z_data.append(u[i, j, k])

        plt.close("all")
        plt.tricontourf(x_data, y_data, z_data, levels=levels, cmap="RdBu")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(f"iteration_{iter}_.png")
