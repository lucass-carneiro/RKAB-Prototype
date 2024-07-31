import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt


font_size = 18
line_thickness = 2.0
line_color = "black"

mpl.rcParams["mathtext.fontset"] = "cm"
mpl.rcParams["font.family"] = "Latin Modern Roman"
mpl.rcParams["xtick.labelsize"] = font_size
mpl.rcParams["ytick.labelsize"] = font_size

f = h5py.File("output.h5", "r")

x = f["grid/x_coords"][:]
y = f["state/Phi_0010"][:]

plt.close("all")
plt.plot(x, y)
plt.ylim(-1.0, 1.0)
plt.tight_layout()
plt.show()
