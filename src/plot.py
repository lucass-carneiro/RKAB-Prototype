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

plt.close("all")
# plt.plot(x, f["rhs/Phi_rhs_0001"][:], color="red")
# plt.plot(x, f["rhs/Pi_rhs_0001"][:], color="black")
# plt.plot(x, f["rhs/Dx_rhs_0001"][:], color="blue")
plt.plot(x, f["state/Phi_0080"][:], color="red")
plt.plot(x, f["state/Pi_0080"][:], color="black")
plt.plot(x, f["state/Dx_0080"][:], color="blue")
plt.xlim(-1.0, 1.0)
# plt.ylim(-8.0, 8.0)
plt.tight_layout()
plt.show()
