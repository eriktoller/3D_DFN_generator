from DFN_classes import *

import matplotlib.pyplot as plt
import math
import pickle

from matplotlib.axes import Axes
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import datetime



def find_intersections(dfn0: DFN) -> List:
    """
    Program that computs the intersections of circular fractures. It also adds them to the fractures in the DFN-class
    file.

    :param dfn0: DFN with alla fractures
    :return: list with intersecting fractures and there locations
    """
    # Get som reoccurring variables
    nf = dfn0.num_frac()
    frac = dfn0.fractures

    # make arrays with all possible intersections
    range_temp = np.array(range(nf))

    # Remove self from list
    int_list = []
    for row in range(nf):
        temp = np.delete(range_temp, row)
        int_list.append(temp)

    # Remove fracture where the circular boundaries will not intersect
    cnt = -1
    for frac_it in frac:
        cnt = cnt + 1
        for col in int_list[cnt]:
            dist = (frac_it.zc3d[0] - frac[col].zc3d[0]) ** 2 + (frac_it.zc3d[1] - frac[col].zc3d[1]) ** 2 + (
                    frac_it.zc3d[2] - frac[col].zc3d[2]) ** 2
            dist = math.sqrt(dist)
            dist = frac_it.r + frac[col].r - dist

            if dist < 0:
                int_list[cnt] = np.delete(int_list[cnt], int_list[cnt] == col)

    # Find the intersections
    cnt = -1
    z_int0 = []
    for frac_it in frac:
        cnt = cnt + 1
        z_int_temp = []
        for col in int_list[cnt]:
            n_1 = frac_it.n3
            n_2 = frac[col].n3
            p_1 = frac_it.zc3d
            p_2 = frac[col].zc3d
            n_int = np.cross(n_1, n_2)

            a = np.matrix(np.array([[2, 0, 0, n_1[0], n_2[0]],
                                    [0, 2, 0, n_1[1], n_2[1]],
                                    [0, 0, 2, n_1[2], n_2[2]],
                                    [n_1[0], n_1[1], n_1[2], 0, 0],
                                    [n_2[0], n_2[1], n_2[2], 0, 0]]))
            b4 = p_1[0] * n_1[0] + p_1[1] * n_1[1] + p_1[2] * n_1[2]
            b5 = p_2[0] * n_2[0] + p_2[1] * n_2[1] + p_2[2] * n_2[2]
            b = np.matrix(np.array([[2.0 * p_1[0]], [2.0 * p_1[1]], [2.0 * p_1[2]], [b4], [b5]]))

            x = np.linalg.solve(a, b)
            p0 = np.squeeze(np.asarray(x[0:3]))
            b = np.matrix(np.array([[2 * p_2[0]], [2 * p_2[1]], [2 * p_2[2]], [b4], [b5]]))
            x = np.linalg.solve(a, b)
            p01 = np.squeeze(np.asarray(x[0:3]))

            z_0 = point_3d_to_2d(p0 - n_int, frac_it)
            z_1 = point_3d_to_2d(p0 + n_int, frac_it)

            # Check if line intersects
            z0 = point_3d_to_2d(p0, frac_it)
            z01 = point_3d_to_2d(p01, frac[col])
            dr0 = np.abs(z0) - frac_it.r
            dr01 = np.abs(z01) - frac[col].r
            dr02 = np.linalg.norm(p01 - p0) - frac_it.r - frac[col].r

            if any([dr0 > 0, dr01 > 0, dr02 > 0]):
                int_list[cnt] = np.delete(int_list[cnt], int_list[cnt] == col)
                int_list[col] = np.delete(int_list[col], int_list[col] == cnt)
            else:
                if np.real(z_1 - z_0) == 0:
                    x1 = np.real(z_0)
                    x2 = x1
                    b0 = frac_it.r ** 2 - x1 ** 2
                    y1 = np.sqrt(b0)
                    y2 = - y1
                    zi0 = x1 + 1j * y1
                    zi1 = x2 + 1j * y2
                    z_int_temp.append([zi0, zi1])
                else:
                    m = np.imag(z_1 - z_0) / np.real(z_1 - z_0)
                    x0 = np.imag(z_0) - m * np.real(z_0)

                    # Calculate the intersection points
                    a0 = (1 + m ** 2)
                    a1 = 2 * m * x0
                    a2 = x0 ** 2 - frac_it.r ** 2
                    a3 = (a1 / (2 * a0)) ** 2 - a2 / a0

                    if a3 < 0:
                        int_list[cnt] = np.delete(int_list[cnt], int_list[cnt] == col)
                        int_list[col] = np.delete(int_list[col], int_list[col] == cnt)
                    else:
                        x1 = - a1 / (2 * a0) + np.sqrt(a3)
                        x2 = - a1 / (2 * a0) - np.sqrt(a3)
                        y1 = m * x1 + x0
                        y2 = m * x2 + x0
                        zi0 = x1 + 1j * y1
                        zi1 = x2 + 1j * y2
                        z_int_temp.append([zi0, zi1])

        z_int0.append(z_int_temp)

    # Go through the intersections to select the smallest lines
    cnt = -1
    for frac_it in frac:
        cnt = cnt + 1
        z_int_temp = []
        cnt1 = -1
        for col2 in range(len(int_list[cnt])):
            cnt1 = cnt1 + 1
            col = int_list[cnt][cnt1]
            pos = np.where(int_list[col] == cnt)[0][0]
            # Get teh corresponding 3D points to all intersection points
            z3d0_it = point_2d_to_3d(z_int0[cnt][cnt1][0], frac_it)
            z3d1_it = point_2d_to_3d(z_int0[cnt][cnt1][1], frac_it)
            z3d0 = point_2d_to_3d(z_int0[col][pos][0], frac[col])
            z3d1 = point_2d_to_3d(z_int0[col][pos][1], frac[col])

            # Calculate their position in the fracture plane for the first fracture
            z_00 = z_int0[cnt][cnt1][0]  # original intersection point
            z_11 = z_int0[cnt][cnt1][1]  # original intersection point
            z_0 = point_3d_to_2d(z3d0, frac_it)  # other fractures intersection point
            z_1 = point_3d_to_2d(z3d1, frac_it)  # other fractures intersection point
            # Calculate the distance to the center point
            dz_00 = np.abs(z_00)
            dz_11 = np.abs(z_11)
            df1 = min([dz_00, dz_11])
            dz_0 = np.abs(z_0)
            dz_1 = np.abs(z_1)

            # Calculate their position in the fracture plane for the second fracture
            z_20 = z_int0[col][pos][0]
            z_31 = z_int0[col][pos][1]
            z_2 = point_3d_to_2d(z3d0_it, frac[col])
            z_3 = point_3d_to_2d(z3d1_it, frac[col])
            dz_20 = np.abs(z_20)
            dz_31 = np.abs(z_31)
            df2 = min([dz_20, dz_31])
            dz_2 = np.abs(z_2)
            dz_3 = np.abs(z_3)
            if all([dz_0 < df1, dz_1 < df1]):
                z_int0[cnt][cnt1][0] = z_0
                z_int0[cnt][cnt1][1] = z_1
            elif all([dz_2 < df2, dz_3 < df2]):
                z_int0[col][pos][0] = z_2
                z_int0[col][pos][1] = z_3
            else:
                if dz_0 <= df1:
                    z_int0[cnt][cnt1][0] = z_0
                    z_int0[col][pos][0] = z_20
                if dz_1 <= df1:
                    z_int0[cnt][cnt1][0] = z_1
                    z_int0[col][pos][0] = z_31
                if dz_2 <= df2:
                    z_int0[cnt][cnt1][1] = z_00
                    z_int0[col][pos][1] = z_2
                if dz_3 <= df2:
                    z_int0[cnt][cnt1][1] = z_11
                    z_int0[col][pos][1] = z_3
            if all([dz_2 > df2, dz_3 > df2, dz_0 > df1, dz_1 > df1]):
                # Remove if circle do not intersect
                int_list[cnt] = np.delete(int_list[cnt], int_list[cnt] == col)
                int_list[col] = np.delete(int_list[col], int_list[col] == cnt)
                del z_int0[cnt][cnt1]
                del z_int0[col][pos]
                cnt1 = cnt1 - 1

    # Add intersection data to DFN
    cnt = -1
    for frac_it in frac:
        cnt = cnt + 1
        dfn0.fractures[cnt].int_list = int_list[cnt]
        dfn0.fractures[cnt].z_int = z_int0[cnt]

    return [int_list, z_int0]


def plot_circle_3d(zc: Circle3D, num: int, ax: Axes) -> None:
    """
    Program that plots a circular plane in 3D for a given center and Cartesian coordinate system.

    :param zc: 3D circle class
    :param num: number of plotting points
    :param ax: (sub-)plot reference
    :return: plot on ax
    """

    # Make array from 0 to 2*pi with num number of points
    theta = np.linspace(0, 2 * math.pi, num)

    # Calculate the 3D points
    xyz = zc.zc3d + np.array([np.cos(theta)]).T * zc.r * zc.n1 + np.array([np.sin(theta)]).T * zc.r * zc.n2

    # Plot the circle in 3D
    ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], linestyle=zc.ls, marker=zc.marker, color=zc.color, linewidth=zc.ls_width)

    # fill the circles
    verts = [list(zip(xyz[:, 0], xyz[:, 1], xyz[:, 2]))]
    collection = Poly3DCollection(verts, linewidths=1, alpha=zc.alpha)
    collection.set_facecolor(zc.face_color)
    ax.add_collection3d(collection)


def plot_frac_coord(zc: Circle3D, scale: float, ls_width: float, ax: Axes) -> None:
    """
    Program that plots the local Cartesian coordinate system for a fracture.

    :param zc: 3D circle class
    :param scale: scale of orientation vectors (1 = unit vector)
    :param ls_width: line width
    :param ax: (sub-)plot reference
    :return: plot on ax
    """
    # first unit vector
    zcp0 = zc.zc3d
    zcp1 = zc.zc3d + zc.n1 * scale
    ax.plot([zcp0[0], zcp1[0]], [zcp0[1], zcp1[1]], [zcp0[2], zcp1[2]], '-', color=zc.color, linewidth=ls_width)
    # second unit vector
    zcp0 = zc.zc3d
    zcp1 = zc.zc3d + zc.n2 * scale
    ax.plot([zcp0[0], zcp1[0]], [zcp0[1], zcp1[1]], [zcp0[2], zcp1[2]], '-', color=zc.color, linewidth=ls_width)
    # third unit vector
    zcp0 = zc.zc3d
    zcp1 = zc.zc3d + zc.n3 * scale
    ax.plot([zcp0[0], zcp1[0]], [zcp0[1], zcp1[1]], [zcp0[2], zcp1[2]], '-', color=zc.color, linewidth=ls_width)


def point_3d_to_2d(pnt: np.ndarray, circ: Circle3D) -> np.ndarray:
    """
    Program that maps a point in 3D to a point in the 2D plane of a fracture.

    :param pnt: point in 3D sapce
    :param circ: circle in 3d class
    :return: 2D location of 3D point in the plane of the circle
    """
    if all(pnt == circ.zc3d):
        z = 0 + 1j
    else:
        dist = pnt - circ.zc3d
        d = np.linalg.norm(dist)
        theta = np.arccos(np.dot(dist, circ.n1) / d)

        if np.dot(dist, circ.n2) < 0:
            theta = -theta

        x = d * np.cos(theta)
        y = d * np.sin(theta)

        z = x + 1j * y

    return z


def point_2d_to_3d(pnt: np.ndarray, zc: Circle3D) -> np.ndarray:
    """
    Program that maps a point in 2D plane of a fracture to the 3D space.

    :param pnt: point or array in 2D plane (complex)
    :param zc: circle in 3d class
    :return: 3D location
    """
    pnt2 = np.where((pnt == 0 + 0j), 1, pnt)
    theta = np.log(pnt2).imag
    dist = np.abs(pnt)
    xyz = zc.zc3d + np.array([np.cos(theta) * dist]).T * zc.n1 + np.array([np.sin(theta) * dist]).T * zc.n2

    return xyz


def plot_point_3d(pnt: np.ndarray, zc: Circle3D, ax: Axes, ls0: str, marker0: str, color0: str, ls_width: str,
                  f_style: str) -> None:
    """
    Program that plots a circular plane in 3D for a given center and Cartesian coordinate system.

    :param pnt: point or array in 2D plane (complex)
    :param zc: circle in 3d class
    :param ax: (sub-)plot reference
    :param ls0: line style
    :param marker0: marker style
    :param color0: line and marker color
    :param ls_width: line width
    :param f_style: fill style of marker
    :return: plot on ax
    """
    pnt2 = np.where((pnt == 0 + 0j), 1, pnt)
    theta = np.log(pnt2).imag
    dist = np.abs(pnt)
    xyz = zc.zc3d + np.array([np.cos(theta) * dist]).T * zc.n1 + np.array([np.sin(theta) * dist]).T * zc.n2

    ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], linestyle=ls0, marker=marker0, color=color0, linewidth=ls_width,
            fillstyle=f_style)


def axis_equal3d(ax: Axes) -> None:
    """
    Make the axis equal for 3D plot

    :param ax: (sub-)plot reference
    :return: 3D plot with equal axis
    """
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:, 1] - extents[:, 0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize / 2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


def generate_fractures(bbox: Box, num_frac: int, r_min: float, r_max: float, dfn_name: str) -> DFN:
    """
    Generator for DFN

    :param bbox: bounding box for generator Box-class
    :param num_frac: number of fractures to generate
    :param r_min: minimum radius of fractures
    :param r_max: maximum radius of fractures
    :param dfn_name: name of the DFN
    :return: DFN as a DFN dataclass
    """
    # Generate points
    x1 = np.random.uniform(low=bbox.x1[0], high=bbox.x1[1], size=num_frac)
    x2 = np.random.uniform(low=bbox.x2[0], high=bbox.x2[1], size=num_frac)
    x3 = np.random.uniform(low=bbox.x3[0], high=bbox.x3[1], size=num_frac)
    zc = np.c_[x1, x2, x3]

    # Generate normals
    n11 = np.random.uniform(low=-1, high=1, size=num_frac)
    n21 = np.random.uniform(low=-1, high=1, size=num_frac)
    n31 = np.random.uniform(low=-1, high=1, size=num_frac)
    n1 = np.c_[n11, n21, n31]
    for n in range(num_frac):
        mag = np.linalg.norm(n1[n])
        n1[n][0] = n1[n][0] / mag
        n1[n][1] = n1[n][1] / mag
        n1[n][2] = n1[n][2] / mag

    n11 = np.random.uniform(low=-1, high=1, size=num_frac)
    n21 = np.random.uniform(low=-1, high=1, size=num_frac)
    n31 = np.random.uniform(low=-1, high=1, size=num_frac)
    n2temp = np.c_[n11, n21, n31]
    n2 = np.cross(n1, n2temp)
    for n in range(num_frac):
        mag = np.linalg.norm(n2[n])
        n2[n][0] = n2[n][0] / mag
        n2[n][1] = n2[n][1] / mag
        n2[n][2] = n2[n][2] / mag

    n3 = np.cross(n1, n2)
    for n in range(num_frac):
        mag = np.linalg.norm(n3[n])
        n3[n][0] = n3[n][0] / mag
        n3[n][1] = n3[n][1] / mag
        n3[n][2] = n3[n][2] / mag

    # Generate radii
    r = np.random.uniform(low=r_min, high=r_max, size=num_frac)

    # Store fractures in DFN dataclass
    zc_it = []
    for jj in range(num_frac):
        zc_it.append(Circle3D(zc3d=zc[jj], r=r[jj], n1=n1[jj], n2=n2[jj], n3=n3[jj], int_list=[], z_int=[]))
    dfn_gen = DFN(fractures=zc_it, name=dfn_name)

    return dfn_gen


print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print('Program started')
print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

# ct stores current time
ct = datetime.datetime.now()
print("Program started at: ", ct)
print(" ")

# Generate DFN
bbox = Box(x1=[-1, 1], x2=[-1, 1], x3=[-1, 1])
num_frac = 3
r_min = 1.0
r_max = 1.9
DFN_rand = generate_fractures(bbox, num_frac, r_min, r_max, 'DFN_random')

# Find the intersections
[int_list, z_int] = find_intersections(DFN_rand)


# Print the DFN information
name = str(ct.date()) + '_' + DFN_rand.name + '_data_0.pkl'
print('DFN information:')
print('  Name: ' + str(DFN_rand.name))
print('  Data file: ' + name)
print('  Number of fracture: ' + str(DFN_rand.num_frac()))
print('  Number of intersections:' + str(DFN_rand.num_int()))

# Initiate the figure
fig = plt.figure()
ax0 = fig.add_subplot(projection='3d')
ax0.set_xlabel('$x_1$-axis')
ax0.set_ylabel('$x_2$-axis')
ax0.set_zlabel('$x_3$-axis')
ax0.set_aspect('equal')

# Plot the DFN
for zcc in DFN_rand.fractures:
    plot_circle_3d(zcc, 100, ax0)
    plot_frac_coord(zcc, .2, .5, ax0)


# Plot intersections
for zcc in DFN_rand.fractures:
    for zi in zcc.z_int:
        plot_point_3d(zi, zcc, ax0, '-', 'o', 'black', '.5', 'none')

ct2 = datetime.datetime.now()
print("")
print("Program finished at: ", ct2)
exe_time = ct2 - ct
print("           and took: ", exe_time)


# Save DFN as pickle
import glob

files = glob.glob(('data/'+name))
cnt_name = 0
while len(files) > 0:
    cnt_name = cnt_name + 1
    m = 4 + len(str(cnt_name-1))
    name = name[0:-m] + str(cnt_name) + '.pkl'
    files = glob.glob(('data/'+name))

output = open(('data/'+name), 'wb')
data1 = {'DFN': DFN_rand,
         'fig': fig,
         'start_time': ct,
         'end_time': ct2,
         'exe_time': exe_time}
pickle.dump(data1, output)

# Set the axis as equal and show the plot
axis_equal3d(ax0)
plt.show()
# plt.savefig('figure1.pdf')
