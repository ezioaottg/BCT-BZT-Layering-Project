import os
import temul.api as tml
import atomap.api as am
import hyperspy.api as hs
import matplotlib.pyplot as plt
from temul.topotem import (combine_atom_deviations_from_zone_axes, plot_polarisation_vectors,get_average_polarisation_in_regions_square)

import atomap
import numpy as np

def GetLatticeA(image, separation):
    atom_positions_A = am.get_atom_positions(image, separation)
    sublattice_A = am.Sublattice(atom_positions_A, image = image.data)

    sublattice_A.find_nearest_neighbors()
    sublattice_A.refine_atom_positions_using_center_of_mass()
    sublattice_A.refine_atom_positions_using_2d_gaussian()

    sublattice_A.construct_zone_axes()

    sublattice_A_image_float64 = sublattice_A.image.astype(np.float64)

    mat_noA = atomap.tools.remove_atoms_from_image_using_2d_gaussian(sublattice_A_image_float64, sublattice_A)
    mat_noA = mat_noA.astype(np.uint8)

    image_noA = hs.signals.BaseSignal(mat_noA)

    return sublattice_A, image_noA

#/Users/dizhang/Desktop/Polarization map fitting/Polarization map fitting/IFFT of Non-regid alignment 002-2_bw copy.png
# Load the image
path = os.getcwd()
file = r'IFFT of Non-regid alignment 002-2_bw copy.png'
image = hs.load(os.path.join(path, file))

# Resize the image to desired dimensions
image_resized = image.rebin((934, 934))

# Visualize and filter the image
tml.visualise_dg_filter(image_resized)
filtered_image = tml.double_gaussian_fft_filter(image_resized, 30, 100)

# Shift all values to be positive
min_val = filtered_image.data.min()
if min_val < 0:
    filtered_image.data -= min_val

# Plot the images
image_resized.plot()
filtered_image.plot()
plt.show()

image = filtered_image
separation = 29 # set based on lower limit of good fit from peaks^, 38 is being used for 2-2 image

sublattice_A, image_noA = GetLatticeA(image, separation)
sublattice_A.plot()
sublattice_A.plot_planes()
# plt.show()

zone_axis_001 = sublattice_A.zones_axis_average_distances[3]
B_positions = sublattice_A.find_missing_atoms_from_zone_vector(zone_axis_001, vector_fraction= 0.5)

sublattice_A_image_float64 = sublattice_A.image.astype(np.float64)

mat_noA = atomap.tools.remove_atoms_from_image_using_2d_gaussian(sublattice_A_image_float64, sublattice_A, percent_to_nn= 0.6) #0.3 used for 2-2, .75
mat_noA = mat_noA.astype(np.uint8)

image_noA = hs.signals.BaseSignal(mat_noA)
sublattice_B = am.Sublattice(B_positions, mat_noA, color='blue')
sublattice_B.construct_zone_axes()
sublattice_B.refine_atom_positions_using_center_of_mass()
sublattice_B.refine_atom_positions_using_2d_gaussian()

atom_lattice = am.Atom_Lattice(image=image.data, name='test', sublattice_list=[sublattice_A, sublattice_B], fix_negative_values=True)
atom_lattice.plot()
# plt.show()

atom_lattice = am.Atom_Lattice(image = image, name = 'test', sublattice_list=[sublattice_A, sublattice_B], fix_negative_values=True)
atom_lattice.plot()

sublattice_A.construct_zone_axes()
za0, za1 = sublattice_A.zones_axis_average_distances[0:2]
s_p = sublattice_A.get_polarization_from_second_sublattice(za0, za1, sublattice_B, color='blue')
vector_list = s_p.metadata.vector_list
# x, y = [i[0] for i in vector_list], [i[1] for i in vector_list]
# u, v = [i[2] for i in vector_list], [i[3] for i in vector_list]
#
# magnitude = np.sqrt(np.array(u)**2 + np.array(v)**2)
# u_normalized = np.array(u) / magnitude
# v_normalized = np.array(v) / magnitude


#
# sublatticeA = atom_lattice.sublattice_list[0]
# sublatticeA.find_nearest_neighbors()
# _ = sublatticeA.refine_atom_positions_using_center_of_mass()
sublattice_A.construct_zone_axes()


image = sublattice_A.image[0:930]
x, y, u, v = combine_atom_deviations_from_zone_axes(sublattice_B,   save=None)
ax = plot_polarisation_vectors(x, y, u, v, image=image_resized, save=None, color='y', overlay=True, monitor_dpi=50, title='Actual Vector Arrows')
ax.plot()
plt.show()
coords = get_average_polarisation_in_regions_square(x, y, u, v, image=image, divide_into=24)
x_new, y_new, u_new, v_new = coords
ax = plot_polarisation_vectors(x_new, y_new, u_new, v_new, image=image_resized, color='y', overlay=True, monitor_dpi=50,title='Averaged Vector Arrows', save=None)
ax.plot()
plt.show()