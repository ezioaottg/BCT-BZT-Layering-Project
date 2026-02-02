import os
import temul.api as tml
import temul
import atomap.api as am
import hyperspy.api as hs
import matplotlib.pyplot as plt
from temul.topotem import (combine_atom_deviations_from_zone_axes, plot_polarisation_vectors,get_average_polarisation_in_regions_square)
from atomap.tools import remove_atoms_from_image_using_2d_gaussian
import atomap
import numpy as np
from matplotlib import ticker
import matplotlib.colors as colors



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

path = os.getcwd()
file = r'IFFT of 1064A_200kV_15MX_JEOL  ADF_0076-2 cropped 2_bw-1.png'
image = hs.load(os.path.join(path, file))

# Resize the image to desired dimensions & perform a FFT filter.
image_resized = image.rebin((934, 934))
filtered_image = tml.double_gaussian_fft_filter(image_resized, 88, 225)


# Shift all values to be positive
min_val = filtered_image.data.min()
if min_val < 0:
    filtered_image.data -= min_val


image = filtered_image


separation = 15

sublattice_A, image_noA = GetLatticeA(image, separation)

zone_axis_001 = sublattice_A.zones_axis_average_distances[3]


B_positions = sublattice_A.find_missing_atoms_from_zone_vector(zone_axis_001)
image_without_A = remove_atoms_from_image_using_2d_gaussian(sublattice_A.image, sublattice_A)

sublattice_B = am.Sublattice(B_positions, image_without_A, color='blue', fix_negative_values = True)
sublattice_B.construct_zone_axes()
sublattice_B.refine_atom_positions_using_center_of_mass()
sublattice_B.refine_atom_positions_using_2d_gaussian()


atom_lattice = am.Atom_Lattice(image=filtered_image.data, name='test', sublattice_list=[sublattice_A, sublattice_B])


za0, za1 = sublattice_A.zones_axis_average_distances[0:2]
s_p = sublattice_A.get_polarization_from_second_sublattice(
    za0, za1, sublattice_B)

vector_list = s_p.metadata.vector_list
x, y = [i[0] for i in vector_list], [i[1] for i in vector_list]
u, v = [i[2] for i in vector_list], [i[3] for i in vector_list]
sampling, units =  0.05, 'nm'

#tml.plot_polarisation_vectors(x, y, u, v, image=atom_lattice.image,
                          #sampling=sampling, units=units,
                          #unit_vector=False, save=None, scalebar=False,
                          #plot_style='vector', color='r',
                          #overlay=False, monitor_dpi=45)


#x, y, u, v = combine_atom_deviations_from_zone_axes(sublattice_B,   save=None)

#ax = plot_polarisation_vectors(x, y, u, v, image=image_resized, save=None, color='y', overlay=True, monitor_dpi=50, title='Actual Vector Arrows')
#ax.plot()



#x, y, u, v = tml.atom_deviation_from_straight_line_fit(
    #sublattice_A, 0,14)


tml.plot_polarisation_vectors(x, y, u, v, image=atom_lattice.image,
                          sampling=3.0321, units='nm', monitor_dpi=50,
                          unit_vector=False, plot_style='colormap',
                          overlay=True, save=None, cmap='viridis',
                          scalebar=True, ticks = ticker.MaxNLocator(nbins=5))

ratio_map = tml.ratio_of_lattice_spacings(sublattice_B, 0, 1,
                units="nm", sampling=0.1)

plt.show()
