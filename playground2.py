import os
import temul.api as tml
import atomap.api as am
import hyperspy.api as hs
import matplotlib.pyplot as plt
from temul.topotem import (combine_atom_deviations_from_zone_axes, plot_polarisation_vectors,get_average_polarisation_in_regions_square)
import temul.dummy_data as dd

import atomap
import numpy as np



sublattice = dd.get_polarised_single_sublattice()

sublatticeA = sublattice


sublatticeA.construct_zone_axes(atom_plane_tolerance=1)
sublattice.plot_planes()





x,y,u,v = combine_atom_deviations_from_zone_axes(
    sublatticeA, save=None)
#ax = plot_polarisation_vectors(x, y, u, v, save=None,
    #image=sublatticeA.image)

#ax.plot()
#x, y, u, v = tml.atom_deviation_from_straight_line_fit(
    #sublattice, 3,4)

tml.plot_polarisation_vectors(x, y, u, v, image=sublatticeA.image,
                              unit_vector=False, save=None,
                              plot_style='vector', color = 'r',
                              overlay=True, monitor_dpi=50)

plt.show()