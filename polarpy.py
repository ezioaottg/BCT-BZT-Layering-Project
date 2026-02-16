import hyperspy
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import ScalarMappable
from decimal import Decimal
import colorcet as cc
from matplotlib_scalebar.scalebar import ScaleBar
from temul.signal_plotting import (
    get_polar_2d_colorwheel_color_list,
    _make_color_wheel)
from temul.topotem import polarisation



def plot_polarisation_maps(
        x, y, u, v, image, sampling=None, units='pix',
        plot_style='vector', overlay=True, unit_vector=False,
        vector_rep='magnitude', degrees=False, angle_offset=None,
        save='polarisation_image', title="", color='yellow',
        cmap=None, cmap2 = None, alpha=1.0, image_cmap='gray', monitor_dpi=96,
        no_axis_info=True, invert_y_axis=True, ticks=None, scalebar=False,
        antialiased=False, levels=20, remove_vectors=False,
        quiver_units='width', pivot='middle', angles='xy',
        scale_units='xy', scale=None, headwidth=3.0, headlength=5.0,
        headaxislength=4.5, width=None, minshaft=1, minlength=1):
    
    if isinstance(image, np.ndarray):
        pass
    elif isinstance(image, hyperspy._signals.signal2d.Signal2D):
        sampling = image.axes_manager[-1].scale
        units = image.axes_manager[-1].units
        image = image.data
    else:
        raise ValueError("``image`` must be a 2D numpy array or 2D Hyperspy "
                         "Signal")
    
    u, v = np.array(u), np.array(v)

    if sampling is not None:
        u, v = u * sampling, v * sampling

    # get the magnitude or angle representation
    if vector_rep == "magnitude":
        vector_rep_val = polarisation.get_vector_magnitudes(u, v)
        map_rep_val = polarisation.get_angles_from_uv(u, -v, degrees=degrees, angle_offset=angle_offset)
    elif vector_rep == "angle":
        # -v because in STEM the origin is top left
        vector_rep_val = polarisation.get_angles_from_uv(u, -v, degrees=degrees,
                                            angle_offset=angle_offset)
        map_rep_val = polarisation.get_vector_magnitudes(u,v)
        
    vector_label = polarisation.angle_label(
        vector_rep=vector_rep, units=units, degrees=degrees)
    
    if vector_rep == "magnitude":
        map_label = polarisation.angle_label(vector_rep = 'angle', units = units, degrees = degrees)
    elif vector_rep == 'angle':
        map_label = polarisation.angle_label(vector_rep = 'magnitude', units = units, degrees = degrees)

    if plot_style == "polar_colorwheel":
        color_list = get_polar_2d_colorwheel_color_list(u, -v)

    # change all vector magnitudes to the same size
    if unit_vector:
        u_norm = u / np.sqrt(u ** 2.0 + v ** 2.0)
        v_norm = v / np.sqrt(u ** 2.0 + v ** 2.0)
        u = u_norm
        v = v_norm

    # setting up norm and cmap for colorbar scalar mappable
    if vector_rep == "angle":
        min_map_val = np.min(map_rep_val)
        max_map_val = np.max(map_rep_val) + 0.00000001
        if degrees:
            min_val, max_val = -180, 180 + 0.0001  # fix display issues
        elif not degrees:
            min_val, max_val = -np.pi, np.pi
    elif vector_rep == "magnitude":
        min_val = np.min(vector_rep_val)
        max_val = np.max(vector_rep_val) + 0.00000001
        if degrees:
            min_map_val, max_map_val = -180, 180+0.001
        elif not degrees:
            min_map_val, max_map_val = -np.pi, np.pi
         
    norm = colors.Normalize(vmin=min_val, vmax=max_val)
    map_norm = colors.Normalize(vmin = min_map_val, vmax = max_map_val)

    if monitor_dpi is not None:
        fig, ax = plt.subplots(figsize=[image.shape[1] / monitor_dpi,
                                        image.shape[0] / monitor_dpi])
    else:
        fig, ax = plt.subplots()
    ax.set_title(title, loc='left', fontsize=20)

    # plot_style options
    if plot_style == "vector":
        Q = ax.quiver(
            x, y, u, v, units=quiver_units, color=color, pivot=pivot,
            angles=angles, scale_units=scale_units, scale=scale,
            headwidth=headwidth, headlength=headlength, minshaft=minshaft,
            headaxislength=headaxislength, width=width, minlength=minlength)
        length = np.max(np.hypot(u, v))
        ax.quiverkey(Q, 0.8, 1.025, length,
                     label='{:.2E} {}'.format(Decimal(length), units),
                     labelpos='E', coordinates='axes')
        
    elif plot_style == "colormap":

        if cmap is None:
            cmap = 'viridis'
        ax.quiver(
            x, y, u, v, vector_rep_val, color=color, cmap=cmap, norm=norm,
            units=quiver_units, pivot=pivot, angles=angles,
            scale_units=scale_units, scale=scale, headwidth=headwidth,
            alpha=alpha, headlength=headlength, headaxislength=headaxislength,
            width=width, minshaft=minshaft, minlength=minlength)

        # norm = colors.Normalize(vmin=min_val, vmax=max_val)
        # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        # sm.set_array([])
        # cbar = plt.colorbar(mappable=sm, fraction=0.046, pad=0.04,
        #                     drawedges=False)
        # cbar.set_ticks(ticks)
        # cbar.ax.set_ylabel(vector_label)

    elif plot_style == "colorwheel":

        if vector_rep != "angle":
            raise ValueError("`vector_rep`` must be set to 'angle' when "
                             "`plot_style`` is set to 'colorwheel'.")
        if cmap is None:
            cmap = cc.cm.colorwheel

        Q = ax.quiver(
            x, y, u, v, vector_rep_val, cmap=cmap, norm=norm, alpha=alpha,
            pivot=pivot, angles=angles, scale_units=scale_units,
            scale=scale, headwidth=headwidth, headlength=headlength,
            headaxislength=headaxislength, units=quiver_units, width=width,
            minshaft=minshaft, minlength=minlength)
        
    elif plot_style == "contour":

        if cmap is None:
            cmap = 'viridis'

        if isinstance(levels, list):
            levels_list = levels
        elif isinstance(levels, int):
            if vector_rep == "angle":
                levels_list = np.linspace(min_val, max_val, levels)
            elif vector_rep == "magnitude":
                levels_list = np.linspace(min_val, max_val, levels)

        plt.tricontourf(
            x, y, vector_rep_val, cmap=cmap, norm=norm, alpha=alpha,
            antialiased=antialiased, levels=levels_list)

        if not remove_vectors:
            ax.quiver(
                x, y, u, v, color=color, pivot=pivot, units=quiver_units,
                angles=angles, scale_units=scale_units,
                scale=scale, headwidth=headwidth, width=width,
                headlength=headlength, headaxislength=headaxislength,
                minshaft=minshaft, minlength=minlength)

        # cbar = plt.colorbar(mappable=contour_map, fraction=0.046, pad=0.04,
        #                     drawedges=False)
        # cbar.ax.tick_params(labelsize=14)
        # cbar.set_ticks(ticks)
        # cbar.ax.set_ylabel(vector_label, fontsize=14)

    elif plot_style == "polar_colorwheel":

        ax.quiver(
            x, y, u, v, color=color_list, pivot=pivot, units=quiver_units,
            angles=angles, scale_units=scale_units, scale=scale,
            headwidth=headwidth, width=width, headlength=headlength,
            headaxislength=headaxislength, minshaft=minshaft,
            minlength=minlength)
        
    elif plot_style == "contour_vector_map":
         if cmap is None:
            cmap = 'viridis'
         if cmap2 is not None: 
             cmap2 = 'cet_colorwheel'

         if isinstance(levels, list):
            levels_list = levels
         elif isinstance(levels, int):
            if vector_rep == "angle":
                levels_list = np.linspace(min_map_val, max_map_val, levels)
            elif vector_rep == "magnitude":
                levels_list = np.linspace(min_map_val, max_map_val, levels)


         plt.tricontourf(
            x, y, map_rep_val, cmap=cmap2, norm=map_norm, alpha=alpha,
            antialiased=antialiased, levels=levels_list)
         
         if not remove_vectors:
            ax.quiver(
            x, y, u, v, vector_rep_val, color=color, cmap=cmap, norm=norm,
            units=quiver_units, pivot=pivot, angles=angles,
            scale_units=scale_units, scale=scale, headwidth=headwidth,
            alpha=alpha, headlength=headlength, headaxislength=headaxislength,
            width=width, minshaft=minshaft, minlength=minlength)

    
    else:
        raise NameError("The plot_style you have chosen is not available.")
    
    if invert_y_axis:
        ax.set(aspect='equal')
        ax.set_xlim(0, image.shape[1])
        ax.set_ylim(image.shape[0], 0)

    if overlay:
        plt.imshow(image, cmap=image_cmap)

    if no_axis_info:
        plt.gca().axes.get_xaxis().set_visible(False)
        plt.gca().axes.get_yaxis().set_visible(False)

    if scalebar is True:
        scbar = ScaleBar(sampling, units, location="lower left", box_alpha=0.0,
                         color="white", scale_loc="top")
        plt.gca().add_artist(scbar)
    elif isinstance(scalebar, dict):
        scbar = ScaleBar(**scalebar)
        plt.gca().add_artist(scbar)

    # colorbars
    if (plot_style == "colormap" or plot_style == "colorwheel" or
            plot_style == "contour" or plot_style == "contour_vector_map"):

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(mappable=sm, fraction=0.046, pad=0.04,
                            drawedges=False, ax = ax)
        cbar.set_ticks(ticks)
        cbar.ax.set_ylabel(vector_label)
        if plot_style == "contour_vector_map":
            sm = plt.cm.ScalarMappable(cmap=cmap2, norm = map_norm)
            sm.set_array([])
            cbar2 = plt.colorbar(mappable=sm, fraction=0.046, pad=0.04,
                            drawedges=False, ax = ax)
            cbar2.ax.set_ylabel(map_label)
    elif plot_style == "polar_colorwheel":
        ax2 = fig.add_subplot(444)
        _make_color_wheel(ax2, rotation=None)
        ax2.set_axis_off()

    # plt.tight_layout()
    if isinstance(save, str):
        plt.savefig(fname=save + '_' + plot_style + '.png',
                    transparent=True, frameon=False, bbox_inches='tight',
                    pad_inches=None, dpi=300, labels=False)
    return ax
