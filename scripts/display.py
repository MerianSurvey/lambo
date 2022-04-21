import os
import copy

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy import wcs
from astropy.convolution import convolve
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clip, sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import (AsymmetricPercentileInterval,
                                   ZScaleInterval, make_lupton_rgb)
from matplotlib import colors, rcParams
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.ticker import (AutoMinorLocator, FormatStrFormatter,
                               MaxNLocator, NullFormatter)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from palettable.colorbrewer.sequential import (Blues_9, Greys_9, OrRd_9,
                                               Purples_9, YlGn_9)


def display_merian_cutout_rgb(images, filters=list('griz') + ['N708'],
                              ax=None, half_width=None,
                              minimum=-0.15,
                              channel_map=None,
                              N708_strength=0.7,
                              stretch=1.2, Q=3,
                              color_norm=None
                              ):
    import scarlet
    from scarlet.display import AsinhMapping

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))

    # Crop
    if half_width is not None:
        if half_width * 2 < min(images.shape[1], images.shape[2]):
            cen = (images.shape[1] // 2, images.shape[2] // 2)
            images = images[:, cen[0] - half_width:cen[0] +
                            half_width, cen[1] - half_width:cen[1] + half_width]

    # Norm color
    if color_norm is None:
        color_norm = {'g': 1.9, 'r': 1.2, 'i': 1.0,
                      'z': 0.85, 'y': 0.5, 'N708': 1.2, 'N540': 1.0}

    if channel_map is None:
        channel_map = scarlet.display.channels_to_rgb(len('griz'))
        _map = np.zeros((3, len(filters)))
        _map[:, :4] = channel_map
        if 'N708' in filters:
            _map[0, 4] = N708_strength  # N708 for red
        # _map[0, 5] = 0.2
        _map /= _map.sum(axis=1)[:, None]
    else:
        _map = channel_map

    f_c = np.array([color_norm[filt] for filt in filters])
    _images = images * f_c[:, np.newaxis, np.newaxis]

    # Display
    norm = AsinhMapping(minimum=minimum, stretch=stretch, Q=Q)

    img_rgb = scarlet.display.img_to_rgb(
        _images, norm=norm, channel_map=_map)

    ax.imshow(img_rgb, origin='lower')
    ax.axis('off')

    if ax is None:
        return fig
    return ax, img_rgb
