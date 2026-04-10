"""
Make nice JWST plots of the galaxies!

Created: Wednesday 8th October 2025.
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import glob
from shapely.geometry import Point, Polygon
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from astropy.cosmology import Planck18 as cosmo
from matplotlib.patches import Rectangle
from astropy.wcs import WCS


# Cycling dithering pattern with small size.
position_index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
x_offset = [-0.2328, 0.1292, -0.0259, -0.0776, 0.0000, -0.1552, 0.0000, 0.0517, -0.1552, -0.0629, 0.0581, -0.0157]
y_offset = [-0.0774, 0.1855, -0.1333, 0.2415, 0.0000, 0.0553, 0.1604, 0.1063, -0.1585, -0.0010, 0.1156, -0.0753]




def findPlotLimits(data: np.ndarray) -> tuple:

    mean = np.mean(data)
    std_dev = np.std(data)

    # Sigma clip the data
    data = data[(data < mean + 4 * std_dev)]

    # Recalculate mean and std_dev
    mean = np.mean(data)
    std_dev = np.std(data)

    lower = mean - 2 * std_dev
    upper = mean + 6 * std_dev

    return lower, upper


cutout_dir = Path.cwd() / 'cutouts'

# Get all files in the cutout directory
files = sorted(glob.glob(str(cutout_dir / '*f444w.fits')))

# Remove 615679
files = [f for f in files if '615679' not in f]
# Remove 178396
#files = [f for f in files if '178396' not in f]

# Replace the 178396_f444w.fits with 178396_f115w.fits
#files = [f.replace('178396_f444w.fits', '178396_f50w.fits') if '178396_f444w.fits' in f else f for f in files]

# Replace the 178396 with my own JWST imaging

redshifts = {
    '387777': 6.7614,
    '878786': 7.0975,
    '178396': 6.7014,
}

# Set up subplot figure
fig, axes = plt.subplots(1, len(files), figsize=(4 * len(files), 4))
if len(files) == 1:
    axes = [axes]  # Ensure axes is iterable

for ax, file in zip(axes, files):

    with fits.open(file) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        pix_scale = np.abs(header['CDELT1']) * 3600
        wcs = WCS(header)

        # This file is 10x10 arcsec, cut down to 2x2 arcsec
        cut_size = int(6 / pix_scale / 2)
        center = np.array(data.shape) // 2
        #data = data[center[0]-cut_size:center[0]+cut_size, center[1]-cut_size:center[1]+cut_size]

        # Print RA DEC of the original center
        ra, dec = wcs.wcs_pix2world(center[1], center[0], 0)
        print(f"File: {file}, Original center RA: {ra}, DEC: {dec}")

        # ---- Apply manual center shift for a specific ID ----
        if "878786" in file:
            # shift the center DOWN (positive in y index)
            # adjust the number of pixels to taste (e.g., +10 = down by 10 pixels)
            center[0] -= 28   # move down by 10 pixels
            center[1] -= 11   # optionally shift right if needed

        if "178396" in file:
            # shift the center DOWN (positive in y index)
            # adjust the number of pixels to taste (e.g., +10 = down by 10 pixels)
            center[0] -= 10   # move down by 10 pixels
            center[1] -= 0   # optionally shift right if needed

        if "387777" in file:
            # shift the center DOWN (positive in y index)
            # adjust the number of pixels to taste (e.g., +10 = down by 10 pixels)
            center[0] += 0   # move down by 10 pixels
            center[1] += 6   # optionally shift right if needed
        # ------------------------------------------------------

        # For each file, now print the RA DEC of the new center

        ra, dec = wcs.wcs_pix2world(center[1], center[0], 0)
        print(f"File: {file}, New center RA: {ra}, DEC: {dec}")
        continue

        data = data[
            center[0]-cut_size:center[0]+cut_size,
            center[1]-cut_size:center[1]+cut_size
        ]

    print(center, data.shape, pix_scale)
    vmin, vmax = findPlotLimits(data)

    # Display the image
    ax.imshow(data, origin='lower', cmap='gist_yarg', vmin=vmin, vmax=vmax)

    # Draw 3x3 arcsec box on centre of image
    box_size = int(3 / pix_scale / 2)
    rect = Rectangle((data.shape[1]//2 - box_size, data.shape[0]//2 - box_size),
                     2*box_size, 2*box_size,
                     linewidth=1, edgecolor='limegreen', facecolor='none', lw=2, alpha=0.4)
    ax.add_patch(rect)

    # Now loop through the dither positions and plot the boxes
    for dx, dy in zip(x_offset, y_offset):
        # Convert offsets from arcsec to pixels
        dx_pix = dx / pix_scale
        dy_pix = dy / pix_scale

        # Calculate the new center position
        new_center_x = data.shape[1]//2 + dx_pix
        new_center_y = data.shape[0]//2 + dy_pix

        # Create a rectangle at the new position
        rect = Rectangle((new_center_x - box_size, new_center_y - box_size),
                         2*box_size, 2*box_size,
                         linewidth=0.5, edgecolor='limegreen', facecolor='none', alpha=0.4, lw=2)
        ax.add_patch(rect)


    # --- Add a scale bar ---   
    #scalebar_length_arcsec = 1  # choose the physical length to show, e.g. 0.5"
    #scalebar_length_pix = scalebar_length_arcsec / pix_scale

    # Replace: scalebar_length_arcsec = 1
    z = redshifts.get(Path(file).stem.split('_')[0], 7.0)  # Default to z=7.0 if not found
    arcsec_per_kpc = cosmo.arcsec_per_kpc_proper(z).value
    scalebar_length_arcsec = 2 * arcsec_per_kpc  # arcsec corresponding to 1 kpc
    scalebar_length_pix = scalebar_length_arcsec / pix_scale


    fontprops = fm.FontProperties(size=15)
    scalebar = AnchoredSizeBar(
        ax.transData,
        scalebar_length_pix,
        "2 kpc",
        'lower right',
        pad=0.3,
        color='black',
        frameon=False,
        size_vertical=1,
        fontproperties=fontprops
    )
    ax.add_artist(scalebar)


    # # Add another scale bar at bottom left with 0.3 arcsec
    # scalebar_length_arcsec = 0.3
    # scalebar_length_pix = scalebar_length_arcsec / pix_scale
    # scalebar = AnchoredSizeBar(
    #     ax.transData,
    #     scalebar_length_pix,
    #     "0.3\"",
    #     'lower left',
    #     pad=0.3,
    #     color='black',
    #     frameon=False,
    #     size_vertical=1,
    #     fontproperties=fontprops
    # )
    # ax.add_artist(scalebar)

    ax.axis('off')

plt.tight_layout()
plt.savefig('JWST_stamps.pdf', dpi=100, bbox_inches='tight')

plt.show()



