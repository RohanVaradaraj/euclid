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
files = [f for f in files if '178396' not in f]

# Replace the 178396_f444w.fits with 178396_f115w.fits
#files = [f.replace('178396_f444w.fits', '178396_f50w.fits') if '178396_f444w.fits' in f else f for f in files]

# Replace the 178396 with my own JWST imaging

redshifts = {
    '387777': 6.7614,
    '878786': 7.0975,
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

        # This file is 10x10 arcsec, cut down to 2x2 arcsec
        cut_size = int(4 / pix_scale / 2)
        center = np.array(data.shape) // 2
        #data = data[center[0]-cut_size:center[0]+cut_size, center[1]-cut_size:center[1]+cut_size]

        # ---- Apply manual center shift for a specific ID ----
        if "878786" in file:
            # shift the center DOWN (positive in y index)
            # adjust the number of pixels to taste (e.g., +10 = down by 10 pixels)
            center[0] -= 25   # move down by 10 pixels
            center[1] -= 5   # optionally shift right if needed
        # ------------------------------------------------------

        data = data[
            center[0]-cut_size:center[0]+cut_size,
            center[1]-cut_size:center[1]+cut_size
        ]

    vmin, vmax = findPlotLimits(data)

    # Display the image
    ax.imshow(data, origin='lower', cmap='gist_yarg', vmin=vmin, vmax=vmax)

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
    ax.axis('off')

    # Draw 0.3x0.3 arcsec box  
    box_size = int(3 / pix_scale / 2)
    box = Polygon([ (center[1]-box_size, center[0]-box_size),
                    (center[1]+box_size, center[0]-box_size),
                    (center[1]+box_size, center[0]+box_size),
                    (center[1]-box_size, center[0]+box_size) ])
    x, y = box.exterior.xy
    #ax.plot(x - (center[1]-cut_size), y - (center[0]-cut_size), color='red', lw=1.5)

    # Extract ID from filename
    id_str = Path(file).stem.split('_')[0]
    #ax.set_title(f'ID {id_str}', fontsize=14)


plt.tight_layout()
plt.savefig('JWST_stamps.pdf', dpi=100, bbox_inches='tight')



plt.show()



