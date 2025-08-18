from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
import sep
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Circle

import sys
cutout_path = Path.cwd().parents[0] / 'cutouts'
sys.path.append(str(cutout_path))
from cutout_codes import findPlotLimits

# Load catalogue
t = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates' /
               'Euclid_UltraVISTA_z7_sample.fits')
t.sort('Muv')
#print(t)

# Test with one candidate
# i = 1
# t = t[i-1:i]
# t = t[1:2]

# Test with specific ID
ID_needed = 825188
t = t[t['ID'] == ID_needed]
print(t)

ra = t['RA']
dec = t['DEC']
Muv = t['Muv']

# Load the images
image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'
images = {
    'Y': image_dir / 'UVISTA_Y_DR6.fits',
    'J': image_dir / 'UVISTA_J_DR6.fits',
    'YJHK': image_dir / 'UVISTA_YJHK_DR6.fits',
    'H': image_dir / 'UVISTA_H_DR6.fits',
    'K': image_dir / 'UVISTA_K_DR6.fits',
    'YE': image_dir / 'Euclid_Y_vista_matched.fits',
    'JE': image_dir / 'Euclid_J_vista_matched.fits',
    'HE': image_dir / 'Euclid_H_vista_matched.fits',
}

psf_corrections = {
    'Y': 0.69,
    'J': 0.72,
    'YJHK': 0.70,  # Average of Y, J, H, K
    'H': 0.76,
    'K': 0.78,
    'YE': 0.67,
    'JE': 0.67,
    'HE': 0.68,
}


# Collect cutout results for plotting later
plot_data = []

for band, image_path in images.items():
    with fits.open(image_path) as hdul:
        image_data = hdul[0].data
        wcs = WCS(hdul[0].header)

        for j in range(len(ra)):
            # Pixel coordinates
            x, y = wcs.world_to_pixel_values(ra[j], dec[j])
            print(x,y)
            position = (x, y)
            cutout_size = (64, 64)
            cutout = Cutout2D(image_data, position, cutout_size, wcs=wcs)
            print(cutout.data.shape)
            vmin, vmax = findPlotLimits(cutout.data)

            # Make C-contiguous and subtract background
            data_c = np.ascontiguousarray(cutout.data, dtype=float)
            bkg = sep.Background(data_c)
            data_sub = data_c - bkg

            # # Object detection
            objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)

            # kronrad, krflag = sep.kron_radius(data_sub, x, y, a, b, theta, 6.0)
            # flux, fluxerr, flag = sep.sum_ellipse(data_sub, x, y, a, b, theta, 2.5*kronrad,
            #                                     subpix=1)

            # Store for later plotting
            plot_data.append((band, data_sub, objects))

# === Plot all bands in one figure ===
n_bands = len(plot_data)
fig, axes = plt.subplots(1, n_bands, figsize=(3*n_bands, 3), squeeze=False)

for ax, (band, data_sub, objects) in zip(axes[0], plot_data):
    m, s = np.mean(data_sub), np.std(data_sub)
    ax.imshow(data_sub, interpolation='nearest', cmap='gray',
            vmin=m-s, vmax=m+s, origin='lower')

    # Blue circle at centre: 1.8 arcsec diameter, 0.15 arcsec/pix
    circle_radius = 1.8 / 2 / 0.15
    circle = Circle(xy=(64 // 2, 64 // 2), radius=circle_radius,
                    facecolor='none', edgecolor='blue')
    ax.add_artist(circle)

    # Red ellipses for each object
    for k in range(len(objects)):
        e = Ellipse(xy=(objects['x'][k], objects['y'][k]),
                    width=6*objects['a'][k],
                    height=6*objects['b'][k],
                    angle=objects['theta'][k] * 180. / np.pi,
                    facecolor='none', edgecolor='red')
        ax.add_artist(e)

    ax.set_title(band)
    ax.axis('off')

plt.tight_layout()
plt.show()

