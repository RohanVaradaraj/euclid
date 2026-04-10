"""
Check the coverage of the Euclid iDR1 tiles against the UltraVISTA pointing.
Adapted from fornax_footprints.py

Created: Wednesday 5th November 2025
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from pathlib import Path
import glob
import os
from matplotlib.patches import Circle
from astropy.coordinates import SkyCoord
import astropy.units as u
from shapely.geometry import Polygon
import pickle
import numpy as np

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

# === Directories ===
fornax_dir = Path.home() / 'euclid' / 'Y' / 'COSMOS'
vista_tile = Path.home().parents[1] / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS' / 'UVISTA_YJHK_DR6.fits'  

plt.figure(figsize=(10,10))
ax = plt.gca()

# === Euclid tiles ===
pattern = os.path.join(fornax_dir, f'*BGSUB-MOSAIC-NIR-Y*')
euclid_files = glob.glob(pattern)
tile_names = [file_name.split('-')[3].split('Y_TILE')[1] for file_name in euclid_files]

footprints = []
for i, file in enumerate(euclid_files):
    print(f"file {i+1} of {len(euclid_files)} → {os.path.basename(file)}")
    with fits.open(file, ignore_missing_simple=True) as hdu:
        wcs = WCS(hdu[0].header)
    footprints.append(wcs.calc_footprint())

# === Function to plot footprints ===
def plot_footprints(footprints, color='blue', lw=2, alpha=0.8, tilename=False):
    for i, fp in enumerate(footprints):
        fp = fp.tolist()
        fp.append(fp[0])  # close polygon
        ra, dec = zip(*fp)
        plt.plot(ra, dec, '-', color=color, lw=lw, alpha=alpha)
        if tilename:
            plt.text(np.mean(ra), np.mean(dec), tile_names[i],
                     fontsize=6, ha='center', va='center', color=color)
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.title("Euclid vs. UltraVISTA Footprints")

# === Plot Euclid ===
plot_footprints(footprints, color='green', lw=2, alpha=0.6, tilename=True)

# === Plot single UltraVISTA tile ===
with fits.open(vista_tile) as hdu:
    wcs_vista = WCS(hdu[0].header)
    vista_footprint = wcs_vista.calc_footprint()

plot_footprints([vista_footprint], color='black', lw=3, alpha=1)

# === Find overlaps ===
vista_poly = Polygon(vista_footprint)
euclid_polygons = [Polygon(fp) for fp in footprints]

euclid_within_vista = [
    tile_names[i]
    for i, ep in enumerate(euclid_polygons)
    if ep.intersects(vista_poly)
]

print("\nEuclid tiles overlapping with VISTA:")
if euclid_within_vista:
    for t in euclid_within_vista:
        print(f"  - {t}")
else:
    print("  None")

# === Highlight overlapping Euclid tiles ===
for t in euclid_within_vista:
    idx = tile_names.index(t)
    fp = footprints[idx].tolist()
    fp.append(fp[0])
    ra, dec = zip(*fp)
    plt.plot(ra, dec, '-', color='red', lw=2.5, alpha=0.9)

# === Final plot ===
plt.plot([], [], color='black', lw=3, label='UltraVISTA')
plt.plot([], [], color='green', lw=2, label='Euclid DR1')
plt.plot([], [], color='red', lw=2.5, label='Overlap')
plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()
plt.show()

# === Save results ===
with open('euclid_within_ultravista.pkl', 'wb') as f:
    pickle.dump(euclid_within_vista, f)
