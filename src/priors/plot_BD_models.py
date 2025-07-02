"""
Plotting Rebecca's BD model distributions.

Created: Tuesday 1st July 2025.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

# --- Plotting style ---
plt.rcParams.update({
    'axes.linewidth': 2.5,
    'font.size': 15,
    'figure.dpi': 100
})

# --- Switches ---
field_name = 'XMM'   # Choose from: 'COSMOS', 'XMM', 'CDFS'
bd_type = 'T'        # Choose from: 'M', 'L', 'T'

# --- File paths ---
cat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'
catalogues = {
    'COSMOS': cat_dir / 'COSMOS_dwarfCounts_300.0kpc.ascii',
    'XMM': cat_dir / 'SXDS_dwarfCounts_300.0kpc.ascii',
    'CDFS': cat_dir / 'EDFF_dwarfCounts_300.0kpc.ascii'
}

# --- Read and process tables ---
tables = {name: Table.read(path, format='ascii') for name, path in catalogues.items()}

# Read actual brown dward numbers
t_bd = Table.read('brown_dwarfs.fits')

# --- Select band names based on BD type ---
band_dict = {
    'M': [f'M{i}-{i+1}' for i in range(3, 9)],
    'L': [f'L{i}-{i+1}' for i in range(9)],
    'T': [f'T{i}-{i+1}' for i in range(8)]
}

# Check bd_type validity
if bd_type not in band_dict:
    raise ValueError("bd_type must be one of 'M', 'L', or 'T'")

selected_bands = band_dict[bd_type]

# --- Sum columns for each BD type ---
for tbl in tables.values():
    tbl[f'{bd_type}_total'] = sum(tbl[band] for band in selected_bands)

# --- Filter real sources by type ---
type_filters = {
    'M': [('M' in str(s)) for s in t_bd['stellar_type']],
    'L': [('L' in str(s)) for s in t_bd['stellar_type']],
    'T': [('T' in str(s)) for s in t_bd['stellar_type']],
}
t_bd_filtered = t_bd[type_filters[bd_type]]

# --- Area and binning ---
area_deg2 = 4.33
mapp = np.array(tables[field_name]['mapp'])

# --- Histogram of real sources ---
counts, bins = np.histogram(t_bd_filtered['mJ'], bins=mapp)
counts_per_deg2 = counts / area_deg2
bin_centers = 0.5 * (bins[1:] + bins[:-1])

# --- Plotting ---
plt.figure(figsize=(12, 8))

# Plot real sources
plt.step(bin_centers, counts_per_deg2, where='mid', color='black', lw=2, label=f'{bd_type} dwarfs (data)')

# Plot BD model for selected field
plt.plot(
    tables[field_name]['mapp'],
    tables[field_name][f'{bd_type}_total'],
    label=f'{field_name} model',
    color='tab:purple',
    marker='o'
)

plt.xlabel(r'$m_{\rm{AB}}$')
plt.ylabel('Number per square degree')
plt.title(f'{bd_type} dwarfs in {field_name}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()