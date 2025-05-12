"""
Using Rebecca's BD models to predict BD counts in COSMOS

Created: Saturday 8th February 2025.
"""

from astropy.table import Table
from astropy.io import ascii
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

# Load the table
table_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'
t = ascii.read(table_dir / 'COSMOS_dwarfCounts_300.0kpc.ascii')

# Plot number counts as function of mag
# plt.plot(bd_model['mapp'], np.sum([bd_model[col] for col in bd_model.colnames[1:]], axis=0))
# plt.show()

# Muv_range = np.array([-22.6, -20.2])
# survey_area = 1.78 # sq. deg.

#survey_area *= 0.83 # Covering fraction

# Euclid
Muv_range = np.array([-22.2, -20.2])
survey_area = 0.65 # sq. deg.

#survey_area *= 0.83 # Covering fraction

# Mean redshift of sample
z_mean = 6.80

# Define the cosmology
H = 70
omegaM = 0.3
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)

# And then get the luminosity distance
DL = cosmo.luminosity_distance(z_mean) * 10 ** 6 # Convert to pc

# Get the column names
column_names = t.colnames

#column_names = ['mapp', 'T7-8', 'L6-7', 'L7-8'] # Most frequent
column_names = ['mapp', 'T1-2', 'T2-3', 'T7-8', 'T3-4', 'M5-6', 'M6-7', 'M7-8', 'L6-7', 'L7-8']
# column_names = ['mapp', 'M3-4', 'M4-5', ]

#t = t[column_names]
#print(t)
# We are primarily interested in L and T dwarfs, so remove M types.
# substring = 'M'
# columns_to_remove = [col for col in column_names if substring in col]

# Convert my absolute UV magnitude range into apparent magnitudes.
m_app = Muv_range + 5 * np.log10(DL.value / 10)  - 2.5 * np.log10(1 + z_mean)

# Find values within this range
idx = np.where(np.logical_and(t['mapp'] >= m_app[0], t['mapp'] <= m_app[1]))[0]
#print(idx)

# Restrict the table. Remove the 'mapp' so we can simply sum over remaining table elements.
t = t[idx]
#print(t)
columns_to_sum = [col for col in t.colnames if col != 'mapp']

# Convert the desired columns to a NumPy array
data_array = np.array([t[col] for col in columns_to_sum])

# Sum the rows of the array
row_sums = np.sum(data_array, axis=0)

# Calculate the sum of the row sums
total_sum = np.sum(row_sums)

# Multiply by the area: this takes number of BDs per square degree to number of BDs in our survey
N_BD = survey_area * total_sum 

print(f"Expected number of brown dwarfs in 1.72 deg² (depth 26.1): {N_BD:.2f}")

