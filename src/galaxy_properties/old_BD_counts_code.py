#!/usr/bin/env python3

"""

Compute expected number of brown dwarfs in our survey based on work by Rebecca's masters student.

Created: Wednesday 24th May 2023

"""

import numpy as np
from pathlib import Path
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM

# Directories
codeDir = Path.cwd()
catDir = codeDir.parent.parent / 'ref_catalogues'
catName = 'CDFS_dwarfCounts_v1.ascii'

########### SAMPLE SET-UP ##############

# Choose the field
field = 'XMM'

# Define the survey area here
if field == 'CDFS':
    survey_area = 3.89 # sq. deg. We expect CDFS and XMM to be similar based on their position on the sky. Scale height of disk assumed to be 300pc.
else:
    survey_area = 4.33

# Define absolute UV magnitude range of interest.

# CDFS
if field == 'CDFS':
    Muv_range = np.array([-23.5, -22.4]) # Range of objects
    Muv_range = np.array([-23.5, -22.4]) # Down to depth
# XMM
else:
    Muv_range = np.array([-22.9, -21.6]) # Range of objects
    Muv_range = np.array([-23.5, -21.8]) # Down to depth

# Deal with covering fraction
survey_area *= 0.83

# Mean redshift of the sample
z_mean = 6.74 # XMM = 6.73, CDFS = 6.75.

# Define the cosmology
H = 70
omegaM = 0.3
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)

# And then get the luminosity distance
DL = cosmo.luminosity_distance(z_mean) * 10 ** 6 # Convert to pc

########### TABLE SET-UP ##############

# Open the brown dwarf number density table
t = ascii.read(catDir / catName)

# Get the column names
column_names = t.colnames

column_names = ['mapp', 'T7-8', 'L6-7', 'L7-8'] # Most frequent
column_names = ['mapp', 'T1-2', 'T2-3', 'T7-8', 'T3-4', 'M5-6', 'M6-7', 'M7-8', 'L6-7', 'L7-8']
column_names = ['mapp', 'M3-4', 'M4-5', ]

t = t[column_names]
print(t)
# We are primarily interested in L and T dwarfs, so remove M types.
substring = 'M'
columns_to_remove = [col for col in column_names if substring in col]

# Remove M-type columns
#or col in columns_to_remove:
#    t.remove_column(col)

#print(t.colnames)


########### BEGIN CALCULATION ###########

# Convert my absolute UV magnitude range into apparent magnitudes.
m_app = Muv_range + 5 * np.log10(DL.value / 10)  - 2.5 * np.log10(1 + z_mean)

# Find values within this range
idx = np.where(np.logical_and(t['mapp'] >= m_app[0], t['mapp'] <= m_app[1]))

# Restrict the table. Remove the 'mapp' so we can simply sum over remaining table elements.
t = t[:][idx]

# Exclude 'mapp' column from column names
columns_to_sum = [col for col in t.colnames if col != 'mapp']

# Convert the desired columns to a NumPy array
data_array = np.array([t[col] for col in columns_to_sum])

# Sum the rows of the array
row_sums = np.sum(data_array, axis=0)

# Calculate the sum of the row sums
total_sum = np.sum(row_sums)

# Multiply by the area: this takes number of BDs per square degree to number of BDs in our survey
N_BD = survey_area * total_sum 

print(f'Number of brown dwarfs expected in {field}: {int(N_BD)}.')










