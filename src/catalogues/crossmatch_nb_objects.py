"""
Crossmatch the HSC narrowband catalogues.

Created: Wendesday 24th July 2024.
"""

from astropy.io import ascii
from astropy.table import Table
import numpy as np
from pathlib import Path

hsc_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'
hsc_cat ='table4.dat'

#hsc_cat = ascii.read(hsc_dir / hsc_cat, format = 'basic')

# Define the column names and types
column_names = ['HSC', 'RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs', 'zspec', 'g', 'r', 'i', 'z', 'y', 'nb', 'which_nb', 'feature', 'ref']
column_formats = ['S20', 'i4', 'i4', 'f4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'S20', 'S20', 'S50']

def ra_to_degrees(rah, ram, ras):
    """Convert RA from hours, minutes, seconds to degrees."""
    return 15.0 * (rah + ram / 60.0 + ras / 3600.0)

def dec_to_degrees(decd, decm, decs, sign):
    """Convert DEC from degrees, minutes, seconds to degrees."""
    dec = abs(decd) + decm / 60.0 + decs / 3600.0
    return -dec if sign == '-' else dec

# Read the data from the .dat file
data = []
with open(str(hsc_dir/hsc_cat), 'r') as f:
    for line in f:
        # Split the line into components
        parts = line.strip().split()
        
        # Extract the relevant parts of the line
        hsc_id = parts[0] + ' ' + parts[1]
        rah = int(parts[2])
        ram = int(parts[3])
        ras = float(parts[4])
        dec_sign = parts[5][0]
        ded = int(parts[5][1:])
        dem = int(parts[6])
        des = float(parts[7])
        if dec_sign == '-':
            des = -des
        zspec = float(parts[8])
        g = float(parts[9])
        r = float(parts[10])
        i = float(parts[11])
        z = float(parts[12])
        y = float(parts[13])
        nb = float(parts[14])
        nb_type = parts[15]
        feature = parts[16] if parts[16] != '' else 'None'
        ref = ' '.join(parts[17:]) if len(parts) > 17 else 'None'
        
        # Convert RA and DEC to degrees
        ra_deg = ra_to_degrees(rah, ram, ras)
        dec_deg = dec_to_degrees(ded, dem, des, dec_sign)
        
        # Append the extracted data to the list, including RA and DEC in degrees
        data.append((hsc_id, ra_deg, dec_deg, zspec, g, r, i, z, y, nb, nb_type, feature, ref))

# Define the new column names and formats including RA and DEC in degrees
column_names = ['HSC', 'RA_deg', 'DEC_deg', 'zspec', 'g', 'r', 'i', 'z', 'y', 'nb', 'nb_type', 'feature', 'ref']
column_formats = ['S20', 'f8', 'f8', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'S6', 'S20', 'S50']

# Convert the list to a structured NumPy array
dtype = [(name, fmt) for name, fmt in zip(column_names, column_formats)]
structured_data = np.array(data, dtype=dtype)

# Create an Astropy Table from the structured array
table = Table(structured_data)

table.write(hsc_dir / 'hsc_nb_catalogue.fits', format='fits', overwrite=True)