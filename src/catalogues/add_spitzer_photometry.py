#!/usr/bin/env python3

"""
Crossmatch my VISTA/Euclid-selected catalogues with the big GIKs catalogues to get the IRAC photometry.
Running custom photometry seems to be slow, so this is a quicker way of recovering it.
A quick TOPCAT crossmatch reveals ~50k sources are not recovered. For these, I will use source extractor python.

Created: Wednesday 16th October 2024.
"""

import numpy as np
from astropy.table import Table, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u 
from astropy.table import hstack
from astropy.table import vstack
from astropy.table import unique
from astropy.table import join
import os
from pathlib import Path
import sep

def grid_depths(gridTable, x, y, faster = True, verbose = False, nearby = False):
    
   ''' Code to find the closest depth measurement from my previous analysis. Faster than truly local depths '''
   
   xgrid = gridTable['x']
   ygrid = gridTable['y']
   keys = gridTable.colnames
   
   depthsOverField = gridTable['depths']
   
   ## Make an output array
   depthArray = np.zeros(x.size)
   depthArray[:] = -99.0
   
   if faster:
       if verbose:
           print("Using faster method.")
           print("Input array size is ", x.size)
       deltay = np.min(ygrid)
       deltax = np.min(xgrid)
      
       for xi in range(xgrid.size):
           
           xmin = xgrid[xi] - deltax
           xmax = xgrid[xi] + deltax
           ymin = ygrid[xi] - deltay
           ymax = ygrid[xi] + deltay

           ii = (x > xmin) & (x <= xmax) & (y > ymin) & (y <= ymax)
           
           depthArray[ii] = depthsOverField[xi]
               
   else:
       
        ## Find the closest point to the objects x and y positions
        ## Loop!
       for xi in range(x.size):
           
       ## make a radius array
           deltax = (xgrid - x[xi])
           deltay = (ygrid - y[xi])
           radius = np.sqrt(deltax*deltax + deltay*deltay)
           mini = np.argmin(radius)
           
       ## try using argpartition
           numpoints = 10
           idx = np.argpartition(radius, numpoints)
           
           if nearby:
               
               mini = idx[0:numpoints]
               print("The nearby depths are = ", depthsOverField[mini])
               print("Before = ", depthsOverField[mini][0])      

           depthArray[xi] = depthsOverField[mini][0]

   
   return depthArray

# Read in the VISTA/Euclid-selected catalogue, the one we need to add the IRAC photometry to!
vista_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS' / 'other'
#vista_file = 'COSMOSFULL_DR3_UNMASKED_Ks_2024_11_06_2.0as_IRAC_2.8as_ALL.fits'
vista_file = 'COSMOSFULL_DR3_MASKVISTADET_HSC-Z_DR3_2025_06_05_1.8as_IRAC_2.8as_ALL.fits' # HSC-Z selected
vista_cat = Table.read(vista_dir / vista_file)

#! LAE
# vista_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
# vista_file = 'LAE.fits'
# vista_cat = Table.read(vista_dir / vista_file)

print(len(vista_cat))
print(vista_cat.colnames)

# Read in the Ks catalogue
cat_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS'
k_cat = Table.read(cat_dir / 'COSMOSFULL_DR3_UNMASKED_Ks_2023_03_11_2.0as_IRAC_2.8as_ALL.fits')

print(k_cat.colnames)

# Crossmatch the VISTA/Euclid-selected catalogue with the GIKs catalogue to get the IRAC photometry
coord1 = SkyCoord(ra=vista_cat['RA'], dec=vista_cat['DEC'])
coord2 = SkyCoord(ra=k_cat['RA'], dec=k_cat['DEC'])

idx_vista, idx_k, d2d, d3d = coord1.search_around_sky(coord2, 0.5*u.arcsec)

# Ensure indices are within bounds
valid_idx = (idx_vista < len(vista_cat)) & (idx_k < len(k_cat))
idx_vista = idx_vista[valid_idx]
idx_k = idx_k[valid_idx]
d2d = d2d[valid_idx]

print(len(idx_vista))
print(len(idx_k))

# Print the number of matches
print('Fraction of matches:', len(idx_vista)/len(vista_cat))
print('Number of matches:', len(idx_vista))
print('Length of VISTA catalogue:', len(vista_cat))

# Add the flux_ch1cds, flux_ch2cds, err_ch1cds, err_ch2cds columns to the VISTA catalogue
vista_cat['flux_ch1cds'] = np.nan
vista_cat['flux_ch2cds'] = np.nan
vista_cat['err_ch1cds'] = np.nan
vista_cat['err_ch2cds'] = np.nan

# Add a separation column to the VISTA catalogue
vista_cat['separation'] = np.nan

# Add the ra dec of the crossmatched sources to the VISTA catalogue
vista_cat['ra_irac'] = np.nan
vista_cat['dec_irac'] = np.nan

print(idx_vista)
print(idx_k)

# Add the IRAC photometry to the VISTA catalogue
vista_cat['flux_ch1cds'][idx_vista] = k_cat['flux_ch1cds'][idx_k]
vista_cat['flux_ch2cds'][idx_vista] = k_cat['flux_ch2cds'][idx_k]
vista_cat['err_ch1cds'][idx_vista] = k_cat['err_ch1cds'][idx_k]
vista_cat['err_ch2cds'][idx_vista] = k_cat['err_ch2cds'][idx_k]

# Add the separation column to the VISTA catalogue
vista_cat['separation'][idx_vista] = d2d.arcsec

# Add the ra dec of the crossmatched sources to the VISTA catalogue
vista_cat['ra_irac'][idx_vista] = k_cat['RA'][idx_k]
vista_cat['dec_irac'][idx_vista] = k_cat['DEC'][idx_k]

# Now get the RA,DEC where the IRAC photometry is missing
vista_cat_miss = vista_cat[np.isnan(vista_cat['flux_ch1cds'])]
x_missing, y_missing = vista_cat_miss['X_IMAGE'], vista_cat_miss['Y_IMAGE']

print(f'There are {len(vista_cat_miss)} sources missing IRAC photometry.')

# First open ch1 image
image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'
with fits.open(image_dir / 'COSMOS_ch1_COMPLETE_microJy.fits') as hdu:

    data = hdu[0].data
    data = data.byteswap().newbyteorder()
    header = hdu[0].header
    print('ch1 image opened')

    flux, fluxerr, flag = sep.sum_circle(data, x_missing, y_missing, 1.4/0.15, subpix = 5)
    flux = np.array(flux)
    flux = 10 ** (-(48.6+23.9)/2.5) * flux
    flux = flux/0.585
    print('flux extracted')

# Replace nans with flux values
print(flux)
vista_cat_miss['flux_ch1cds'] = flux


# now do ch2 image
with fits.open(image_dir / 'COSMOS_ch2_COMPLETE_microJy.fits') as hdu:

    data = hdu[0].data
    data = data.byteswap().newbyteorder()
    header = hdu[0].header
    print('ch2 image opened')

    flux, fluxerr, flag = sep.sum_circle(data, x_missing, y_missing, 1.4/0.15, subpix = 5)
    flux = np.array(flux)
    flux = 10 ** (-(48.6+23.9)/2.5) * flux
    flux = flux/0.569
    print('flux extracted')

# Replace nans with flux values
print(flux)
vista_cat_miss['flux_ch2cds'] = flux

# Add the missing IRAC photometry to the VISTA catalogue
vista_cat['flux_ch1cds'][np.isnan(vista_cat['flux_ch1cds'])] = vista_cat_miss['flux_ch1cds']
vista_cat['flux_ch2cds'][np.isnan(vista_cat['flux_ch2cds'])] = vista_cat_miss['flux_ch2cds']
print('Added missing IRAC photometry')

# Get errors from the depth tables
depth_dir = Path.cwd().parents[3] / 'data' / 'depths' / 'COSMOS' / 'phot'

print('Getting errors from depth tables')

ch1_depths = Table.read(depth_dir / 'ch1cds_2.8as_gridDepths_300_200.fits')
ch2_depths = Table.read(depth_dir / 'ch2cds_2.8as_gridDepths_300_200.fits')

depths_ch1 = grid_depths(ch1_depths, x_missing, y_missing)
depths_ch2 = grid_depths(ch2_depths, x_missing, y_missing)

print('ch1 depths:', depths_ch1)
print('ch2 depths:', depths_ch2)

error_ch1 = 0.2 * 10 ** (-(48.6+depths_ch1)/2.5)
error_ch2 = 0.2 * 10 ** (-(48.6+depths_ch2)/2.5)

# Replace nans with error values
vista_cat['err_ch1cds'][np.isnan(vista_cat['err_ch1cds'])] = error_ch1
vista_cat['err_ch2cds'][np.isnan(vista_cat['err_ch2cds'])] = error_ch2

print('Added missing IRAC errors')

# Save the VISTA catalogue with the IRAC photometry
vista_cat.write(vista_dir /vista_file, format='fits', overwrite=True)






