#!/usr/bin/env python3

"""
An improved version of compute_zmax.py, which was adapted from my first paper code specifically for COSMOS Euclid.

Created: Tuesday 13th May 2025.
"""
import numpy as np
import os
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, Column
import matplotlib.pyplot as plt
from astropy.constants import c
from scipy.integrate import simps, romb
from pathlib import Path
from astropy.wcs import WCS
import glob
import sys
from scipy.interpolate import interp1d

sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file, filter_files

depth_path = Path.cwd().parents[0] / 'depths'
sys.path.append(str(depth_path))
from new_depth_codes import grid_depths

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100



def mag_to_flux(m):
	'''Convert mags to flux'''
	flux = 10**(-0.4*(m+48.6))
	return flux



def flux_to_mag(flux):
	'''Convert flux to mag'''
	mag = -2.5*np.log10(flux)-48.6
	return mag



# Define the cosmology
H = 70
omegaM = 0.3
omegaV = 0.7
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)

#! Field name
field_name = 'COSMOS'

#! Det/non-det filters
filters = {
    'HSC-Z_DR3': {'type': 'detection', 'value': 5},
    'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
}

#! Run type
run_type = ''

#! Redshifting parameters
dz = 0.01
zmin = 5.5
zmax = 6.5

#! Vmax parameters
if field_name == 'COSMOS':
    field_area = 1.7201431257141546 # COSMOS
if field_name == 'XMM':
    field_area = 4.335562874179884 # XMM
covering_fraction = 0.8
field_area *= covering_fraction

#! Photometry params
aperture_size = 2.0 # arcsec

#! SED Fitting folder
# Generate name of the directory we want to use to make the catalogue
det_filters = [f for f, t in filters.items() if t['type'] in ['detection', 'stacked-detection']]
det_filter_str = '_'.join(det_filters)
if run_type != '':
    folder = f'det_{det_filter_str}_{run_type}_z7'
else:
    folder = f'det_{det_filter_str}_z7'

#! Get the list of candidates
obj_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / 'best_fits'
obj_list = glob.glob(str(obj_dir / folder / '*.spec'))

# Get the IDs
IDs = [spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0] for spec_file in obj_list]
IDs = [int(ID) for ID in IDs]

#! Candidate catalogue
if field_name == 'COSMOS':
    if run_type == '':
        cat_name = 'COSMOS_5sig_det_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14.fits' # just vista
        cat_name = 'COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06.fits' # z=6
    if run_type == 'with_euclid':
        cat_name = 'COSMOS_5sig_det_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits' # with euclid

if field_name == 'XMM': 
    #cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_13.fits'
    cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits'

#! Read in catalogue of candidates
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
t = Table.read(cat_dir / cat_name)


#! Add a column for zmax if it doesnt exist already
if 'zmax' not in t.colnames:
    t['zmax'] = Column(np.zeros(len(t)), name='zmax', dtype=float)
    t['Vmax'] = Column(np.zeros(len(t)), name='Vmax', dtype=float)

#! Read in depths. Based on detection filter(s).
det_filter = [f for f, t in filters.items() if t['type'] in ['detection', 'stacked-detection']]
if len(det_filter) == 1:
    det_filter = det_filter[0]
#* TO DO: write separate clause if there are multiple detection filters. But in principle just make a stacked image.

depth_dir = Path.cwd().parents[3] / 'data' / 'depths'

depth_files = []

# If in the VIDEO fields, need all three depth files for each VIDEO tile
if field_name == 'XMM' or field_name == 'CDFS':

    # Generate tile names 
    tile_names = [field_name+x for x in ['1', '2', '3']]

    tile1_depth = Table.read(depth_dir / tile_names[0] / 'phot' / f'{det_filter}_{aperture_size}as_gridDepths_300_200.fits')
    tile2_depth = Table.read(depth_dir / tile_names[1] / 'phot' / f'{det_filter}_{aperture_size}as_gridDepths_300_200.fits')
    tile3_depth = Table.read(depth_dir / tile_names[2] / 'phot' / f'{det_filter}_{aperture_size}as_gridDepths_300_200.fits')
    depth_files.append(tile1_depth)
    depth_files.append(tile2_depth)
    depth_files.append(tile3_depth)

if field_name == 'COSMOS':
    tile1_depth = Table.read(depth_dir / field_name / 'phot' / f'{det_filter}_{aperture_size}as_gridDepths_300_200.fits')
    depth_files.append(tile1_depth)

#! Get the VIDEO WCS, if necessary
data_dir = Path.cwd().parents[3] / 'data'

wcs_tiles = []
if field_name != 'COSMOS':

    # Generate tile names 
    tile_names = [field_name+x for x in ['1', '2', '3']]

    for i, tile in enumerate(tile_names):

        if field_name == 'XMM':
            with fits.open(data_dir / tile / 'HSC-G_DR3.fits') as hdu:
                w = WCS(hdu[0].header)

        if field_name == 'CDFS':
            with fits.open(data_dir / tile / 'HSC-G.fits') as hdu:
                w = WCS(hdu[0].header)
        wcs_tiles.append(w)

if field_name == 'COSMOS':
    with fits.open(data_dir / 'COSMOS' / 'UVISTA_Y_DR6.fits') as hdul:
        w= WCS(hdul[0].header)
        wcs_tiles.append(w)


#! Mapping from tiles to index
tile_to_index = {'XMM1': 0, 'XMM2': 1, 'XMM3': 2, 'CDFS1': 0, 'CDFS2': 1, 'CDFS3': 2, 'COSMOS': 0}

#! Get the filter transmission curves
filter_dir = Path.home() / 'lephare' / 'lephare_dev' / 'filt'
filt_files = filter_files() 


# Remove filters not in the field
if field_name == 'XMM':
    filt_files.pop('CFHT-u')
    filt_files.pop('CFHT-g')
    filt_files.pop('CFHT-r')
    filt_files.pop('CFHT-z')
    filt_files.pop('f115w')
    filt_files.pop('f150w')
    filt_files.pop('f277w')
    filt_files.pop('f444w')
    filt_files.pop('VIS')
    filt_files.pop('Ye')
    filt_files.pop('Je')
    filt_files.pop('He')
    filt_files.pop('ch1cds')
    filt_files.pop('ch2cds')
if run_type == '':
    filt_files.pop('f115w')
    filt_files.pop('f150w')
    filt_files.pop('f277w')
    filt_files.pop('f444w')
    filt_files.pop('VIS')
    filt_files.pop('Ye')
    filt_files.pop('Je')
    filt_files.pop('He')
    filt_files.pop('ch1cds')
    filt_files.pop('ch2cds')

# Find the index of the detection filter in the filter files
det_index = list(filt_files.keys()).index(det_filter)

det_filter_path = filter_dir / filt_files[det_filter]
det_filter_curve = np.genfromtxt(det_filter_path)

# Get the filter components
det_filter_wlen = det_filter_curve[:, 0]
det_filter_trans = det_filter_curve[:, 1]

# Normalise the transmission curve, divide by maximum value
det_filter_trans /= np.max(det_filter_curve)

# Compute filter area
det_filter_area = simps(det_filter_trans, det_filter_wlen)

# Now loop through the candidates
for i, ID in enumerate(IDs):

    print('Redshifting object', ID, ' which is object', i+1, 'out of', len(IDs))

    #! Properties of the object from the table and the SED

    # Find the row index in the table with this ID
    row_index = np.where(t['ID'] == ID)[0]

    # Open the SED solution file
    file_name = obj_list[i]
    spec_data = parse_spec_file(file_name)

    # Get SED table
    sed = spec_data.get('sed')

    # High-z model is the first component
    sed = sed[0]

    # Get model photometry
    model_phot = spec_data.get('phot')
    
    # Get the detection filter mag
    det_mag_lephare = float(model_phot['col1'][det_index])

    # Convert to flux
    det_flux_lephare = mag_to_flux(det_mag_lephare)

    # Get the RA,DEC 
    ra = t['RA'][row_index]
    dec = t['DEC'][row_index]

    # Get Muv
    Muv = t['Muv'][row_index]
    print('Muv =', round(Muv[0], 2))

    # Get the redshift
    zphot = t['Zphot'][row_index]

    # Get the tile
    if field_name == 'XMM' or field_name == 'CDFS':
        tile = t['video_tile'][row_index][0]
    else:
        tile = 'COSMOS'

    # Get the correct wcs based on the tile
    wcs = wcs_tiles[tile_to_index[tile]]

    # Get the x,y position of the object from the WCS
    x, y = wcs.all_world2pix(ra, dec, 0)
    
    #! Depth at this position
    this_tile_depth = depth_files[tile_to_index[tile]]
    det_depth_here = grid_depths(this_tile_depth, x, y)[0]

    # Convert to flux
    det_flux_here = mag_to_flux(det_depth_here)

    print('HSC-Z depth here =', round(det_depth_here, 2))

    #! Get SED components and compute initial flux

    wlen = sed['lambda']
    sed = sed['flux']
    wlen = np.array([float(w) for w in wlen])
    sed = np.array([float(s) for s in sed])
    sed = mag_to_flux(sed)

    # Interpolate SED to same grid as filters
    det_interpol = np.interp(det_filter_wlen, wlen, sed)

    # Compute initial flux
    det_flux_init = simps((det_interpol * det_filter_trans), det_filter_wlen) / det_filter_area

    # get initial magnitude
    det_mag_init = flux_to_mag(det_flux_init)

    # Print information about this integrated flux compared to real flux
    print('Integrated flux =', det_flux_init)
    print('Real flux =', det_flux_lephare)
    print('Ratio =', det_flux_init/det_flux_lephare)

    #! Begin the redshifting of the SED
    z_init = zphot
    DL_init = cosmo.luminosity_distance(z_init)
    det_mag = det_mag_init
    z = z_init[0]

    while (det_mag < det_depth_here) and (z <= zmax):

        # Increment the redshift
        z += dz

        # Compute the luminosity distance
        DL = cosmo.luminosity_distance(z)

        # Shift the wavelength grid and SED
        wlen_shift = wlen * ((1 + z) / (1 + z_init))
        sed_mag = flux_to_mag(sed)
        sed_shift = sed_mag - (5 * np.log10(DL_init / DL)) + (2.5 * np.log10((1 + z) / (1 + z_init)))
        sed_shift = mag_to_flux(sed_shift)

        # Interpolate and calculate fluxes
        sed_interpol = np.interp(det_filter_wlen, wlen_shift, sed_shift)

        # Compute the flux
        det_flux = simps((sed_interpol * det_filter_trans), det_filter_wlen) / det_filter_area

        # Convert to magnitude
        det_mag = flux_to_mag(det_flux)

    print('zphot =', round(zphot[0], 2))
    print('zmax =', round(z, 2))
    t['zmax'][row_index] = z

    #! Compute Vmax

    # Convert field area into steradians
    field_area_ster = field_area * (np.pi / 180)**2

    # Compute the volume
    V = field_area_ster/3 * (cosmo.comoving_distance(z)**3 - cosmo.comoving_distance(zmin)**3)
    maximum_V = field_area_ster/3 * (cosmo.comoving_distance(zmax)**3 - cosmo.comoving_distance(zmin)**3)

    Vmax = min(V, maximum_V)
    t['Vmax'][row_index] = Vmax
    print('Vmax =', Vmax)
    print('Maximum Vmax =', maximum_V)

print(t)
t.write(cat_dir / cat_name, overwrite=True)
print(f'Saved catalogue to {cat_dir / cat_name}')
