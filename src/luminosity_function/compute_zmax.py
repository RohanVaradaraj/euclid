"""
Compute zmax and vmax.

Created: Wednesday 13th November 2024
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

run_type = ''

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

#! SED Fitting folder
# Name of the directory we want to use to make the catalogue
folder = f'det_Y_J_{run_type}_z7' if run_type != '' else 'det_Y_J_z7'

# Get the list of objects that made it through the SED fitting
obj_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits'
obj_list = glob.glob(str(obj_dir / folder / '*.spec'))

# Get the IDs
IDs = [spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0] for spec_file in obj_list]
IDs = [int(ID) for ID in IDs]

#! Catalogue of above objects
# Parent catalogue from which to get fluxes
#cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_28_with_euclid.fits'
cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits'

# Read in the parent catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
t = Table.read(cat_dir / cat_name)

# If there isn't a completeness column, add it
if 'completeness' not in t.colnames:
    t['completeness'] = Column(np.zeros(len(t)), name='completeness', dtype=float)

#! Add a column for zmax if it doesnt exist already
if 'zmax' not in t.colnames:
    t['zmax'] = Column(np.zeros(len(t)), name='zmax', dtype=float)
    t['Vmax'] = Column(np.zeros(len(t)), name='Vmax', dtype=float)

#! Get the COSMOS depth grid
depth_dir = Path.home().parents[1] / 'vardy' / 'vardygroupshare' / 'data' / 'depths' / 'COSMOS' / 'phot'
Y_depths = Table.read(depth_dir / 'Y_1.8as_gridDepths_300_200.fits')
J_depths = Table.read(depth_dir / 'J_1.8as_gridDepths_300_200.fits')

#! Get the filter transmission curves
filter_dir = Path.home() / 'lephare' / 'lephare_dev' / 'filt'
filt_files = filter_files() 
Y_filter_path = filter_dir / filt_files['Y']
J_filter_path = filter_dir / filt_files['J']
Y_filter = np.genfromtxt(Y_filter_path)
J_filter = np.genfromtxt(J_filter_path)

# Get the filter components
Y_filter_wlen = Y_filter[:, 0]
J_filter_wlen = J_filter[:, 0]

# Convert to frequency
# Y_filter_freq = c.value / (Y_filter_wlen * 1e-10)
# J_filter_freq = c.value / (J_filter_wlen * 1e-10)

Y_filter_trans = Y_filter[:, 1]
J_filter_trans = J_filter[:, 1]

# Normalise the transmission curves, divide by 
Y_filter_trans /= np.max(Y_filter)
J_filter_trans /= np.max(J_filter)

# Filter areas
Y_filter_area = simps(Y_filter_trans, Y_filter_wlen)
J_filter_area = simps(J_filter_trans, J_filter_wlen)
# Y_filter_area = simps(Y_filter_trans, Y_filter_freq)
# J_filter_area = simps(J_filter_trans, J_filter_freq)

print('Y filter area:', Y_filter_area)
print('J filter area:', J_filter_area)

#! Open one COSMOS image to get the WCS
with fits.open(Path.home().parents[1] / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS' / 'UVISTA_Y_DR6.fits') as hdul:
    wcs = WCS(hdul[0].header)

#! Redshifting parameters
dz = 0.01
zmin = 6.50
zmax = 7.50

#! Vmax parameters
cosmos_area = 1.7214634933517867 # UVISTA 
#cosmos_area = 0.6536 # Euclid
#covering_fraction = 0.8 # UVISTA
#covering_fraction = 0.84 # Euclid 
covering_fraction = 1.0 # =1 since accounted for in injection recovery
cosmos_area *= covering_fraction

ratios = []

# Loop through the objects in the table and get its SED properties
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
    
    # Get the LePhare Y and J band mags
    Y_mag_lephare = float(model_phot['col1'][-4])
    J_mag_lephare = float(model_phot['col1'][-3])

    # Convert to fluxes
    Y_flux_lephare = mag_to_flux(Y_mag_lephare)
    J_flux_lephare = mag_to_flux(J_mag_lephare)


    # Get the RA,DEC 
    ra = t['RA'][row_index]
    dec = t['DEC'][row_index]

    # Get Muv
    Muv = t['Muv'][row_index]
    print('Muv =', round(Muv[0], 2))

    # Get the redshift
    zphot = t['Zphot'][row_index]

    # Get the x,y position of the object from the WCS
    x, y = wcs.all_world2pix(ra, dec, 0)

    # Get the Y,J depths at this position
    Y_depth_here = grid_depths(Y_depths, x, y)
    J_depth_here = grid_depths(J_depths, x, y)

    # Convert to fluxes and then get the YJ depth
    Y_depth_here = mag_to_flux(Y_depth_here)
    J_depth_here = mag_to_flux(J_depth_here)

    # ! Add depths in quadrature!
    YJ_depth_here = np.sqrt(Y_depth_here**2 + J_depth_here**2)
    YJ_depth_here = flux_to_mag(YJ_depth_here)

    print('YJ depth =', round(YJ_depth_here[0], 2))

    # Get components of the SED model
    wlen = sed['lambda']
    sed = sed['flux']
    wlen = np.array([float(w) for w in wlen])
    sed = np.array([float(s) for s in sed])
    sed = mag_to_flux(sed)

    #! Interpolate the SED to the same grid as the filters
    Y_interpol = np.interp(Y_filter_wlen, wlen, sed)
    J_interpol = np.interp(J_filter_wlen, wlen, sed)
    
    # Remove duplicate wavelengths
    # unique_wlen, unique_indices = np.unique(wlen, return_index=True)
    # unique_sed = sed[unique_indices]

    # # Perform cubic interpolation
    # Y_interpol = interp1d(unique_wlen, unique_sed, kind='cubic', fill_value='extrapolate')(Y_filter_wlen)
    # J_interpol = interp1d(unique_wlen, unique_sed, kind='cubic', fill_value='extrapolate')(J_filter_wlen)


    #! Initial fluxes
    Y_flux_init = simps((Y_interpol * Y_filter_trans), Y_filter_wlen) / Y_filter_area
    J_flux_init = simps((J_interpol * J_filter_trans), J_filter_wlen) / J_filter_area

    # Y_flux_init = np.trapz(Y_interpol * Y_filter_trans, Y_filter_wlen) / Y_filter_area
    # J_flux_init = np.trapz(J_interpol * J_filter_trans, J_filter_wlen) / J_filter_area

    print('Y flux from integration:', Y_flux_init)
    print('Real flux:', Y_flux_lephare)
    print('ratio: ', Y_flux_init / Y_flux_lephare)
    ratios.append(Y_flux_init / Y_flux_lephare)

    # Get mYJ by adding in flux space, then converting
    YJ_flux_init = Y_flux_init + J_flux_init
    mYJ_init = flux_to_mag(YJ_flux_init)
    print('mYJ: ', mYJ_init)

    # Check the local depth is less than object magnitude
    if mYJ_init > YJ_depth_here[0]:
        print('Object is too faint for the depth')
        YJ_depth_here[0] = 26.2 # Replace with global depth
    # Check again
    if mYJ_init > YJ_depth_here[0]:
        print('Object is still too faint for the depth')

    #! Set up the redshifting
    z_init = zphot
    DL_init = cosmo.luminosity_distance(z_init)
    mYJ = mYJ_init
    z = z_init[0]

    # Initialize verbose flag for printing
    verbose = True

    #! Start redshifting the SED
    while (mYJ < YJ_depth_here[0]) and (z <= zmax):

        # Increment the redshift
        z += dz

        # Get the luminosity distance
        DL = cosmo.luminosity_distance(z)

        # Shift the wavelength grid and SED
        wlen_shift = wlen * ((1 + z) / (1 + z_init))
        sed_mag = flux_to_mag(sed)
        sed_shift = sed_mag - (5 * np.log10(DL_init / DL)) + (2.5 * np.log10((1 + z) / (1 + z_init)))
        sed_shift = mag_to_flux(sed_shift)

        # Interpolate and calculate fluxes
        Y_interpol = np.interp(Y_filter_wlen, wlen_shift, sed_shift)
        J_interpol = np.interp(J_filter_wlen, wlen_shift, sed_shift)

        Y_flux = simps((Y_interpol * Y_filter_trans), Y_filter_wlen) / Y_filter_area
        J_flux = simps((J_interpol * J_filter_trans), J_filter_wlen) / J_filter_area

        # Compute combined flux and magnitude
        YJ_flux = Y_flux + J_flux
        mYJ = flux_to_mag(YJ_flux)

    # Add the determined zmax to the table
    print('zphot =', round(zphot[0], 2))
    print('zmax =', round(z, 2))
    t['zmax'][row_index] = z

    #! Determine the Vmax for this galaxy.
    
    # Convert field area to steradians
    field_area = cosmos_area * (np.pi/180)**2

    V = field_area/3 * (cosmo.comoving_distance(z)**3 - cosmo.comoving_distance(6.5)**3)
    maximum_V = field_area/3 * (cosmo.comoving_distance(7.5)**3 - cosmo.comoving_distance(6.5)**3)

    Vmax = min(V, maximum_V)
    t['Vmax'][row_index] = Vmax
    print('Vmax =', Vmax)

# Write the table
t.write(cat_dir / cat_name, overwrite=True)

plt.hist(ratios, bins=np.arange(0, 2, 0.05))
plt.show()






