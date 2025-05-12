"""
Measure the non-parametric sizes of euclid galaxies!

Created: Tuesday 11th February 2025.
"""
import numpy as np
import glob
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from pathlib import Path
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.segmentation import detect_threshold, detect_sources

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

def aperture_30pkpc(z, cosmo):
    '''
    Computes the aperture radius for an object at redshift z such that 30kpc is enclosed.
    '''
    
    ASperkpc = cosmo.arcsec_per_kpc_proper(z)
    aperture_size = (30 * ASperkpc.value) / 2
    return aperture_size

#! Switches
plot = False

#! Directories
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
stamp_dir = Path.cwd().parents[1] / 'data' / 'stamps'

#! Euclid filters to measure sizes on
filter_names = ['Y', 'J', 'H']

psf_fwhms = {'Y': 0.48, 'J': 0.51, 'H': 0.55}

# Define cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Pixel scale
pixel_scale = 0.15 # arcsec/pixel


# Get all the Y-band images in the stamps dir
stamp_names = glob.glob(str(stamp_dir / 'Y*.fits'))

# From this, get all the IDs
IDs = [int(stamp.split('/')[-1].split('Y_')[-1].split('.fits')[0]) for stamp in stamp_names]

# Load the catalogue, brightest to faintest
t = Table.read(cat_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_01_31.fits')
t.sort('Muv')

# Testing
#t = t[t['ID'] == 631099]

# Make a size column and previous ID column, if they don't exist
if 'rY_nonparam' not in t.colnames:
    t['rY_nonparam'] = -99 * np.ones(len(t))
if 'rJ_nonparam' not in t.colnames:
    t['rJ_nonparam'] = -99 * np.ones(len(t))
if 'rH_nonparam' not in t.colnames:
    t['rH_nonparam'] = -99 * np.ones(len(t))

if 'ID_prev' not in t.colnames:
    t['ID_prev'] = ['none'] * len(t)  # Initialize with a longer string format
    t['ID_prev'] = t['ID_prev'].astype('<U50')  # Ensure it can hold up to 50 characters

# Sort IDs from the file names to match the order of the catalogue. the IDs are a subset of the catalogue
cat_IDs = [ID for ID in t['ID']]
IDs = [ID for ID in cat_IDs if ID in IDs]

# And reference catalogue for existing objects in the sample
t_ref = Table.read(cat_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_01_31_XMATCH_WITH_LITERATURE.fits')

sizes_in_Y = []
sizes_in_J = []
sizes_in_H = []

# List of objects to apply custom aperture sizes to, and to allow additional segments to remain unmasked
custom_IDs = []
custom_apertures = []
additional_segments = []



# Loop through IDs
for i, ID in enumerate(IDs):

    # Get the redshift of this ID from the table
    z = t['Zphot'][t['ID'] == ID][0]
    print(f'ID: {ID}, z: {z}')

    # If this ID is in the reference table, get the object name and add it to prev ID
    if ID in t_ref['ID']:
        t['ID_prev'][t['ID'] == ID] = t_ref['Object Name'][t_ref['ID'] == ID][0]

    # Aperture corresponding to 30pkpc, at this redshift
    aperture_size_pix = aperture_30pkpc(z=z, cosmo=cosmo) / pixel_scale
    print(f'Aperture size: {aperture_size_pix} pixels')

    # Loop through filters
    for j, filter_name in enumerate(filter_names):

        # Open the stamp
        stamp = fits.open(stamp_dir / f'{filter_name}_{ID}.fits')[0].data

        # Show the stamp
        # if plot:
        #     plt.imshow(stamp, origin='lower', vmin=-0.1, vmax=0.8)
        #     plt.title(f'{ID}')
        #     aperture = CircularAperture((stamp.shape[0] / 2, stamp.shape[1] / 2), r=aperture_size_pix)
        #     aperture.plot(color='red')
        #     plt.show()
        #     plt.close()
        
        # Detect sources in the image
        threshold = detect_threshold(stamp, nsigma=2.0, background=0.)

        # Detect sources
        segm = detect_sources(stamp, threshold, npixels=5)

        if segm is None:
            r_kpc = -99
            continue

        # Show segmentation map
        if plot:
            plt.title(f'{ID}, {filter_name}')
            plt.imshow(segm, origin='lower')    
            aperture = CircularAperture((stamp.shape[0] / 2, stamp.shape[1] / 2), r=aperture_size_pix)
            aperture.plot(color='red')   
            plt.show()
            plt.close()

        # Find value of segment in middle of image
        segm_val = segm.data[int(stamp.shape[0] / 2), int(stamp.shape[1] / 2)]
        print(segm_val)

        # If it's zero, check the pixels around central one for a non-zero value
        if segm_val == 0:
            for i in range(-3, 3):
                for j in range(-3, 3):
                    segm_val = segm.data[int(stamp.shape[0] / 2) + i, int(stamp.shape[1] / 2) + j]
                    if segm_val != 0:
                        break
        print(segm_val)

        # mask the other segments
        mask = (segm.data != segm_val)

        # And mask things outside aperture
        xx, yy = np.meshgrid(np.arange(stamp.shape[0]), np.arange(stamp.shape[1]))
        dist = np.sqrt((xx - stamp.shape[0] / 2)**2 + (yy - stamp.shape[1] / 2)**2)
        mask = mask | (dist > aperture_size_pix)

        # Show the stamp
        # if plot:
        #     plt.imshow(stamp, origin='lower', vmin=-0.1, vmax=0.8)
        #     plt.title(f'{ID}')
        #     aperture = CircularAperture((stamp.shape[0] / 2, stamp.shape[1] / 2), r=aperture_size_pix)
        #     aperture.plot(color='red')
        #     plt.show()
        #     plt.close()

        # Make copy of stamp with mask applied, setting masked pixels to np.nan
        stamp_copy = stamp.copy()
        stamp_copy[mask] = np.nan

        # Get flux array of valid pixels in aperture
        fluxes = np.sort(stamp_copy.flatten())[::-1]
        fluxes = fluxes[~np.isnan(fluxes)]

        # Get the total flux from aperture photometry
        aperture = CircularAperture((stamp.shape[0] / 2, stamp.shape[1] / 2), r=aperture_size_pix)
        phot_table = aperture_photometry(stamp, aperture, mask=mask)
        flux_in_aperture = phot_table['aperture_sum'][0]

        # Now go through the pixels, counting flux until we reach 50%
        cum_sum = 0
        num_pixels = 0

        for flux in fluxes:
            cum_sum += flux
            num_pixels += 1
            if cum_sum > 0.5 * flux_in_aperture:
                break
        
        # Get the radius of the 50% flux
        radius_50 = np.sqrt(num_pixels / np.pi)

        # Draw the 50% radius on the stamp
        if plot:
            plt.imshow(stamp, origin='lower', vmin=-0.1, vmax=0.8)
            plt.title(f'{ID}, {filter_name}')
            aperture = CircularAperture((stamp.shape[0] / 2, stamp.shape[1] / 2), r=aperture_size_pix)
            aperture.plot(color='orange')

            aperture = CircularAperture((stamp.shape[0] / 2, stamp.shape[1] / 2), r=radius_50)
            aperture.plot(color='red', linewidth=2)

            # Plot black region corresponding to mask
            plt.imshow(mask, origin='lower', alpha=0.5, cmap='Greys')

            plt.show()
            plt.close()

        # Compute the effective radius
        area = num_pixels * pixel_scale**2 # in arcsec^2
        r_arcsec = np.sqrt(area / np.pi)

        # Unbroaden by PSF
        r_arcsec = np.sqrt(r_arcsec**2 - (psf_fwhms[filter_name]/2)**2)

        # Convert to kpc
        as_per_kpc = cosmo.arcsec_per_kpc_proper(z).value
        r_kpc = r_arcsec / as_per_kpc

        # Print the effective radius
        print(f'Effective radius: {r_kpc} kpc')

        if j==0:
            sizes_in_Y.append(r_kpc)
            # fill in correspinding column in table
            t['rY_nonparam'][t['ID'] == ID] = r_kpc
        elif j==1:
            sizes_in_J.append(r_kpc)
            t['rJ_nonparam'][t['ID'] == ID] = r_kpc
        elif j==2:
            sizes_in_H.append(r_kpc)
            t['rH_nonparam'][t['ID'] == ID] = r_kpc


# Save table
t.write(cat_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_01_31_nonparam_sizes.fits', overwrite=True)

sizes_in_Y = np.array(sizes_in_Y)
sizes_in_J = np.array(sizes_in_J)
sizes_in_H = np.array(sizes_in_H)

#plt.hist(sizes_in_Y[sizes_in_Y>0], bins=np.arange(0, 6, 0.25), alpha=0.5, color='blue', label=r'$Y_E$', histtype='step', linewidth=3)
plt.hist(sizes_in_J[sizes_in_J>0], bins=np.arange(0, 15, 0.25), alpha=0.5, color='green', label=r'$J_E$', histtype='step', linewidth=3)
#plt.hist(sizes_in_H[sizes_in_H>0], bins=np.arange(0, 6, 0.25), alpha=0.5, color='red', label=r'$H_E$', histtype='step', linewidth=3)
plt.xlabel(r'$r_e$ (kpc)')
plt.ylabel('N')
plt.legend()
plt.show()
