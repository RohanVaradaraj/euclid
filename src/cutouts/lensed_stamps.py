"""
lensed_stamps.py

Extracts stamps for high-redshift LBGs that have a small angular separation to DEVILS galaxies, thus may experience strong lensing.

Created: Friday 15 December 2023

"""

import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.nddata.utils import Cutout2D
from matplotlib.backends.backend_pdf import PdfPages
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel

import sys
sys.path.append(str(Path.cwd().parent))
from catalogues.devils_cats_classes import *
from luminosity_function.lf_classes import *

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100
plotDir = Path.cwd().parent.parent / 'plots' / 'stamps'
imageDir = Path.home() / 'OneDrive - Nexus365' / 'data' / 'images' / 'ground'
tableDir = Path.cwd().parent.parent / 'tables'
pickleDir = Path.cwd().parent.parent / 'data' / 'pkl'

# --- Set up all the variables ---

# If we want to avoid overwriting the existing PDFs, just set this flag to save to 'tmp.pdf'.
tmp = False

# Centre the stamp on DEVILS?
centre_on_lens = True

# Look at a subsample corresponding to those which have HST stamps
hst_stamps = True
t_hst = Table.read(tableDir / 'example_lenses.fits')
ID_HST = t_hst['ID_source']

fields = ['D10'] #, 'D02'] # D10 is COSMOS, D02 is XMM-LSS

redshifts = [3] #, 4, 5]

separation_angle = 12 # arcsec

# How many multiples of Re around the DEVILS galaxy do we want to mask?
foreground_mask_factor = 0

# Cut at high mass
mass_thresh = 10

stampSize = 24 # arcsec on a side

catDir = Path.cwd().parent.parent / 'catalogues'

## Loop through the fields
for i, fieldName in enumerate(fields):
    print(fieldName)

    # Read in the images
    if fieldName == 'D10':
        hdu = fits.open(imageDir / 'HSC-Z_DR3_COSMOS.fits', memmap=True)
    if fieldName == 'D02':
        hdu1 = fits.open(imageDir / 'HSC-Z_DR3_XMM1.fits', memmap=True)
        hdu2 = fits.open(imageDir / 'HSC-Z_DR3_XMM2.fits', memmap=True)

    # Read DEVILS catalogue containing stellar mass
    devils = Catalogue(Catalogue.catalogue_names[fieldName])
    devils.read_catalogue()

    # Add in coords to the XMM spec catalogue.
    if fieldName != 'D10':
        devils.get_columns_from_other_cat(fieldName, ['RAcen', 'DECcen', 'zBest'])

    # Add in sizes to the spec catalogue
    devils.get_sizes(fieldName)

    # Take the catalogue
    cat = devils.cat

    # Get the correct mass column and redshift column name
    massColumn = Catalogue.mass_colnames[fieldName]

    # Restrict to z_Lens where it is most probable
    cat = cat[cat['zBest'] < 2]

    # Cut at high mass
    massCut = 10 ** mass_thresh
    cat = cat[cat[massColumn] > massCut]

    # Get sizes of DEVILS galaxies
    size_mask = cat['R50_Y'] # already in arcsec

    # Create SkyCoord object
    cDevils = SkyCoord(ra=cat['RAcen']*u.degree, dec=cat['DECcen']*u.degree)

    # Loop through redshifts
    for j, zbin in enumerate(redshifts):

        print(f'z={zbin}')

        # Read in LBG sample in redshift bin
        sample = NathanCat(redshift=zbin)
        sample.read_catalogue(field=devils.fieldDict[fieldName])

        for colname in sample.cat.colnames:
            if ',' in colname:
                new_colname = colname.replace(',', '')
                sample.cat.rename_column(colname, new_colname)

        # Mask to DEVILS region, adding a buffer.
        maskNathan = sample.mask_to_DEVILS(fieldName, buffer=0)

        # Create SkyCoord object
        cNathan = SkyCoord(ra=maskNathan['RA'], dec=maskNathan['DEC'])

        ### Check if crossmatching has already run and saved to pickle ###

        # Pickle file name based on whether we are masking the DEVILS galaxy or not.
        if foreground_mask_factor == 0:
            pickle_filename = pickleDir / f'{fieldName}_z{zbin}_sep{separation_angle}as_size{stampSize}as.pkl'
        else:
            pickle_filename = pickleDir / f'{fieldName}_z{zbin}_sep{separation_angle}as_size{stampSize}as_masked{foreground_mask_factor}Re.pkl'

        if Path(pickle_filename).is_file():
            # If the pickle file exists, load the data from the file
            with open(pickle_filename, 'rb') as file:
                close_idx, devils_idx = pickle.load(file)

        else:

            # Initialize an empty list to store cNathan coordinates within each circle
            close_idx = []

            # Get the devils index to extract its properties
            devils_idx = []

            # Counter for following loop
            counter = 1
            ############## Get all closely separated objects ##############
            for k, devils_coord in enumerate(cDevils):

                counter += 1
                if counter % 1000 == 0:
                    print(f'{k+1} of {len(cDevils)}')
                    counter = 1
                
                # Get the coordinates of the current DEVILS object
                devils_ra, devils_dec = devils_coord.ra.value, devils_coord.dec.value

                # Calculate the separation to all cNathan coordinates
                sep_to_cNathan = cNathan.separation(SkyCoord(ra=devils_ra*u.deg, dec=devils_dec*u.deg))

                # Create a mask for cNathan coordinates within the circle.
                in_circle_mask = (sep_to_cNathan.arcsec <= separation_angle) & ((sep_to_cNathan.arcsec >= foreground_mask_factor * size_mask[k]))

                # Collect lensing devils galaxies
                if np.sum(in_circle_mask) > 0:
                    devils_idx.append(k)

                # Append the indices within the circle to the list
                indices_within_circle = np.where(in_circle_mask)[0]
                close_idx.append(indices_within_circle)

            # Convert the list of indices to a single NumPy array
            close_idx = np.unique(np.concatenate(close_idx))

        # Get the sources and lenses
        sources = maskNathan[close_idx]
        lenses = cat[devils_idx]

        if ~Path(pickle_filename).is_file():
            # Save the data to a pickle file for future use
            with open(pickle_filename, 'wb') as file:
                pickle.dump((close_idx, devils_idx), file)

        ################# Extract stamps ###################
        if tmp:
            plot_filename = plotDir / 'tmp.pdf'
        else:
            if foreground_mask_factor == 0:
                plot_filename = plotDir / f'z{zbin}_{fieldName}_sep_{separation_angle}as.pdf'
            else:
                plot_filename = plotDir / f'z{zbin}_{fieldName}_sep_{separation_angle}as_masked2Re_stamps.pdf'

        if hst_stamps:
            plot_filename = plotDir / f'HST_z{zbin}_{fieldName}_sep_{separation_angle}as.pdf'

        with PdfPages(plot_filename) as pdf:

            for k, obj in enumerate(sources):

                print(f'{k+1} of {len(sources)}')

                ra = sources[k]['RA']
                dec = sources[k]['DEC']
                id = sources[k]['UID']
                zphot = sources[k]['ZPhot']

                # If running the HST stamps, skip over everything else.
                if hst_stamps:
                    if id not in ID_HST:
                        continue

                # Print the source object ID.
                print(f'Object {id}')

                # Take correct XMM tile depending on RA
                if fieldName == 'D02':
                    if ra <  35.0270833:
                        hdu = hdu1
                    else:
                        hdu = hdu2

                ############################## Centering on source high-z galaxy ################################
                if centre_on_lens == False:
                    position = SkyCoord(ra=ra, dec=dec, frame='fk5', unit='deg')

                    imageHeader = hdu[0].header
                    imageData = hdu[0].data

                    w = WCS(imageHeader)

                    pixScale = abs(imageHeader['CD1_1'])*3600.0

                    x, y = skycoord_to_pixel(position, wcs=w)

                    width = int(np.round(stampSize / pixScale))

                    cutout = Cutout2D(data=imageData, position=(x,y), size=width, fill_value=np.nan, mode='partial', wcs=w)

                    # Plotting limits
                    flatStamp = cutout.data.flatten()
                    stampSum = sum(flatStamp)

                    mean = np.mean(flatStamp)

                    median = np.median(flatStamp)
                    MAD = np.median(np.abs(flatStamp - median))

                    stdev = 1.4826 * MAD

                    if np.nanmax(flatStamp) > 0.0025:
                        flatStamp = flatStamp[flatStamp < 2*stdev]
                        flatStamp = flatStamp[flatStamp > -2*stdev]

                    mean = np.mean(flatStamp)
                    median = np.median(flatStamp)
                    MAD = np.median(np.abs(flatStamp - median))

                    stdev = 1.4826 * MAD

                    # Find the index of the closest object in lenses
                    lens_positions = SkyCoord(ra=lenses['RAcen'], dec=lenses['DECcen'], unit='deg')
                    separations = position.separation(lens_positions)
                    closest_lens_idx = np.argmin(separations)

                    zLens = lenses['zBest'][closest_lens_idx]
                    ang_sep = np.min(separations)
                    lens_mass = lenses[massColumn][closest_lens_idx]

                    ra_lens = lenses['RAcen'][closest_lens_idx]
                    dec_lens = lenses['DECcen'][closest_lens_idx]

                    lens_coords = SkyCoord(ra_lens, dec_lens, frame='fk5', unit='deg')
                    x_lens, y_lens = skycoord_to_pixel(lens_coords, wcs=cutout.wcs)

                    # Print the coords of the lens
                    print(f'Coords of lens: {ra_lens, dec_lens}')

                    fig, ax = plt.subplots(figsize=(8, 8))

                    # Annotations
                    zphot_annotation = f'z={zphot:.3}'
                    zlens_annotation = f'z={zLens:.3}'
                    mass_annotation = f'log(M)={round(np.log10(lens_mass), 1)}' + r'$M_{\odot}$'

                    # Check relative positions of the red and blue crosses
                    offset = 25
                    fontsize = 11
                    if y_lens > (width-1)/2:
                        # Blue cross above the red one
                        ax.text((width-1)/2, (width-1)/2 + offset, zphot_annotation, color='red', fontsize=fontsize, ha='center', va='top')
                        ax.text(x_lens, y_lens - offset, f'{zlens_annotation}\n{mass_annotation}', color='deepskyblue', fontsize=fontsize, ha='center', va='bottom')
                    else:
                        # Red cross is to the right of the blue one
                        ax.text(x_lens, y_lens - offset, f'{zlens_annotation}\n{mass_annotation}', color='deepskyblue', fontsize=fontsize, ha='center', va='bottom')
                        ax.text((width-1)/2, (width-1)/2 + offset, zphot_annotation, color='red', fontsize=fontsize, ha='center', va='top')

                    # Plot each stamp on a new page in the PDF
                    ax.imshow(cutout.data, origin='lower', cmap='viridis', vmin=mean - 2*stdev, vmax=mean + 20*stdev)
                    ax.set_title(f'{id}, separation={round(ang_sep.arcsec, 1)}as')
                    ax.scatter((width-1)/2, (width-1)/2, marker='x', color='red', s=50)  # Draw an 'x' at the RA, DEC of the LBG
                    ax.scatter(x_lens, y_lens, marker='x', color='deepskyblue', s=50) # Draw x on the lens
                    ax.axis('off')
                    pdf.savefig()
                    plt.close()

                ############################# Centering on the lens, like the HST stamps ################################
                if centre_on_lens:
                    
                    position = SkyCoord(ra=ra, dec=dec, frame='fk5', unit='deg')

                    imageHeader = hdu[0].header
                    imageData = hdu[0].data

                    w = WCS(imageHeader)

                    pixScale = abs(imageHeader['CD1_1'])*3600.0

                   # Find the index of the closest object in lenses
                    lens_positions = SkyCoord(ra=lenses['RAcen'], dec=lenses['DECcen'], unit='deg')
                    separations = position.separation(lens_positions)
                    closest_lens_idx = np.argmin(separations)

                    zLens = lenses['zBest'][closest_lens_idx]
                    ang_sep = np.min(separations)
                    lens_mass = lenses[massColumn][closest_lens_idx]

                    ra_lens = lenses['RAcen'][closest_lens_idx]
                    dec_lens = lenses['DECcen'][closest_lens_idx]

                    lens_coords = SkyCoord(ra_lens, dec_lens, frame='fk5', unit='deg')
                    x_lens, y_lens = skycoord_to_pixel(lens_coords, wcs=w)

                    # Cut out stamp
                    width = int(np.round(stampSize / pixScale))

                    cutout = Cutout2D(data=imageData, position=(x_lens,y_lens), size=width, fill_value=np.nan, mode='partial', wcs=w)

                    # Pixel coordinate of lens
                    x, y = skycoord_to_pixel(position, wcs=cutout.wcs)

                    # Plotting limits
                    flatStamp = cutout.data.flatten()
                    stampSum = sum(flatStamp)

                    mean = np.mean(flatStamp)

                    median = np.median(flatStamp)
                    MAD = np.median(np.abs(flatStamp - median))

                    stdev = 1.4826 * MAD

                    if np.nanmax(flatStamp) > 0.0025:
                        flatStamp = flatStamp[flatStamp < 2*stdev]
                        flatStamp = flatStamp[flatStamp > -2*stdev]

                    mean = np.mean(flatStamp)
                    median = np.median(flatStamp)
                    MAD = np.median(np.abs(flatStamp - median))

                    stdev = 1.4826 * MAD
 
                    # Print the coords of the lens
                    print(f'Coords of lens: {ra_lens, dec_lens}')

                    fig, ax = plt.subplots(figsize=(8, 8))

                    # Annotations
                    zphot_annotation = f'z={zphot:.3}'
                    zlens_annotation = f'z={zLens:.3}'
                    mass_annotation = f'log(M)={round(np.log10(lens_mass), 1)}' + r'$M_{\odot}$'

                    # Check relative positions of the red and blue crosses
                    offset = 25
                    fontsize = 11
                    if y_lens > (width-1)/2:
                        # Blue cross above the red one
                        ax.text(x, y-offset, zphot_annotation, color='red', fontsize=fontsize, ha='center', va='top')
                        ax.text((width-1)/2, (width-1)/2 + offset, f'{zlens_annotation}\n{mass_annotation}', color='deepskyblue', fontsize=fontsize, ha='center', va='bottom')
                    else:
                        # Red cross is to the right of the blue one
                        ax.text((width-1)/2, (width-1)/2 + offset, f'{zlens_annotation}\n{mass_annotation}', color='deepskyblue', fontsize=fontsize, ha='center', va='bottom')
                        ax.text(x, y-offset, zphot_annotation, color='red', fontsize=fontsize, ha='center', va='top')

                    # Plot each stamp on a new page in the PDF
                    ax.imshow(cutout.data, origin='lower', cmap='viridis', vmin=mean - 2*stdev, vmax=mean + 20*stdev)
                    ax.set_title(f'{id}, separation={round(ang_sep.arcsec, 1)}as')
                    ax.scatter(x, y, marker='x', color='red', s=50)  # Draw an 'x' at the RA, DEC of the LBG
                    ax.scatter((width-1)/2, (width-1)/2, marker='x', color='deepskyblue', s=50) # Draw x on the lens
                    ax.axis('off')
                    pdf.savefig()
                    plt.close()






