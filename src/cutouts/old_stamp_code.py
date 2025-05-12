#!/usr/bin/env python3

"""
cut_stamps.py

Cut postage stamps of objects, and save pdf.

Make interactive? So that I can do visual checks quickly...

Created: Wednesday 1st June 2022

"""

''' Import modules '''

import matplotlib as mpl
#mpl.use('Agg') # This is now done automatically based on the 'search' clause.
mpl.rcParams['figure.dpi'] = 100

import os
from os.path import exists
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
import astropy.units as u
from new_catalogue_codes import cut_out, which_field, label_ct
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.ndimage as ndimage
from scipy.stats import norm
from scipy import stats
from matplotlib.patches import Circle

''' Set-up '''

# Read catlaogue name and other things in from lephare input txt file
lephareFile = 'lephare_params_DR3.txt' # XMM Y+J
#lephareFile = 'lephare_params.txt' # CDFS Y+J
#lephareFile = 'COS_YJ_lephare_params.txt' # COSMOS Y+J

#lephareFile = 'COS_Z_lephare_params.txt' # COSMOS HSC-Z

lephareParams = Table.read(lephareFile, format='ascii.commented_header')

# Catalogue
catName = str(lephareParams['catalogue'][0])
print(catName)
# Field
field = str(lephareParams['field'][0])

# Get name for file from the source directory
sourceDir = str(lephareParams['dir_source'][0])
# Or maybe the target dir
#sourceDir = str(lephareParams['dir_target'][0])
targetDir = str(lephareParams['dir_target'][0])

title = sourceDir.split('/')[-2]
print(title)
print(field)
# Plot dir from this dir too.

plotDir = '/mnt/zfsusers/varadaraj/sed_fitting/{0}/plots/stamps/'.format(field) # Primary objects
#plotDir = '/mnt/zfsusers/varadaraj/sed_fitting/{0}/plots/insecure/stamps/'.format(field) # Inclusive candidates
#plotDir = '/mnt/zfsusers/varadaraj/sed_fitting/{0}/plots/lya/stamps/'.format(field) # Lyman alpha objects

# Catalogue directory
if field == 'XMM' or field == 'CDFS':
    catDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/cutting/'.format(field)
if field == 'COSMOS':
    catDir = '/mnt/vardy/vardygroupshare/data/catalogues/finalCOSMOS/cutting/'

#if field == 'XMM':
#    catDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/cutting/'.format(field)
# Catalogue
#if field == 'CDFS':
#    catName = 'CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_04_2_lephare_cut.fits'
#if field == 'XMM':
    #catName = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_lephare_cut.fits'
    #catName = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_irac_z6.5.fits'
    #catName = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_irac_cut_2.fits'
    #catName = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_z6.5_6.fits'
    #catName = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_7.fits'
    #catName = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_BD_TOCHECK.fits'
#    catName = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_fin1_8.fits'

# Full catalogues
if field == 'CDFS':
    catFull = '../CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_07_25.fits'

if field == 'XMM':
    catFull = '../XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_10_26.fits'

if field == 'COSMOS':
    catFull = '../../COSMOS/det_HSC-Z_DR3/COSMOS_DR3_MASKVISTADET_HSC-Z_DR3_2.0as_IRAC2.8as_2023_02_09.fits'


# Maybe object directory
#maybeDir = '/mnt/hoy/temporaryFilesROHAN/lephare/vis_checked/maybe_{0}/'.format(field)
maybeDir = sourceDir + '../2maybe/'

# Image directory
imDir = '/mnt/vardy/vardygroupshare/data/'

# Where to save stamps
outDir = '/mnt/hoy/temporaryFilesROHAN/stamps/'

# Get detection filter
#det = catName.split('MASKVISTADET_')[1].split('_1.8as')[0].split('_2.0as')[0]
det ='YJ'
# Looping through directory
loop = True

# Choose whether to plot. Turn off if we've already generated a pdf, might be quicker.
plot = True

# Smooth the images with a Gaussian kernel
smooth = False
smoothVal = 2 # Sigma for smoothing, pretty sure this is in pixels.

# Switch to do visual check
interactive = False

# Doesn't do anything right now
overwrite = False

# Search for specific objects by ID
search = False

# Saving individual stamps
indiv = False

# Save fits in the search?
save_fits = False

# Set the matplotlib backend to turn off if search is true.
#if search == False:
#    mpl.use('Agg')

# Min/max saturation for plotting.
noise_sat = 2.0
max_sat = 4.0

# Desired filters
if field == 'XMM':
    reqFilters = ['YJ', 'GRI', 'HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1servs', 'ch2servs']
    reqFilters = ['YJ', 'GRI', 'HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1', 'ch2']
    reqFilters = ['GRI', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y', 'Y', 'WIRDS_J', 'J', 'H', 'Ks', 'ch1servs', 'ch2servs'] # 1610530
    #names = ['Y+J', 'GRI', 'G', 'R', 'I', 'NB0816', 'Z', 'NB0921', 'y', 'Y', 'J', 'H', 'Ks', '[3.6]', '[4.5]']
    names = ['GRI', 'NB0816', 'Z', 'NB0921', 'y', 'Y', 'WIRDS J', 'J', 'H', 'Ks', '[3.6]', '[4.5]'] # 1610530

if field == 'CDFS':
    #reqFilters = ['u', 'g', 'r', 'i', 'HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
    reqFilters = ['YJ', 'grGI', 'u', 'g', 'HSC-G', 'r', 'i', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
    names = ['Y+J', 'grGI', 'u', 'g', 'G', 'r', 'i', 'I', 'Z', 'Y', 'J', 'H', 'Ks', '[3.6]', '[4.5]']
    reqFilters = ['Y', 'J']
    names = ['Y', 'J']


if field == 'COSMOS':

    if det == 'YJ':
        # Y+J z=7
        reqFilters = ['YJ', 'GRI_DR3', 'HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
        names = ['Y+J', 'GRI', 'G', 'R', 'I', 'NB0816', 'Z', 'NB0921', 'y', 'Y', 'J', 'H', 'Ks', '[3.6]', '[4.5]']

    if det == 'HSC-Z_DR3':
        # HSC-Z z=6
        reqFilters = ['GRI_DR3', 'HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
        names = ['GRI', 'G', 'R', 'I', 'NB0816', 'Z', 'NB0921', 'y', 'Y', 'J', 'H', 'Ks', '[3.6]', '[4.5]']


# Stamp size
stampSize = 10.0 # Arcseconds on each side

''' EXTRACT IMAGES '''

# Flag crosstalk artefacts. Only need to run once for a catalogue
#if field == 'XMM' or field == 'CDFS':
#    label_ct(catDir+catName, '{0}FULL'.format(field))

#if field == 'COSMOS':
#    label_ct(catDir+catName, 'COSMOS')

#exit()

# Some pixel scale stuff for drawing apertures.
#pixScale = 0.200000000000016*3600
#arcsecperPixel = pixScale * 3600
#pixelperArcsec = 1 / arcsecperPixel
#print(pixScale)


# Loop through objects
if loop:

    # Read in
    hdu = fits.open(catDir+catName)
    header = hdu[0].header
    cat = hdu[1].data

    #with PdfPages(plotDir + '6primary'  + 'fixed_DUD_smooth.pdf') as pdf:
    with PdfPages('tmp.pdf') as pdf:

        for i, obj in enumerate(cat['ID']):

            print('OBJECT NUMBER: {0}'.format(i+1))
            print('OBJECT ID: {0}'.format(obj))

            # Extract the VISTA tile used so we know which image to open.
            if field != 'COSMOS':
                tile = cat['VISTA_tile_used'][i]
            if field == 'COSMOS':
                tile = 'COSMOS'

            print(tile)


#            name = cat[i]['name']

            # If there is UDS IRAC data, use it!
#            if cat['flux_ch1'][i] == -99. or cat['flux_ch2'][i] == -99.:
#                reqFilters[-2] = 'ch1'
#                reqFilters[-1] = 'ch2'

            # Create empty array to save each filter image, for subplots after.
            filterData = []

            # Create subplots
            if det == 'YJ':
                fig, axs = plt.subplots(1, len(reqFilters), figsize=(24,7))
            if det == 'HSC-Z_DR3':
                fig, axs = plt.subplots(1, len(reqFilters), figsize=(24,7))

            # Loop through the filters
            for j, filterName in enumerate(reqFilters):

                # If we have DR3 in the filter name, remove it.
                if filterName.find('DR3') > -1:
                    filterTitle = filterName.split('_')[0]
                else:
                    filterTitle = filterName


                # Get directory
                tileDir = imDir + tile + '/'

                # Open the images.lis to get the image directory
                im = ascii.read(tileDir + 'images.lis', format='commented_header')

                # Get current filter, image location and name
                filt_idx = np.where(im['Name'] == filterName)
                imageLoc = im['directory'][filt_idx][0]
                imageName = im['Image'][filt_idx][0]

                if plot:

                    #print(imageLoc, filterName, imageName)
                    # Open image
                    if imageLoc == 'here':
                        hdu_image = fits.open(tileDir + imageName, memmap=True)
                    if imageLoc != 'here' and filterName[0:2] !='ch':
                        hdu_image = fits.open(imageLoc + imageName, memmap=True)
                    if imageLoc != 'here' and (filterName == 'ch1') or (filterName == 'ch2'):
                        hdu_image = fits.open('/mnt/hoy/temporaryFilesROHAN/movedData/depths/{0}/images/{1}_bgsub.fits'.format(tile, filterName), memmap=True)

                    hdu_im = hdu_image[0]
                    imageHeader = hdu_im.header
                    imageData = hdu_im.data
                    if filterName == 'WIRDS_J':
                        imageHeader = hdu_image[1].header
                        imageData = hdu_image[1].data
                    

                    # Convert RA, DEC into x,y pixels
                    w = WCS(imageHeader)
                    x, y = w.all_world2pix(cat['RA'], cat['DEC'], 1)

                    # Get width of image and cut out!
                    if filterName != 'WIRDS_J':
                        pixScale = abs(imageHeader['CD1_1'])*3600.0
                    else:
                        pixScale = abs(imageHeader['PC1_1'])*3600.0
                    width = int(np.round(stampSize/pixScale))

                    arcsecperPixel = pixScale * 3600
                    pixelperArcsec = 1 / arcsecperPixel
                    #print(pixelperArcsec)

                    # Ensure width is even
                    if width % 2 == 0:
                        halfWidth = width / 2
                    else:
                        halfWidth = (width - 1) / 2

                    halfWidth = int(halfWidth)

                    # Cut
                    xcentre = int(np.round(x[i]))
                    ycentre = int(np.round(y[i]))

                    # Get corner coordinates.
                    xmin = xcentre - halfWidth
                    xmax = xcentre + halfWidth+1
                    ymin = ycentre - halfWidth
                    ymax = ycentre + halfWidth+1

                    # Display image in each filter sequntially, fo checking.
                    data = imageData[ymin:ymax, xmin:xmax]

                    # Draw on Spitzer apertures
#                    if imageName[0:2] == 'ch':

                    #c1 = Circle((width/2, width/2), 1.4*pixelperArcsec*3600)
                    #c1.set_facecolor('none')
                    #c1.set_edgecolor('red')

                    #c2 = Circle((width/2, width/2), 1.4*pixelperArcsec*3600)
                    #c2.set_facecolor('none')
                    #c2.set_edgecolor('red')

                    #axs[-1].add_patch(c1)
                    #axs[-2].add_patch(c2)

                    # Try smoothing?
                    if smooth == True:
                        #data = cv2.GaussianBlur(data, (0.001,0.001), 0) # Doesn' work
                        data = ndimage.gaussian_filter(data, sigma=(smoothVal,smoothVal), order=0)

                    #plt.imshow(data, cmap='gray')

                    #if str(obj) == '393636': why tf did I do this? Oh yeah, for Fig. 4 in my paper.
                    axs[j].set_title(names[j], size=14)


                    #plt.show()

                    # Find stdev of data to set ploting limits. Use median absolute deviation.
                    tmp = data.flatten()
                    #stdev = np.std(tmp)
                    mean =  np.mean(tmp)

                    median = np.median(tmp)
                    MAD = np.median(np.abs(tmp - median))

                    stdev = 1.4826 * MAD

                    # Sigma clip?
                    if np.max(tmp) > 0.0025:
                        tmp = tmp[tmp < 2*stdev]
                        tmp = tmp[tmp > -2*stdev]

                    # Re-fit
                    mean = np.mean(tmp)
                    median = np.median(tmp)
                    MAD = np.median(np.abs(tmp - median))
                    stdev = 1.4826 * MAD

                    filterData = np.append(filterData, data)
                    axs[j].imshow(data, vmin=mean-stdev, vmax=mean+3*stdev, origin ='lower', cmap='binary')
                    axs[j].axis('off')

            # Plot
 #           if plot:

                #plt.text(10, -10, tile + ' ' + str(obj))
                # Label crosstalk live.
#                if cat['CLASS'][i] == 5.0:
#                    plt.text(10, -20, 'Crosstalk flag', color='red')

            # Check if I marked an object as maybe in the visual inspection, and flag this.
#            if exists(maybeDir + 'Id' + str(obj).zfill(9) + '.spec'):
#                plt.text(10, -30, 'Maybe object!', color='orange')

            # ID on side of stamps, for my goldrush BD plot
            #plt.text(-870, 2.7, 'ID '+str(obj), rotation='vertical', size=14)

            plt.show()

#            pdf.savefig(fig)

            if indiv:
                print('Saving to ', plotDir+'{0}.pdf'.format(obj))
#                print('Saving to ', plotDir+'{0}.pdf'.format(name))
                plt.savefig(plotDir+'GOLDRUSH_{0}.pdf'.format(obj), bbox_inches='tight')
#                plt.savefig(plotDir+'{0}.pdf'.format(name), bbox_inches='tight')
            plt.close(fig)

            # Interactive removal of objects
            if interactive:

                # Ask user
                keep = input("Keep this object? y, n, or m for maybe. Type stop to end the code. ")

                # If we keep the object, move to vischeck folder
                if keep == 'y':
                    strID = 'Id' + str(obj).zfill(9) + '.spec'
                    os.system('cp ' + sourceDir + '/' + strID + ' ' + targetDir)
                    #os.system('mv /mnt/hoy/temporaryFilesROHAN/lephare/vis_checked/precheck_{0}/{1} /mnt/hoy/temporaryFilesROHAN/lephare/vis_checked/{2}/'.format(field, strID, field))

                    # Visual check after IRAC and z>6.5
#                    os.system('cp /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{0}_irac_z6.5/{1} /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{2}_vis/'.format(field, strID, field))
#                    os.system('cp /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/2{0}_irac/{1} /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/4{2}_vis/'.format(field, strID, field))

                    #cat['vis_check'][i] = 1.0
                    print('Object kept')

                # Maybe objects go in a separate folder
                if keep == 'm':
                    strID = 'Id' + str(obj).zfill(9) + '.spec'
                    #os.system('mv /mnt/hoy/temporaryFilesROHAN/lephare/vis_checked/precheck_{0}/{1} /mnt/hoy/temporaryFilesROHAN/lephare/vis_checked/{2}/'.format(field, strID, field))
                    os.system('cp ' + sourceDir + '/' + strID + ' ' + maybeDir)

                    # Visual check after IRAC and z>6.5
#                    os.system('cp /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{0}_irac_z6.5/{1} /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{2}_vis/'.format(field, strID, field))
#                    os.system('cp /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/2{0}_irac/{1} /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/4{2}_maybe/'.format(field, strID, field))

                    #cat['vis_check'][i] = 2.0
                    print('Object put into maybe file')

                # If we don't keep it, do nothing.
                if keep == 'n':
                    print('Object not kept')

                # Exit clause
                if keep == 'stop':
                    print('Stopped at object number {0}, ID {1}'.format(i+1, obj))
                    exit()

        # Save final catalogue, with vis_check flag
        #if field == 'CDFS':
            #hdu.writeto('/mnt/vardy/vardygroupshare/data/catalogues/CDFSFULL/cutting/CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_04_2_vischeck_2.fits', overwrite=True)
        #if field == 'XMM':
            #hdu.writeto('/mnt/vardy/vardygroupshare/data/catalogues/XMMFULL/XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_vischeck_2.fits', overwrite=True)

if search:

    tmp = []
    mean = 0.0
    stdev = 0.0

    # Read in catalogue
    hdu = fits.open(catDir+catFull)
    header = hdu[0].header
    cat = hdu[1].data

    ID = input("Type object ID (without leading zeros): ")

    ID = int(ID)

    #strID = ID.zfill(9)

    # Extract the VISTA tile used so we know which image to open.
    if field != 'COSMOS':
        tile = cat['VISTA_tile_used'][cat['ID'] == ID][0]

        print(tile)

    # Create empty array to save each filter image, for subplots after.
    filterData = []

    # Create subplots
    fig, axs = plt.subplots(1, len(reqFilters), figsize=(24,7))
    plt.title(tile + ' ' + str(ID))

    # Loop through the filters
    for j, filterName in enumerate(reqFilters):

        # Get directory
        tileDir = imDir + tile + '/'

        # Open the images.lis to get the image directory
        im = ascii.read(tileDir + 'images.lis', format='commented_header')

        # Get current filter, image location and name
        filt_idx = np.where(im['Name'] == filterName)
        imageLoc = im['directory'][filt_idx][0]
        imageName = im['Image'][filt_idx][0]

        # Open image
        if imageLoc == 'here':
            hdu_image = fits.open(tileDir + imageName, memmap=True)
        if imageLoc != 'here':
            hdu_image = fits.open(imageLoc + imageName, memmap=True)

        hdu_im = hdu_image[0]
        imageHeader = hdu_im.header
        imageData = hdu_im.data

        # Convert RA, DEC into x,y pixels
        w = WCS(imageHeader)
        x, y = w.all_world2pix(cat['RA'][cat['ID']==ID], cat['DEC'][cat['ID']==ID], 1)

        # Get width of image and cut out!
        pixScale = abs(imageHeader['CD1_1'])*3600.0
        width = int(np.round(stampSize/pixScale))

        arcsecperPixel = pixScale
        pixelperArcsec = 1 / arcsecperPixel

        print(arcsecperPixel)

        # Ensure width is even
        if width % 2 == 0:
            halfWidth = width / 2
        else:
            halfWidth = (width - 1) / 2

        halfWidth = int(halfWidth)

        # Cut
        xcentre = int(np.round(x))
        ycentre = int(np.round(y))

        # Get corner coordinates.
        xmin = xcentre - halfWidth
        xmax = xcentre + halfWidth+1
        ymin = ycentre - halfWidth
        ymax = ycentre + halfWidth+1

        # Display image in each filter sequntially, for checking.
        data = imageData[ymin:ymax, xmin:xmax]

        if save_fits:

            cutout = Cutout2D(data=imageData, position=(x,y), size=width, fill_value=np.nan, mode='partial', wcs=w)

            hdu_image[0].data = cutout.data

            hdu_image[0].header.update(cutout.wcs.to_header())

            hdu_image.writeto('VIDEO_z7_28_{0}.fits'.format(filterName))

        # Draw on Spitzer apertures
        hdu_image[0].data = data

        c1 = Circle((width/2, width/2), 1.*pixelperArcsec)
        c1.set_facecolor('none')
        c1.set_edgecolor('red')
        c2 = Circle((width/2, width/2), 1.*pixelperArcsec)
        c2.set_facecolor('none')
        c2.set_edgecolor('red')
        c3 = Circle((width/2, width/2), 1.*pixelperArcsec)
        c3.set_facecolor('none')
        c3.set_edgecolor('red')

        #axs[0].add_patch(c1)
        #axs[1].add_patch(c2)
        #axs[2].add_patch(c3)


        # Try smoothing?
        if smooth == True:
            data = ndimage.gaussian_filter(data, sigma=(smoothVal,smoothVal), order=0)

        #plt.imshow(data, cmap='gray')
        axs[j].set_title(filterName, size=40)
        #plt.show()

        # Find stdev of data to set ploting limits. Use median absolute deviation.
        tmp = data.flatten()
        #stdev = np.std(tmp)
        mean =  np.mean(tmp)

        median = np.median(tmp)
        MAD = np.median(np.abs(tmp - median))

        stdev = 1.4826 * MAD

        # Sigma clip?
        if np.max(tmp) > 0.0025:
            tmp = tmp[tmp < 2*stdev]
            tmp = tmp[tmp > -2*stdev]

        # Re-fit
        mean = np.mean(tmp)
        median = np.median(tmp)
        MAD = np.median(np.abs(tmp - median))
        stdev = 1.4826 * MAD

        filterData = np.append(filterData, data)
        #axs[j].imshow(data, vmin=mean-0.3*stdev, vmax=mean+2*stdev)
        axs[j].imshow(data, vmin=mean-noise_sat*stdev, vmax=mean+max_sat*stdev, origin='lower', cmap='binary')
        axs[j].axis('off')

    # Plot
    #for j, filterName in enumerate(reqFilters):
    #    axs[j].imshow(filterData[j])
    plt.text(10, 70, tile + ' ' + str(ID))
    plt.show()

    plt.clf()
    exit()

    n, bins, patches = plt.hist(tmp, bins=30, density=1)
    plt.title(stdev)
    y = stats.norm.pdf(bins, mean, stdev)
    l = plt.plot(bins, y, 'r--', linewidth=2)
    plt.show()
