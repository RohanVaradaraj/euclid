#!/usr/bin/env python3

"""
depth_algo.py

A much, much better version of my depth code.

Originally created: Wed 27th October 2021
Modified: Monday 25th March 2024.

"""

# Import libraries
import matplotlib as mpl

# Use backend that doesn't display plots to user, to prevent issues when running on the queue.
mpl.use('Agg')

import numpy as np
import os
from astropy.io import fits, ascii
import sep
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy import stats
from astropy.table import Table
from pathlib import Path

'''SETUP'''

# A function for smoothing the distribution of local depths, to stop the mode moving around.
# Convolution acts as a moving average.
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


#! control the field and filter here
fieldName = 'COSMOS'
filterName = 'Y'

# Size of apertures we put down.
apSize = 2.0 # arcsec

# Number of chunks, I am using 61 (=50*(28/23)) or 50 depending on orientation of the image.
n_chunks = 100
# Arcsec speacing between each centre of each 2'' aperture in the grid.
asSpacing = 3.0

# Write results files
write = True

# ! Directories
euclid_dir = Path.home() / 'euclid' / fieldName
calibDir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'bertin_config'
configDir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'HSC_SSP_DR3' / 'config_files'
plotDir = Path.cwd().parent.parent / 'plots'
gridDir = Path.cwd().parent.parent / 'data' / 'grids'
resultsDir = Path.cwd().parent.parent / 'data' / 'results'

# Change to directory containing SExtractor configuration files
os.chdir(configDir)

#! Begin code

#* Read in image info and parse information
info = ascii.read(euclid_dir / 'images.lis', format='commented_header')
filter_info = info[info['Name'] == filterName]

if filter_info['directory'] == 'here':
    imageDir = euclid_dir
    weightDir = euclid_dir
else:
    imageDir = Path(filter_info['directory'][0])
    weightDir = Path(filter_info['directory'][0])

junkDir = imageDir / 'junk'
checkImage =  junkDir / f'{fieldName}_{filterName}_seg.fits'

imageTitle = filter_info['Image'][0]
weightTitle = filter_info['Weight'][0]

zpt = filter_info['zeropoint'][0]

'''RUNNING SEXTRACTOR'''

# Load SExtractor
os.system('module load null sextractor')

# Run SExtractor

# Print the command to check
print('/mnt/users/videouser/sextractor/bin/sex ' + str(imageDir/imageTitle) + f" -CHECKIMAGE_NAME= '{str(checkImage)}' "  + f"-WEIGHT_IMAGE= '{str(weightDir/weightTitle)}' " +'-c rohanparam.sc')

# Run on real image
os.system('/mnt/users/videouser/sextractor/bin/sex ' + str(imageDir/imageTitle) + f" -CHECKIMAGE_NAME= '{str(checkImage)}' "  + f"-WEIGHT_IMAGE= '{str(weightDir/weightTitle)}' " +'-c rohanparam.sc')

# NOTE: I have to run the -WEIGHT_IMAGE from the command line because I was having issues accessing it from rohanparam.sc.
# NOTE: As of 1/12/21, running the checkimage name from the command line to have more control.
# NOTE: Source extractor has been run once, placing catalogues and segmentation maps in the junkDir.

# Read in the segmentation map for masking.
segDir = junkDir
hdulistSeg = fits.open(checkImage, memmap=True)
segHeader = hdulistSeg[0].header
segMap = hdulistSeg[0].data

'''PREPARE IMAGES'''

# Read in and extract image and extract header.
hdulist = fits.open(imageDir/imageTitle, memmap=True)
print('hdulist: ', imageDir/imageTitle)
imageHeader = hdulist[0].header # Header info
image = hdulist[0].data         # Actual image

# Do byte swap operation to prevent SExtractor errors
image = image.byteswap(inplace=True).newbyteorder() # See sep documentation, https://sep.readthedocs.io
segMap = segMap.byteswap(inplace=True).newbyteorder()

# Pass in size of image.
x_size = imageHeader['NAXIS1']
y_size = imageHeader['NAXIS2']

print('Image dimensions are ', x_size, y_size)

# Pass in the degrees per pixel
degperPixel = imageHeader['CD2_2']
print('Degrees per pixel: ', degperPixel)
arcsecperPixel = degperPixel * 3600
print('arcsec per pixel: ', arcsecperPixel)
pixelperArcsec = 1 / arcsecperPixel
print('pixels per arcsec: ', pixelperArcsec)

# Set up appropriately shaped grid.

# Image dimensions ratio
imRatio = y_size / x_size
print('image ratio: ', imRatio)
# Number of gridpoints/apertures along x axis. The number along the y axis will be scaled by the image dimensions ratio.
# Compute based on desired amount of pixel spacing between centres of apertures. 
pixSpacing = asSpacing * pixelperArcsec
print('4.0 Arcsec spacing is {0} pixels'.format(pixSpacing))

n_x = round(x_size / pixSpacing)
print('n_x: ', n_x)

n_y = round(n_x * imRatio)
print('n_y: ', n_y)



'''CREATE APERTURE GRID'''

# Gridpoint spacing in each direction
spacingX = np.linspace(0, x_size-1, n_x)
spacingY = np.linspace(0, y_size-1, n_y)

# Split grid into some number of chunks, and compute fluxes in these!
n_chunksX = n_chunks

n_chunksY = round(n_chunksX * imRatio)

chunksX = np.array_split(spacingX, n_chunksX)
chunksY = np.array_split(spacingY, n_chunksY)

# Create a meshgrid of the chunks, then unravel
depthgridX, depthgridY = np.meshgrid(np.linspace(0, x_size-(x_size/n_chunksX), n_chunksX), np.linspace(0, y_size-(y_size/n_chunksY), n_chunksY))
depthgridX = depthgridX.reshape((np.prod(depthgridX.shape),))
depthgridY = depthgridY.reshape((np.prod(depthgridY.shape),))

# Compute the central coordinate of each chunk.
depthgridX = depthgridX + 0.5*(x_size/n_chunksX)
depthgridY = depthgridY + 0.5*(y_size/n_chunksY)

# Check it worked:
print('DEPTH GRID', depthgridX, depthgridY)



'''CREATE MASK FROM SEGMENTATION MAP AND PLACE APERTURES'''

# Turn segmentation map into a Boolean matrix.
#segMap = np.where(segMap > 0, 1, segMap)
#segMap = segMap.astype(np.uint8)


# List of mu and sigma from chunks
mus = []
sigmas = []
mags = []
numAps = []


'''SPLITTING XMM2'''

# Compute fluxes and Gaussians in chunks
for i in range(len(chunksX)):
    for j in range(len(chunksY)):
        # Create gridpoint mesh in the chunk
        gridpointsX, gridpointsY = np.meshgrid(chunksX[i], chunksY[j])
        # Reshape mesh into 1D arrays
        gridpointsX = gridpointsX.reshape((np.prod(gridpointsX.shape),))
        gridpointsY = gridpointsY.reshape((np.prod(gridpointsY.shape),))

        # Split XMM2 into UDEEP and DEEP based on weight image

        # Index arrays for points satisfying weight condition
        indexX = []
        indexY = []

        # Place apertures and find fluxes
        flux, fluxerr, flag = sep.sum_circle(image, gridpointsX, gridpointsY, apSize/2 * pixelperArcsec, mask=segMap)

        '''REMOVE MASKED APERTURES'''

        # Check for mask flags or truncated apertures and find their indices
        is_masked = (flag & sep.APER_HASMASKED) | (flag & sep.APER_TRUNC) != 0
        flagIndex = [i for i, x in enumerate(is_masked) if x]

        # Delete flagged flux errors
        flux = np.delete(flux, flagIndex)
        gridpointsX = np.delete(gridpointsX, flagIndex)
        gridpointsY = np.delete(gridpointsY, flagIndex)

        '''FITTING'''

        # Print some error array info
        print('Number of apertures is ', len(flux))
        numAps = np.append(numAps, len(flux))


        # Median absolute deviation (MAD)
        median = np.nanmedian(flux)
        MAD = np.nanmedian(np.absolute(flux - median))
        sigma = 1.4826 * MAD

        # We need different sigma clipping regimes for either HSC/VOICE and VIDEO.
        # The VIDEO clipping regime was determined by checking different values.

        # Sigma clipping round 1
        flux = flux[flux <  2 * sigma]
        flux = flux[flux > -2. * sigma]

        # Fit again round 1
        median = np.nanmedian(flux)
        MAD = np.nanmedian(np.absolute(flux - median))
        sigma = 1.4826 * MAD

        # Sigma clipping round 2
        flux = flux[flux <  2 * sigma]
        flux = flux[flux > -2. * sigma]

        # Fit again round 2
        median = np.nanmedian(flux)
        MAD = np.nanmedian(np.absolute(flux - median))
        sigma = 1.4826 * MAD

        # Remove zeros, if there are any
        flux = flux[flux != 0.0]

        # Fit a mean for plot checks
        mu = np.nanmean(flux)

        # Plot histogram with Gaussian as a visual check. Can plot global or local histogram.
        #if (i == 14 and j == 9) or (i == 30 and j == 32):
        #    n, bins, patches = plt.hist(flux, bins=30, density=1)
        #    y = stats.norm.pdf(bins, mu, sigma)
        #    l = plt.plot(bins, y, 'r--', linewidth=2)
        #    plt.title('{0} {1} {2}'.format(round(mu, 2), round(sigma, 2), len(flux)))
        #    plt.legend()
        #    plt.savefig('/mnt/vardy/vardygroupshare/HSC_SSP_DR3/depths/{0}/plots/local_hist_{0}_{1}_chunk_{2}_{3}.pdf'.format(fieldName,  filterName, i, j))
        #    plt.clf()

        #if n_chunks == 1:
        #    n,bins, patches = plt.hist(flux, bins=40, density=1)
        #    y = stats.norm.pdf(bins, mu, sigma)
        #    l = plt.plot(bins, y, 'r--', linewidth=2)
        #    plt.title('{0} {1} {2}'.format(round(mu, 2), round(sigma, 2), len(flux)))
        #    plt.savefig('/mnt/vardy/vardygroupshare/HSC_SSP_DR3/depths/{0}/plots/global_hist_{0}_{1}.pdf'.format(fieldName,  filterName))

        '''COMPUTE LOCAL DEPTH'''

        # The limitng depth is 5 sigma
        limitingCount = 5.0 * sigma

        # Wanted 8sigma depths for cdfs
#        limitingCount = 8.0 * sigma

        # Compute the magnitude!
        limitingMag = -2.5 * np.log10(limitingCount) + zpt

        # If there are not enough apertures, remove the mag. NOTE: This does more harm than good for Y+J.
        if len(flux) < 30:
            limitingMag = np.nan
        print('The limiting magnitude is ', limitingMag)
        # Save the mag
        mags = np.append(mags, limitingMag)


# Stop code if we are running a global depth histogram.
if n_chunks == 1:
    exit()

# Change apSize to a string for formatting file names
apSize = str(apSize)



'''DEPTH MAP'''

# Clear previous plot so that saving figures works!

#plt.cla()
#plt.clf()

# For the depth map colourmap ranges, delete nans
plotMags = mags[np.isfinite(mags)]
# Delete most extreme values to give better dynamic range in final plots.
qlow, qhigh = np.percentile(plotMags, [0.5, 98])
print('PLOTTING LIMS: ', qlow, qhigh)

# Reshape magnitudes array into image rectangle. Rotate and flip to get the right orientation.
mags = np.rot90(mags.reshape(n_chunksX, n_chunksY))
mags = np.flip(mags, 0)

# Depth map
plt.imshow(mags, extent=[0, x_size, 0, y_size], aspect=1, origin='lower', vmin=qlow, vmax=qhigh)
plt.savefig(plotDir / f'{fieldName}_{filterName}_{apSize}as_depth_map.pdf')
plt.colorbar()

#plt.show()



'''SAVE THE DEPTH GRID'''

# Need to do this here because I remove nans from the magnitude array to plot the histogram.

# Grid directory


# Assign NaNs to -99.0
mags[np.isnan(mags)] = -99.0

# Reshape magnitudes into a 1D array
mags = mags.reshape((np.prod(mags.shape),))

# Create a mask column corresponding to mag nans
mask_mag = np.zeros(len(mags))
mask_mag[np.isfinite(mags)] = 1

print('Mag dim: ', mags.shape)
print('X grid dim: ', depthgridX.shape)
print('Y grid dim: ', depthgridY.shape)
print('Mask dim: ', mask_mag.shape)

t = Table([depthgridX, depthgridY, mags,  mask_mag], names=['x', 'y', 'depths', 'mask'], dtype=['f4', 'f4', 'f4', 'f4'])

t.write(gridDir / f'{fieldName}_{filterName}_{apSize}as_grid_depths.fits', overwrite=True)

# Re-assign NaNs
mags[mags == -99.0] = np.nan

'''HISTOGRAM'''

# Remove NaNs so we can plot the histogram
mags = mags[np.isfinite(mags)]

# Re-bin for mag > lim for histogram and remove outliers
lim = 18.
mags = mags[mags > lim]
plotBins = np.arange(mags.min(), mags.max()+0.1, 0.025)

# Median
medianMag = round(np.nanmedian(mags), 2)

# Mean
meanMag = round(np.nanmean(mags), 2)

# Bin data for mode
bins = np.arange(mags.min()-.1, mags.max()+.1, 0.01)
binnedMags = np.digitize(mags, bins)
modeMag = stats.mode(binnedMags, axis=None)

# Mode
modeMag = round(bins[modeMag[0]][0], 2)

# Clear previous plot so that saving figures works!
plt.cla()
plt.clf()

# Histogram
n, binsMode = np.histogram(mags, bins=bins)
plt.step(binsMode[0:-1], n, color='orange', alpha=0.6)

# Smooth the mags for calculation.
smoothMags = smooth(n, 20)

# Split the peaks for UltraVISTA
split = int(0.99*len(smoothMags))

# Histogram and plots of smoothed mags.
vals, valBins = np.histogram(smoothMags[0:split], bins=bins)
plt.plot(binsMode[0:-1], smoothMags, color='black', linestyle='None', marker='.', linewidth=1, markersize=2)

# Mode
modeMag = round(binsMode[np.argmax(smoothMags)], 2)

qlow, qhigh = np.percentile(mags, [2, 98])

plt.axvline(medianMag, color='r', linestyle='dashed', label='median={0}'.format(medianMag))
plt.axvline(meanMag, color='green', linestyle='dashed', label='mean={0}'.format(meanMag))
plt.axvline(modeMag, color='orange', linestyle='dashed', label='mode={0}'.format(modeMag))
#if cosmos == True and video == True:
#    plt.axvline(modeMag2, color='purple', linestyle='dashed', label='ultradeep mode={0}'.format(modeMag2))
plt.title('{0} {1}'.format(fieldName, filterName))
plt.xlim(left=qlow)
plt.xlim(right=qhigh)
plt.xlabel('mag')
plt.ylabel('Count')
plt.legend()

plt.savefig(plotDir / f'{fieldName}_{filterName}_{apSize}as_depth_hist.pdf'.format(imageDir, filterName, fieldName, apSize))

'''WRITE DEPTH TO FILE'''

measure1 = modeMag
measure2 = medianMag

filename = f'{fieldName}_{filterName}_{apSize}as_depths.txt'


if write:
    f = open(str(resultsDir / filename), 'a')

    f.write('#ap	{0}as	type'.format(apSize))
    f.write('\n')
    f.write('full	{0}	median'.format(measure2))
    f.write('\n')
    f.write('full	{0}	mode'.format(measure1))

