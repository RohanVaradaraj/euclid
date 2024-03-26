import numpy as np
from pathlib import Path
import os
from astropy.io import fits
from astropy.table import Table, hstack, join, Column
import time
import sep
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import gridspec
from scipy.stats import norm
from astropy import units as u
import matplotlib.backends.backend_pdf
from new_catalogue_codes import return_instrips, mask_column
from scipy.stats import norm
from scipy.signal import find_peaks
from astropy.wcs import WCS
from typing import Tuple, Union, List
import math

def image_depth(image_name: str, zeropoint: float, ap_diametersAS: np.ndarray = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), 
                wht_name: str = 'NONE', wht_type: str = 'NONE', 
                IRACap_diametersAS: np.ndarray = np.array([2.8, 3.8, 5.8, 9.8, 11.6]), seg_name: str = 'NONE', 
                output_dir: str = 'none', filter_name: str = 'NONE', 
                num_apertures: int = 300, step: int = 200, overwrite: bool = False, 
                input_sex: Path = Path.home().parent.parent / 'vardy' /'vardygroupshare' / 'data' / 'bertin_config' / 'video_mine.sex', 
                strips: bool = False, bgSub: bool = True, mask: str = 'none', gridSepAS: float = 1.0) -> None:
    """
    Compute depths across an image.

    Args:
        image_name (str): Path to the input image.
        zeropoint (float): Zeropoint for magnitude calculations.
        ap_diametersAS (np.ndarray, optional): Array of aperture diameters in arcseconds. Defaults to [1.8, 2.0, 3.0, 4.0, 5.0].
        wht_name (str, optional): Path to the weight file. Defaults to 'NONE'.
        wht_type (str, optional): Type of weight file. Defaults to 'NONE'.
        IRACap_diametersAS (np.ndarray, optional): Array of IRAC aperture diameters in arcseconds. Defaults to [2.8, 3.8, 5.8, 9.8, 11.6].
        seg_name (str, optional): Path to the segmentation map. Defaults to 'NONE'.
        output_dir (str, optional): Output directory path. Defaults to 'none'.
        filter_name (str, optional): Filter name. Defaults to 'NONE'.
        num_apertures (int, optional): Number of apertures. Defaults to 300.
        step (int, optional): Step size. Defaults to 200.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.
        input_sex (Path, optional): Path to the sextractor configuration file. Defaults to vardygroupshare/bertin_config/video_mine.sex.
        strips (bool, optional): Whether to use strips. Defaults to False.
        bgSub (bool, optional): Whether to perform background subtraction. Defaults to True.
        mask (str, optional): Path to the mask file. Defaults to 'none'.
        gridSepAS (float, optional): Grid separation in arcseconds. Defaults to 1.0.

    Returns:
        None
    """


    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    

    os.environ['EXTRACTOR_DIR'] = str(Path.home().parent.parent / 'users' / 'videouser' / 'sextractor' / 'share' / 'sextractor')
    
    if filter_name[0:2] == 'ch':
        print('Using IRAC apertures.')
        ap_diametersAS = IRACap_diametersAS

    # main code where all the action happens!
    # first define the output files

    if output_dir == 'none':
        # put everything here!
        # make a new directory
        output_dir = 'depths/'
        if os.path.isdir(output_dir) == False:
            os.system('mkdir ' + output_dir)

    output_dir = Path(output_dir)
    # for the seg map
    image_dir = output_dir / 'images'
    if os.path.isdir(image_dir) == False:
        os.system('mkdir ' + str(image_dir))

    plotDir = output_dir / 'plots'
    if os.path.isdir(plotDir) == False:
        os.system('mkdir ' + str(plotDir))

    catDir = output_dir / 'catalogues'
    if os.path.isdir(catDir) == False:
        os.system('mkdir ' + str(catDir))

    aperDir = output_dir / 'phot'
    if os.path.isdir(aperDir) == False:
        os.system('mkdir ' + str(aperDir))

    resultsDir = output_dir / 'results'
    if os.path.isdir(resultsDir) == False:
        os.system('mkdir ' + str(resultsDir))
        
    # make a sub-directory here?  no but useful
    parts = image_name.split('/')
    if len(parts) > 1:
        baseName = parts[-1]

    if baseName[0:5] == 'primer':
        pparts = baseName.split('_')
        baseName = pparts[4]
        print("The base name is ", baseName)

    else:
        # also remove the cat bit
        pparts = baseName.split('.')
        baseName = pparts[0]
        print("The base name is ", baseName)
        
    # check all the necessary files exist
    imyes = os.path.isfile(image_name)
    if imyes == False:
        print("ERROR: the image file does not exist, check the path:\n" + image_name)
        exit()

    if wht_type != "NONE":
        whtyes = os.path.isfile(wht_name)
        if whtyes == False:
            print("ERROR: the weight file does not exist, check the path:\n" + wht_name)            
            exit()
            
    ####################################################
    # Get the seg map

    if filter_name == 'NONE':
        filter_name = baseName

    # get the pixel scale
    hdulist = fits.open(image_name)
    imageHeader = hdulist[0].header

    if 'CD1_1' in imageHeader:
        cdone_o = -3600.0*imageHeader['CD1_1']
    else:
        cdone_o = 3600.0*np.abs(imageHeader['CDELT1'])
    pix_scale = round(cdone_o, 5)
    print("The pixel scale is {0:.4f}".format(pix_scale))
    
    if seg_name == 'NONE':
        
        # run source extractor
        ## define aperture sizes
        apStringPix = str(ap_diametersAS[0]/pix_scale) 
        for ai in range(1, ap_diametersAS.size):
            apStringPix = apStringPix + ',' + str(ap_diametersAS[ai]/pix_scale)
            
        # now run sex
        seg_name = image_dir / str(filter_name + '_seg.fits')
        bgSubName = image_dir  / str(filter_name + '_bgsub.fits')
        outputCatalogue = catDir / str('d' + filter_name + '.fits')

        keywordsbase = ' -CATALOG_TYPE FITS_1.0 -CATALOG_NAME ' + \
                       str(outputCatalogue) + \
                       ' -MAG_ZEROPOINT '+ str(zeropoint) + \
                       ' -WEIGHT_TYPE ' + wht_type + \
                       ' -WEIGHT_IMAGE ' + str(wht_name)
        
        keywords = keywordsbase + \
                   ' -CHECKIMAGE_TYPE "-BACKGROUND,SEGMENTATION" '\
                   +'-CHECKIMAGE_NAME "' + \
                   str(bgSubName) + ',' + str(seg_name) + '" -PHOT_APERTURES ' \
                   + apStringPix
    
        command = '/mnt/users/videouser/sextractor/bin/sex '+ str(image_name) +' -c ' + str(input_sex) + keywords
        print(command)
        
        if os.path.isfile(bgSubName) == False or os.path.isfile(seg_name) == False or overwrite:
            print("The SEG and BG subtracted map do not exist.  Running...")
            os.system(command)
        else:
            print("The SEG and BG subtracted map exist.  Moving on...")

    #######################################################################
    # Next step is to place apertures down
    aperPhotFile = aperDir / str(filter_name + '_aperPhot.fits')
    #overwrite = True
    if os.path.isfile(aperPhotFile) == False or overwrite == True:

        # define the grid separation
        gridSepPixels = np.ceil(gridSepAS/pix_scale) # 5'' separation
#        gridSepPixels = 10.0
        
        # make this tunable...
        if bgSub == False:
            bgSubName = image_name
            print("Not using bg subtracted image.")
        print("Measuring the aperture photometry.")

        ii = str(seg_name).rfind('NIRSPEC')
        if ii > 0:
            field = 'NIRSPEC'
        else:
            field = 'NONE'
            
        aperture_photometry_blank(bgSubName, seg_name, wht_name, ap_diametersAS, grid_separation = gridSepPixels, clean = True, output_fits_name = aperPhotFile, image_dir = image_dir, field = field, overwrite = overwrite)

    
    #######################################################################
    # Then calculate the local depths, and make a nice plot
    # if COSMOS, I need to run in strips too
    recalculate = True

    # mask
    regions, globaldepths, meddepths, modedepths = extract_local_depths(aperPhotFile, ap_diametersAS, zeropoint, recalculate = recalculate, num_apertures = num_apertures, step = step, plotDir = str(plotDir), strips = strips, maskreg = mask, refimage = bgSubName) #, plot = True)
    
    ######################################################################
    # make a nice file with the output
    depthFile = resultsDir + '{0}_{1}.txt'.format(filter_name, num_apertures)
    f = open(depthFile, 'w')
    
    apString = ''
    for i in range(ap_diametersAS.size):
#        apString = apString + '{0:.1f}as\t'.format(ap_diametersAS[i])

        apString = apString + '{0:.2f}as\t'.format(ap_diametersAS[i]) # Change to 2sf for jwst apertures.


    print('######## AP STRING ' + apString + '################')

    f.write('#ap\t{0}\t{1}\n'.format(apString, 'type'))

    depthtype = ['median', 'global', 'mode']
    for di, deptht in enumerate(depthtype):
        
        for r, reg in enumerate(regions):
            
            apResultString = ''
            for i in range(ap_diametersAS.size):
                if deptht == 'median':
                    apResultString = apResultString + '{0:.2f}\t'.format(meddepths[r, i])
                elif deptht == 'global':
                    apResultString = apResultString + '{0:.2f}\t'.format(globaldepths[r, i])
                elif deptht == 'mode':
                    apResultString = apResultString + '{0:.2f}\t'.format(modedepths[r, i])
                    
            printString = '{0}\t{1}\t{2}\n'.format(reg, apResultString, deptht)
            f.write(printString)

    f.close()
    print("Output file saved to ", depthFile)
    
    return



def get_depths(field_name: str, req_filters: list, queue: str = 'none',
               ap_diametersAS: np.ndarray = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), 
               data_dir: Path = Path.home() / 'euclid', output_dir: str = 'none', 
               overwrite: bool = False) -> None:
    """
    Runs the depth code on the glamdring queue.

    Args:
        field_name (str): Name of the field.
        queue (str, optional): Queue name for running in parallel. Defaults to 'none'.
        req_filters (list, optional): List of required filters.
        ap_diametersAS (np.ndarray, optional): Array of aperture diameters in arcseconds. Defaults to [1.8, 2.0, 3.0, 4.0, 5.0].
        data_dir (Path, optional): Directory containing data. Defaults to Path.home() / 'euclid'.
        output_dir (str, optional): Output directory path. Defaults to 'none'.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.

    Returns:
        None
    """

    # Default 
    gridSepAS = 3.0

    strips=False # Placeholder - old code has this for UVISTA strips

    # loop through each filter
    # for the queue run, run each as a separate file...

    for fi, filter_name in enumerate(req_filters):

        # Read in the images file
        dirHere = data_dir / filter_name / field_name
        print(dirHere)
        imagedata = read_image_lis(dirHere)
        availableFilters = np.array(imagedata['Name'])
        print("The available filters are ", availableFilters)
        
        #! Loop through each Euclid tile
        for j, tile_name in enumerate(availableFilters):
            # define the images etc to send through
            image_name = imagedata['Image'][j]
            wht_name = imagedata['Weight'][j]
            wht_type = imagedata['wht_type'][j]
            zeropoint = imagedata['zeropoint'][j]
            image_dir = imagedata['directory'][j]
            maskName = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'masks' / field_name / imagedata['mask'][j]
            
            if image_dir == 'here':
                image_dir = data_dir / filter_name / field_name
            
            # Now spawn the depths!
            if queue == 'none':
                print("Running here ")
                        
                image_depth(str(image_dir/image_name), zeropoint=zeropoint, wht_name = str(image_dir/wht_name), wht_type = wht_type, output_dir = output_dir, strips = strips, filter_name = tile_name, overwrite = overwrite, mask = maskName, gridSepAS = gridSepAS, ap_diametersAS = ap_diametersAS)

            else:

                # make an ap diameters string
                ap_diametersAS = np.array(ap_diametersAS)
                ap_diametersASstring = '{0:.2f}'.format(ap_diametersAS[0])
                for i in range(ap_diametersAS.size-1):

                    ap_diametersASstring = ap_diametersASstring + ',{0:.2f}'.format(ap_diametersAS[i+1])
                

                print(ap_diametersASstring)
                print("Spawning in the queue...", queue)
                # make shell script
                tmpName = "tmp_{1}_{0}.sh".format(tile_name, field_name)
                f = open(tmpName, 'w')
                f.write('#!/bin/bash\n')
                f.write('python3 stupid.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(str(image_dir/image_name), str(image_dir/wht_name), wht_type, zeropoint, output_dir, strips, tile_name, overwrite, maskName, gridSepAS, ap_diametersASstring))
                f.close()
                
                # now execute this
                command = "addqueue -c 'tmp_{0}' -m 9 -q {0} -d ./{1}".format(queue, tmpName)
                #print(command)
                os.system('chmod +x {0}'.format(tmpName))
                os.system(command)
        
    return
    
def aperture_photometry_blank(image_name: str, seg_map: str, wht_map: str, ap_size: np.ndarray,
                               grid_separation: float = 100, pix_scale: float = -99.0, next: int = 0,
                               clean: bool = False, output_fits_name: str = 'none', image_dir: str = '',
                               verbose: bool = False, field: str = 'NORMAL', overwrite: bool = False) -> None:
    """
    Place apertures on blank regions of the image, chunked into a grid.

    Args:
        image_name (str): Name of the image.
        seg_map (str): Segmentation map.
        wht_map (str): Weight map.
        ap_size (np.ndarray): Array of aperture diameters.
        grid_separation (float, optional): Grid separation. Defaults to 100.
        pix_scale (float, optional): Pixel scale. Defaults to -99.0.
        next (int, optional): Index of HDU to use from FITS file. Defaults to 0.
        clean (bool, optional): Whether to clean at the aperture phot level. Defaults to False.
        output_fits_name (str, optional): Output FITS file name. Defaults to 'none'.
        image_dir (str, optional): Directory containing images. Defaults to ''.
        verbose (bool, optional): Verbosity. Defaults to False.
        field (str, optional): Field type. Defaults to 'NORMAL'.
        overwrite (bool, optional): Whether to overwrite existing files. Defaults to False.

    Returns:
        None
    """
    
    # first check if output exists
    if os.path.isfile(output_fits_name) and (overwrite == False):

        origTable = Table.read(output_fits_name)
        cols = np.array(origTable.colnames)
        
        # check if all the columns are there
        change = -1
        
        # modify the column names
        for tt, typ in enumerate(['IMAGE', 'SEG', 'WHT']):
            for	ai, apD	in enumerate(ap_size):
                
                oldcolname = '{0}_flux_{1}'.format(typ, ai)
                newcolname = '{0}_flux_{1:.1f}as'.format(typ, apD)
                
                if np.any(oldcolname == cols):
                    origTable.rename_column(oldcolname, newcolname)
                    change = 1
                    print('Renaming column from ', oldcolname, newcolname)
        
        # overwrite the file
        if change > 0:
            print(origTable.colnames)
            origTable.write(output_fits_name, overwrite = True)
            print('aperture phot file updated ', output_fits_name)

            
    # check if the column I want exists yet or not
    if os.path.isfile(output_fits_name) and (overwrite == False):
        
        origTable = Table.read(output_fits_name)
        cols = np.array(origTable.colnames)
        missingAps = np.ones(ap_size.size, dtype = bool)
        
        for ai, apD in enumerate(ap_size):
            reqCol = '{0}_flux_{1:.1f}as'.format('IMAGE', apD)
            kk = np.any(cols == reqCol)
            if kk:
                missingAps[ai] = False
                
        print('Checking the apertures ', missingAps)
        if np.any(missingAps):
            print('I need to re-run adding this aperture', ap_size[missingAps])
        else:
            print('All required aps are present, do not need to run again')
            return

        append = True
        ap_size = ap_size[missingAps]
        
    else:
        append = False
        
    ## Get the pixel scale
    if verbose:
        print(image_name)
    
    hdulist = fits.open(image_name)
    header = hdulist[next].header
    imageData = hdulist[next].data
    
    if pix_scale < 0.0:
        # read from header
        if 'CD1_1' in header:
            cdone_o = -3600.0*header['CD1_1']
        else:
            cdone_o = 3600.0*np.abs(header['CDELT1'])
        pix_scale = round(cdone_o, 5)
        
    
    ## Get the apertures size, in pixels
    ap_sizePix = ap_size/pix_scale
    
    ## First just try a simple grid
    ## grab the dimensions
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    
    ## create arrays of the central coordinates
    numberX = int((naxis1-grid_separation)/grid_separation)
    numberY = int((naxis2-grid_separation)/grid_separation)
    numberApertures = numberX*numberY
    print('The number of apertures is ', numberApertures)
    
    
    #apertureArray = np.zeros([numberApertures, 2])
    xArray = np.zeros(numberApertures)
    yArray = np.zeros(numberApertures)
    halfGrid = grid_separation/2
        
    numberHere = 0
    for xi in range(numberX):
        for yi in range(numberY):

            xArray[numberHere] = halfGrid+xi*grid_separation
            yArray[numberHere] = halfGrid+yi*grid_separation
            numberHere = numberHere + 1

    ## now do aperture photometry on both
    ## setup the apertures
    radii = ap_sizePix/2.0
    if verbose:
        print("Here")

    ## 1) the image
    tic = time.time()
    phot_image = aperture_phot_fast(imageData, xArray, yArray, radii)
    toc = time.time()
    hdulist.close()
    if verbose:
        print("Finished doing the photometry for image in time {0}".format(toc-tic))

    ## 2) the seg
    ## I don't care about interpolation here
    hdulist = fits.open(seg_map)
    segData = hdulist[next].data
    phot_seg = aperture_phot_fast(segData, xArray, yArray, np.array(radii), subpix = 1)    
    hdulist.close()
    if verbose:
        print("Finished doing the photometry for seg in time {0}".format(toc-tic))
    
    ## 3) the wht
    ## to exclude pixels off the edge
    if wht_map[-4:].lower() == 'none':
        print("No weight data. ")
        ## Just use the image instead.
        phot_wht = Table(phot_image, copy = True)
        
        ## absolute these
        for ri, r in enumerate(radii):
            name = 'flux_' + str(ri)
            phot_wht[name] = np.abs(phot_wht[name])
        
    else:
        hdulist = fits.open(wht_map)
        whtData = hdulist[next].data
        phot_wht = aperture_phot_fast(whtData, xArray, yArray, np.array(radii), subpix = 1)
    # centre means a pixel is either in or outside the aperture
        hdulist.close()
        
    ## Save these results to a fits file
    ## I can do cuts etc in another code
    ## to speed this up!!
    
    if output_fits_name == 'none':
        directory = image_dir + 'depths/catalogues/'
        filter_name = seg_map[0:seg_map.rfind('_')]
        output_fits_name = filter_name + '_aperPhot.fits'
        
    # stack the tables
    bigTable = hstack([phot_image, phot_seg, phot_wht], table_names=['IMAGE', 'SEG', 'WHT'], uniq_col_name='{table_name}_{col_name}')
    
    bigTable = Table(bigTable)
    
    # remove columns to make it more streamlined
    bigTable.remove_column('WHT_xcenter')
    bigTable.remove_column('WHT_ycenter')
    bigTable.remove_column('SEG_xcenter')
    bigTable.remove_column('SEG_ycenter')

    if append:
        for ri, r in enumerate(radii):
            bigTable.remove_column('SEG_flux_{0}'.format(ri))
            bigTable.remove_column('WHT_flux_{0}'.format(ri))
    else:
        for ri, r in enumerate(radii):
            if ri != 2:
                bigTable.remove_column('SEG_flux_{0}'.format(ri))
                bigTable.remove_column('WHT_flux_{0}'.format(ri))        

                
        if clean:

            print('Cleaning at the aperture phot level')
            # remove the bad columns here.
            smallNum = 0.0000001
            deltaZero = 1E-13
            apString = '2'
            
            seg_sum = np.array(bigTable['SEG_flux_' + apString])
            wht_sum = np.array(bigTable['WHT_flux_' + apString])
            ap_sum = np.array(bigTable['IMAGE_flux_' + apString])
            
            if field == 'NIRSPEC':
                good_indicies = (seg_sum < 0.5)  & (wht_sum > -1E-28) & (wht_sum < 1E28)
                
            else:
                good_indicies = (seg_sum < 0.5)  & (wht_sum > smallNum) & ((ap_sum > deltaZero) | (ap_sum < -deltaZero))
            bigTable = bigTable[good_indicies]

    # rename the columns with the diameter size
    for tt, typ in enumerate(['IMAGE', 'SEG', 'WHT']):
        for ai, apD in enumerate(ap_size):
            
            oldcolname = '{0}_flux_{1}'.format(typ, ai)
            newcolname = '{0}_flux_{1:.1f}as'.format(typ, apD)
            if oldcolname in np.array(bigTable.colnames):
                print('Renaming column from ', oldcolname, newcolname)
                bigTable.rename_column(oldcolname, newcolname)
            
    bigTable.info

    if append:
        print(bigTable.colnames)
        print(origTable.colnames)
        
        # join with the big table!
        print('Appending to big aper phot table, lengths = {0}, {1}'.format(len(bigTable), len(origTable)))
        bigTable = join(origTable, bigTable, keys = ['IMAGE_xcenter', 'IMAGE_ycenter'])
        print(bigTable.colnames)
        print('After ', bigTable)

    bigTable.write(output_fits_name, overwrite = True)
    print("Aperture table has been saved to ", output_fits_name)
     
    return

def read_image_lis(dirHere: Path) -> Table:
    """
    Read in photometric filters/images/tiles from 'images.lis' file.

    Args:
        dirHere (Path): Directory path.

    Returns:
        Table: Table containing image data.
    """

    # read in filters
    inputFile = dirHere / 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
            
    return Table.read(inputFile, format = 'ascii.commented_header')

def aperture_phot_fast(imageData: np.ndarray, xArray: np.ndarray, yArray: np.ndarray, radii: np.ndarray, subpix: int = 5) -> Table:
    """
    Perform rapid aperture photometry on some predefined positions using sep (SExtractor Python wrapper).

    Args:
        imageData (np.ndarray): Image data.
        xArray (np.ndarray): Array of x coordinates.
        yArray (np.ndarray): Array of y coordinates.
        radii (np.ndarray): Array of aperture radii.
        subpix (int): Subpixel sampling factor (default is 5, which is what SEXtractor uses.).

    Returns:
        Table: Table containing aperture photometry results.
    """
    
    print("In ap phot fast ")
    
    data = imageData.byteswap().newbyteorder()    
    for ri, r in enumerate(radii):
        print("Before phot fast ")

        flux, fluxerr, flag = sep.sum_circle(data, xArray, yArray, r, subpix = subpix)
        print("After ap phot fast ")
        
        if ri < 1:
            # create a table of the results
            phot_apertures = Table([xArray, yArray, flux], names = ['xcenter', 'ycenter', 'flux_0'], dtype = ['f4', 'f4', 'f4'])
            
        else:
            newcolumn = Column(name = 'flux_' + str(ri), data = flux, dtype = 'f4')
            phot_apertures.add_column(newcolumn)
    
    return phot_apertures

def extract_local_depths(inputTableFile: str, ap_diametersAS: np.ndarray, zeropoint: float, step: int = 500, num_apertures: int = 200, strips: bool = False, 
                         plot: bool = True, local: bool = True, recalculate: bool = True, globalplot: bool = True, clean: bool = True, plotDir: str = '', 
                         maskreg: str = 'none', refimage: str = 'none') -> tuple:
    """Extract depths from the input table of aperture photometry.

    Args:
        inputTableFile (str): Path to the input table file.
        ap_diametersAS (np.ndarray): Array of aperture diameters.
        zeropoint (float): Zeropoint value.
        step (int, optional): Step size. Defaults to 500.
        num_apertures (int, optional): Number of apertures. Defaults to 200.
        strips (bool, optional): Whether to use strips. Defaults to False.
        plot (bool, optional): Whether to plot. Defaults to True.
        local (bool, optional): Whether to calculate local depths. Defaults to True.
        recalculate (bool, optional): Whether to recalculate. Defaults to True.
        globalplot (bool, optional): Whether to plot global depths. Defaults to True.
        clean (bool, optional): Whether to clean. Defaults to True.
        plotDir (str, optional): Directory to save plots. Defaults to ''.
        maskreg (str, optional): Mask region. Defaults to 'none'.
        refimage (str, optional): Reference image. Defaults to 'none'.

    Returns:
        tuple: Tuple containing regions, global depth, median local depth, and mode local depth.
    """

    import matplotlib
    matplotlib.use('Agg')

    
    ###########################################
    # important setup
    nirspec = False
    edgebuffer = 0.1 ## 10% of the edge #500
    edgebuffer_full = 0 #3000

    # remove things close to zero!
    deltaZero = 1E-13

    basename = str(inputTableFile)[:-13]
        
    colourArray = ['Blue', 'Green', 'Green', 'Green', 'Green', 'Red', 'Red', 'Red', 'Red', 'Red']
    
    ## Now loop through the different regions of the image
    if strips:
        # for ultravista!
        regions = ['fullimage', 'stripone', 'striptwo', 'stripthree', 'stripfour', 'gap1', 'gap2', 'gap3', 'gap4']
        regions = ['full', 'str1', 'str2', 'str3', 'str4', 'gap1', 'gap2', 'gap3', 'gap4']
        
    else:
        regions = ['full']
        
    global_depth = np.zeros([len(regions), ap_diametersAS.size])
    median_local_depth = np.zeros([len(regions), ap_diametersAS.size])
    mode_local_depth = np.zeros([len(regions), ap_diametersAS.size])
    
    # define a nice figure
    # save the plot somewhere sensible
    # extract the filter_name
    startI = str(inputTableFile).rfind('/') + 1
    endI = str(inputTableFile).rfind('_')
    filter_name = str(inputTableFile)[startI:endI]
    plotName = Path(plotDir) / str(filter_name + '_' + str(num_apertures) + '_{0}.pdf'.format(step))
    pdf = matplotlib.backends.backend_pdf.PdfPages(plotName)
    
    # Loop through the available apertures
    for ai, apDiAS in enumerate(ap_diametersAS):
        
        print("Extracting local and global depths for aperture = ", apDiAS, " as.")
        
        # The apertures are strings _0, _1 etc
        apString = '_{0:.1f}as'.format(apDiAS)
        
        localDepthsFile = basename + str(apDiAS) + 'as_gridDepths_{0}_{1}.fits'.format(num_apertures, step)

        print("File will be saved to ", localDepthsFile)
        
        ## Always read in aperture file
        ## So i can do global depths
        
        ##############################################
        ## Read in the table file
        ## This has the aperture photometry (not the depths)
        inputTable = Table.read(inputTableFile)
        print("Reading file ", inputTableFile)
        
        if clean == False:
            # do the cleaning here
            smallNum = 0.0000001

            ## cut the table depending on apertures
            ## only accept coordinates where the seg map
            ## aperture is blank, and the wht map is > 0
            seg_sum = np.array(inputTable['SEG_flux' + apString])
            wht_sum = np.array(inputTable['WHT_flux' + apString])
            ap_sum = np.array(inputTable['IMAGE_flux' + apString])
            good_indicies = (seg_sum < 0.5)  & (wht_sum > smallNum) & \
                ((ap_sum > deltaZero) | (ap_sum < -deltaZero))
            print("There are ", sum(good_indicies), " good indicies, out of ", len(seg_sum))
            good_indicies = np.where(good_indicies)
            
            ## This table has all the good apertures
            reducedTable = inputTable[good_indicies]
            
        else:
            reducedTable = inputTable
            
        ## Get the x and y coordianates of the apertures
        apX = reducedTable['IMAGE_xcenter'] #/u.pix
        apY = reducedTable['IMAGE_ycenter'] #/u.pix
        
        ## Check if I need to run local depths from the
        ## apertures
        if (os.path.isfile(localDepthsFile) == False) or recalculate:
            ## Files doesn't exist or I want to recalculate it  
            
            #############################################################
            ################## RUNNING LOCAL DEPTHS ####################
            
            xmax = max(apX)
            ymax = max(apY)
            
            numX = np.ceil(xmax/step)
            numY = np.ceil(ymax/step)

            x = min(apX) + np.arange(numX)*step 
            y = min(apY) + np.arange(numY)*step
            print("Step = ", step, " numx, y = ", numX, numY)
            print("Max = ", xmax, ymax, max(apX), min(apX))
            
            # create x, y arrays
            x = np.zeros(1)
            y = np.zeros(1)             
            
            for xi in np.arange(step/2.0, numX*step, step):
                for yi in np.arange(step/2.0, numY*step, step):
                    x = np.append(x, xi)
                    y = np.append(y, yi)
                    
            # I want a constant grid over the image.
            # remove the first elements
            x = x[1:]
            y = y[1:]

            ## Now run local depths at those points
            depthsLocalFull, maskArray = local_depths(reducedTable, apString, x, y, num_apertures, zeropoint = zeropoint, mask = True, sigmaClip = 3.0)
            
            ## Now save these results for faster calculating/plotting in future
            ## Create a table
            localTable = Table([x, y, depthsLocalFull, maskArray], names = ['x', 'y', 'depths', 'mask'], dtype = ['f4', 'f4', 'f4', 'f4'])
            localTable.write(localDepthsFile, overwrite=True)
            print("Local depths saved to ", localDepthsFile)
            
        else:
            # simply restore the results
            localTable = Table.read(localDepthsFile)
               
    # for plotting and median depths, remove negative objects!
        gg = (localTable['mask'] > 0)

        localTable = localTable[gg]
        x = localTable['x']
        y = localTable['y']
        keys = localTable.colnames
        print(keys)
        depthsLocalFull = localTable['depths']
        
    ################################################
    ## Extract median depths etc
    ## And for different sub-regions
    
    ## Loop through the different regions
        for ri, region in enumerate(regions):
        
            print("Calculating depths in region ", region)
                    
            if region != 'full':
                print('Splitting by strip')
                good_indicies = return_instrips(x, y, region = region)
            else:

                if maskreg == 'none':
                    print('PRINTING X: ', x)
                    print('PRINTING Y: ', y)

                    ## put a buffer here
                    maxNumx = max(x)
                    maxNumy = max(y)
                    minNumx = min(x)
                    minNumy = min(y)
                    
                    minNumx = minNumx + edgebuffer*maxNumx
                    minNumy = minNumy + edgebuffer*maxNumy
                    maxNumx = maxNumx - edgebuffer*maxNumx
                    maxNumy = maxNumy - edgebuffer*maxNumy
                    
                    good_indicies = (x > minNumx) & (x < maxNumx) & (y > minNumy) & (y < maxNumy)
                    
                else:
                    # read in a header
                    w = WCS(str(refimage))
                    ra, dec= w.all_pix2world(x,y, 1)

                    print(ra)
                    print(dec)
                    print('Masking with ', maskreg)
                                        
                    hsci = str(refimage).find('HSC')
                    
                    if hsci > -1:
                        hsc = True
                        print('Masking HSC')
                    else:
                        hsc = False


                    if dec[0] < -15.0:

                        # get the directory
                        kk = maskreg.rfind('/')
                        regDir = maskreg[:kk+1]
                        print(regDir)
                        
                        fff = str(refimage).find('HSC-R')
                        if fff > -1:
                            circlesFile = regDir + 'HSC_circle_cdfs_R_xy.reg'
                        else:
                            circlesFile = regDir + 'HSC_circle_cdfs_xy.reg'
                            
                        good_indicies = mask_column(x, y, maskreg, tokeep = True, hsc = hsc, xy = True, circlesFile = circlesFile)
                    else:
                        good_indicies = mask_column(ra, dec, maskreg, tokeep = True, hsc = hsc)
                    

            print('There are {0} good indicies'.format(np.sum(good_indicies)))

            good_indicies = good_indicies & np.logical_not(np.isnan(localTable['depths']))
            
            finalTable = localTable[good_indicies]
            apXregion = finalTable['x'] 
            apYregion = finalTable['y'] 
            depthsLocal = finalTable['depths'] 

            ii = np.logical_not(np.isnan(depthsLocal))
            
            print('Check for NANS', np.isnan(np.sum(depthsLocal)), np.sum(ii), depthsLocal[ii])
            
            # For the median depth of the local depths
            # I want to exclude the edges
            medianLocalDepth = np.median(depthsLocal)

            magmin = medianLocalDepth - 1.2
            magmax = medianLocalDepth + 0.8

            if nirspec:
                magmin = medianLocalDepth - 0.9
                magmax = medianLocalDepth + 1.3
            
            if region == 'full':

                                
                fig = plt.figure(figsize=(6,8))
                gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
                ax = plt.subplot(gs[0])
                plt.axis('equal')
                
                # get colours
                cm = plt.cm.get_cmap('RdYlBu')
                sc = plt.scatter(x, y, s = 5, c = depthsLocalFull, cmap = cm, linewidth = 0.0, vmin = magmin, vmax = magmax)
                
                plt.colorbar(sc)
                plt.title('Local depths for filter {0}\n Aperture diameter is {1:.1f}as'.format(filter_name, apDiAS))
                
                # now make a histogram to go underneath, do this differently if in strips!
                ax = plt.subplot(gs[1])
                binwidth = 0.01
                low = np.round(magmin*100)/100.0
                bins = np.arange(low+binwidth/2.0, magmax+binwidth, binwidth)
                
                if strips:
                    ## plot subsets of the results
                    ## and colour based on strip/gap 
                    
                    ## 1) get the strips
                    ss = return_instrips(apXregion, apYregion)
                    strip_histy, strip_histx, _ = plt.hist(depthsLocal[ss], facecolor= 'blue', alpha = 0.8, density = True, range = [magmin, magmax], bins = bins, label = ' Ultra-deep')

                    gg = return_instrips(apXregion, apYregion, region = 'puregap')# notinstrips       
                    gap_histy, gap_histx, _ = plt.hist(depthsLocal[gg], facecolor= 'red', alpha = 0.8, density = True, range = [magmin, magmax], bins = bins, label = ' Deep')
                    
                    # Make a caption
                    plt.legend(loc = 'upper left', handletextpad = 0.0, fontsize = 8, frameon = False)
                    
                    N  = 10
                    strip_smoothed = np.convolve(strip_histy, np.ones((N,))/N, mode = 'same')
                    plt.plot(bins[:-1]+binwidth/2.0, strip_smoothed, color= 'k', linestyle = ':')
                    gap_smoothed = np.convolve(gap_histy, np.ones((N,))/N, mode = 'same')
                    plt.plot(bins[:-1]+binwidth/2.0, gap_smoothed, color= 'k', linestyle = ':')

                    smoothed = strip_smoothed
                    
                else:


                    print(depthsLocal)
                    print(magmin, magmax)
                    bad = True
                    if bad:
                        
                        histy, hists, _ = plt.hist(depthsLocal, facecolor= colourArray[ri], alpha = 0.8, density= True, range = [magmin, magmax], bins =bins, histtype = 'step')
                        hx, hy, _ = plt.hist(localTable['depths'], facecolor= 'grey', alpha = 0.8, density = True, range = [magmin, magmax], bins = bins, zorder = 1, histtype = 'step')
                        
                        # smooth this
                        # and plot again
                        N = 10
                        smoothx = bins[:-1]+binwidth/2.0
                        smoothed = np.convolve(histy, np.ones((N,))/N, mode = 'same')
                        plt.plot(smoothx, smoothed, color= 'k', linestyle = ':')
                    
            #####################################################
            ######## GLOBAL depths #############
            apertureResults = np.array(reducedTable['IMAGE_flux'+apString])
                        
            ## First clip the data to remove outliers
            ## that will skew the SD
            medianFlux = np.median(apertureResults)
            mad = np.median(abs(apertureResults - medianFlux))
            sigma_mad = 1.4826*mad
            sigmaClip = 2.0
            #print "The mad is ", mad, medianFlux
            
            ## cut the results for the FIRST time
            good_indicies = (apertureResults > medianFlux - sigmaClip*sigma_mad) & \
                            (apertureResults < medianFlux + sigmaClip*sigma_mad)
            
            cutTable = reducedTable[good_indicies]
            apertureResultsTwo = np.array(cutTable['IMAGE_flux'+apString])
            
            medianFlux = np.median(apertureResultsTwo)
            mad = np.median(abs(apertureResultsTwo - medianFlux))
            sigma_mad = 1.4826*mad
            
            ## cut the results for the SECOND time
            good_indicies = (apertureResultsTwo > medianFlux - sigmaClip*sigma_mad) & \
                            (apertureResultsTwo < medianFlux + sigmaClip*sigma_mad)
            
            cutTableTwo = cutTable[good_indicies]
            reducedResults = np.array(cutTableTwo['IMAGE_flux'+apString])

            ## Plot a histogram of the results.. is it a gaussian?!
            n, binsss = np.histogram(reducedResults, 20, density= True)
            
            ## Fit with Gaussian
            (mu, sigma) = norm.fit(reducedResults)

            
            medianFlux = np.median(reducedResults)
            mad = np.median(abs(reducedResults - medianFlux))
            sigma_mad = 1.4826*mad
            global_depth[ri, ai] = return_mag(sigma_mad, zeropoint, sigma = 5.0)

            ####################################################
            # Median local depth
            print("The median local depth is then: MAD = {0:.2f}, Gaussian = {1:.2f}".format(medianLocalDepth, global_depth[ri, ai]))
            median_local_depth[ri, ai] = medianLocalDepth

            # Mode depth
            #            mod = (histy == np.max(histy))
            mod = (smoothed == np.max(smoothed))
            smallbins = bins[:-1]
            mode = smallbins[mod]+binwidth/2.0
            mode = mode[0]
            mode_local_depth[ri, ai] = mode
            
            ## add these depths to the plot
            if region == 'full':
                if strips == False:
                    
                ## plot this global depth over the top!
                   
                    ylim = ax.get_ylim()

                    depthString = 'Mode depth = {0:4.2f}'.format(mode)
                    plt.text(magmin,ylim[1]*0.9, depthString, fontsize = 10, color = 'k')

                    if nirspec == False:
                        plt.plot([mode, mode], ax.get_ylim())
                        plt.plot([global_depth[ri, ai],global_depth[ri, ai]], ax.get_ylim(), 'r') 
                        plt.plot([medianLocalDepth, medianLocalDepth], ax.get_ylim(), 'g-', linewidth = 2.0) 
                        depthString = 'Median local depth = {0:4.2f}'.format(medianLocalDepth)
                        plt.text(magmin, ylim[1]*0.8, depthString, fontsize = 10)
                        depthString = 'Global depth = {0:4.2f}'.format(global_depth[ri, ai])
                        plt.text(magmin,ylim[1]*0.7, depthString, fontsize = 10)
                        # print('Mode = ', mode)
                    else:

                        peaks,_ = find_peaks(smoothed, width = 10)
                        print(peaks)
                        print('Smoothed!', smoothx[peaks])
                        plt.text(magmin+0.8, ylim[1]*0.9, 'Peak depths = {0}'.format(smoothx[peaks]))
                        plt.scatter(smoothx[peaks], smoothed[peaks], color = 'red')
                        
                pdf.savefig(fig)
                
            ## also plot the results from global depths
            if globalplot:
                fig2 = plt.figure(figsize=(6,8))
                n, bins, patches = plt.hist(reducedResults, 20, density = True)
                
                #yll = mlab.normpdf(bins, mu, sigma)
                yll = norm(loc = mu, scale = sigma)
                l = plt.plot(bins, yll.pdf(bins), 'r--', linewidth=1)
                plt.xlabel('Flux/counts')
                plt.ylabel('Number of apertures')
                plt.title('Ap diameter  =' + str(apDiAS) + ' region = '+ region)
                pdf.savefig(fig2)
                                
        plt.close()
        ## loop through the different regions.

    # close the plot
    pdf.close()
    print("Plot saved to ", plotName)
    exit()
    
    ## Before I return the results
    ## Collect the results into a nice table
    ## with ok formatting!
    for x, xi in enumerate(regions): print('{0}, {1:.2f}'.format(xi, global_depth[x,0]))
    
    return regions, global_depth, median_local_depth, mode_local_depth


def local_depths(cleanTable: dict, apString: str, x: np.ndarray, y: np.ndarray, num_apertures: int,
                 zeropoint: float = -99.0, verbose: bool = False, mask: bool = False, sigmaClip: float = 3.0,
                 plot: str = 'none', regFile: str = 'none', fitGauss: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    ''' 
    Code to find the local depth around a given aperture coordinate, or coordinates.
    Using the closest num_apertures apertures.

    Parameters:
    cleanTable (dict): A dictionary containing the required columns from the input table.
    apString (str): A string representing the aperture.
    x (np.ndarray): An array of x-coordinates.
    y (np.ndarray): An array of y-coordinates.
    num_apertures (int): The number of closest apertures to consider.
    zeropoint (float, optional): The zeropoint value. Defaults to -99.0.
    verbose (bool, optional): Whether to print verbose output. Defaults to False.
    mask (bool, optional): Whether to return the mask array. Defaults to False.
    sigmaClip (float, optional): The sigma value for clipping. Defaults to 3.0.
    plot (str, optional): Plotting option. Defaults to 'none'.
    regFile (str, optional): Regular file option. Defaults to 'none'.
    fitGauss (bool, optional): Whether to fit a Gaussian. Defaults to False.

    Returns:
    Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]: An array of local depths or a tuple of local depths and mask array.
    '''

    verbose = True

    ## Check x and y have the same size
    if len(x) != len(y):
        print("Error, the x and y arrays have a different length.")
        print(len(x), len(y))
        exit()
    
    ## extract the required columns from the input table
    apX = np.array(cleanTable['IMAGE_xcenter']) #/u.pix
    apY = np.array(cleanTable['IMAGE_ycenter']) #/u.pix
    number = np.arange(apX.size)
    
    allApertureResults = np.array(cleanTable['IMAGE_flux'+apString])
    
    ## do two lots of sigma clipping
    
    ## Make a nice output array
    localDepths = np.zeros(x.size)
    maskArray = np.zeros(x.size)
    
    print("There are ", x.size, " positions to calculate depths for.")
    
    if regFile != 'none':
        tf = open(regFile, 'w')
        print('Masking with reg file, ', regFile)

    if plot != 'none':
        import matplotlib.backends.backend_pdf
        plotname = plot + '_check_mad.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(plotname)

    # get a minimum separation
    diffx = np.min(np.abs(y - np.roll(y, 1)))
    print('the miminum x separation is ', diffx)
    
    # loop through the positions
    for xi, xpos in enumerate(x):
        ypos = y[xi]

            
        ## calculate the radius        
        deltaX = apX - xpos
        deltaY = apY - ypos
        radius = np.sqrt(deltaX*deltaX + deltaY*deltaY)
        
        ## sort this array, and then the table
        idx = np.argpartition(radius, num_apertures)
        useIndicies = idx[0:num_apertures]
        
        apRedX= apX[useIndicies]
        apRedY= apY[useIndicies]
        
        #print num_apertures
        if regFile != 'none':
            for i in range(apRedX.size):
                tf.write('circle\t{0}\t{1}\t6\n'.format(apRedX[i], apRedY[i]))
        
        
        ## extract the data here
        apertureResults = allApertureResults[useIndicies]
        smallRadius = radius[useIndicies]
        sortedRadius = smallRadius[np.argsort(smallRadius)]
        if verbose:
            print("The largest radii", sortedRadius[-4:])
        
        
        ## do a check to see if there are actually apertures there
        ## or if we are at the edge of the image
        
        #############################################
        ## now calculate the depth!
        ## plot the results
        
            
        ## First clip the data to remove outliers
        ## that will skew the SD
        medianFlux = np.median(apertureResults)
        mad = np.median(abs(apertureResults - medianFlux))
        sigma_mad = 1.4826*mad

        
        
        if (plot != 'none'):
            fig = plt.figure()
            bins = np.arange(-3.0*sigma_mad, 3.0*sigma_mad, sigma_mad/5.0)
            n,bin,patches = plt.hist(apertureResults, bins = bins, facecolor = 'green', alpha = 0.75)
            
            # split by region
            north = (apRedX > np.median(apRedX)) & (apRedY > np.median(apRedY))
            south = (apRedX < np.median(apRedX)) & (apRedY < np.median(apRedY))
            
            # plot the median etc
            plt.plot([medianFlux, medianFlux], [0, max(n)])
            plt.plot([medianFlux-sigma_mad, medianFlux-sigma_mad], [0, max(n)], color = 'k')
            plt.plot([medianFlux+sigma_mad, medianFlux+sigma_mad], [0, max(n)], color = 'k')
            plt.xlim([medianFlux -5.0*sigma_mad, medianFlux + 5.0*sigma_mad])
            
        
        if verbose:
            print("The mad sigma is ", sigma_mad, " after first run, mag = {0:.2f}".format(-2.5*np.log10(5.0*sigma_mad) + zeropoint))

        if sigma_mad > 1E-15:
            
            if sigmaClip > 0.0:
                # First sigma clip
                good_indicies = (apertureResults > medianFlux - sigmaClip*sigma_mad) & \
                    (apertureResults < medianFlux + sigmaClip*sigma_mad)
                
                clippedResults = apertureResults[good_indicies]
                if verbose:
                    print("After the first clipping we have ", np.sum(good_indicies), " good indicies.")
                
        ## Second sigma clip
                medianFlux = np.median(clippedResults)
                mad = np.median(abs(clippedResults - medianFlux))
                sigma_mad = 1.4826*mad
                
                good_indicies_two = (clippedResults > medianFlux - sigmaClip*sigma_mad) & \
                    (clippedResults < medianFlux + sigmaClip*sigma_mad)
                
                if verbose:
                    print("After the second clipping we have ", np.sum(good_indicies_two), " good indicies.")
                    
                finalClippedResults = clippedResults[good_indicies_two]
                if (plot != 'none'):
                    n,bin,patches = plt.hist(finalClippedResults, bins = bins, facecolor = 'red', alpha = 0.75)
                
                    if fitGauss:
                        
                        # fit a gaussian to these points!
                        (mu, sigma) = norm.fit(clippedResults)                
                        yll = mlab.normpdf(bins, mu, sigma)
                        
                        # normalise this
                        plt.plot(bins, np.sum(n)*yll/np.sum(yll), 'k', linewidth=3)
                        print ("The gauss sigma is ", sigma, " after first run, mag = {0:.2f}".format(-2.5*np.log10(5.0*sigma) + zeropoint))

            
            ## Now calculate the final mad sigma value
                medianFlux = np.median(finalClippedResults)
                mad = np.median(abs(finalClippedResults - medianFlux))
                sigma_mad = 1.4826*mad
                if verbose:
                    print("After second clipping... = ", sigma_mad)
                    print("The depth = {0:.2f}".format(return_mag(sigma_mad, zeropoint, sigma = 5.0)))
                    
            
            
        #######################################
        # if not clip then just use this.
        #######################################
        
        if float(zeropoint) > -1.0:
            ## convert to magnitudes
            localDepths[xi] = return_mag(sigma_mad, zeropoint, sigma = 5.0)
            
        else:
            ## keep as fluxes
            localDepths[xi] = sigma_mad


        if np.isnan(localDepths[xi]):
            print('NAN here', sigma_mad)
            #exit()

        if sortedRadius[0] < diffx*0.1:
            ## There are no apertures nearby...
            ## Set a masked value (for plotting)
            maskArray[xi] = 1
            
        if xi % 1000 == 0:
            print("At position ", xi)

        if (plot != 'none'):
            pdf.savefig(fig)
            plt.close()
    
    if regFile != 'none':
        tf.close()

    if plot != 'none':
        pdf.close()
        print('Plot at ', plotname)

    # get rid of nans
    bad_indicies = np.isnan(localDepths)
    if np.any(bad_indicies):
        print('Fixing NANs in the local_depth code')
        localDepths[bad_indicies] = -99.0

    if mask: 
        return localDepths, maskArray
    else:
        print("CONSIDER updating code to use masked depths.")
        return localDepths

def return_mag(flux_counts: float, zeropoint: float, sigma: float = 1.0) -> float:
    ''' 
    Calculate the magnitude from flux counts.

    Parameters:
    flux_counts (float): Flux counts.
    zeropoint (float): The zeropoint value.
    sigma (float, optional): Sigma value. Defaults to 1.0.

    Returns:
    float: The calculated magnitude.
    '''

    if flux_counts < 1E-15:
        return -99.0
    else:
        return -2.5*math.log(sigma*flux_counts, 10.0) + zeropoint
    
def grid_depths(gridTable: dict, x: np.ndarray, y: np.ndarray, faster: bool = True, verbose: bool = False, nearby: bool = False) -> np.ndarray:
    ''' 
    Code to find the closest depth measurement from my previous analysis. Faster than truly local depths

    Parameters:
    gridTable (dict): Dictionary containing grid information.
    x (np.ndarray): Array of x coordinates.
    y (np.ndarray): Array of y coordinates.
    faster (bool, optional): If True, use a faster method. Defaults to True.
    verbose (bool, optional): If True, enable verbose output. Defaults to False.
    nearby (bool, optional): If True, consider nearby pixels. Defaults to False.

    Returns:
    np.ndarray: Array of depth values.
    '''

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

    ## loop through the grid instead of each object
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
        
       
   
def grid_psf(gridTable: dict, x: np.ndarray, y: np.ndarray, faster: bool = True, verbose: bool = False, nearby: bool = False) -> np.ndarray:
    ''' 
    Finds the closest psf measurement from PSFEx map. Adapted from grid_depths, changed a couple things to adapt to PSFEx map.

    Parameters:
    gridTable (dict): Dictionary containing grid information.
    x (np.ndarray): Array of x coordinates.
    y (np.ndarray): Array of y coordinates.
    faster (bool, optional): If True, use a faster method. Defaults to True.
    verbose (bool, optional): If True, enable verbose output. Defaults to False.
    nearby (bool, optional): If True, consider nearby pixels. Defaults to False.

    Returns:
    np.ndarray: Array of PSF values.
    '''

    xgrid = gridTable['x']
    ygrid = gridTable['y']
    keys = gridTable.colnames


    PSFsOverField = gridTable['ef_2.0']

    ## Make an output array
    psfArray = np.zeros(x.size)
    psfArray[:] = -99.0

    if faster:
        if verbose:
            print("Using faster method.")
            print("Input array size is ", x.size)
        deltay = np.min(ygrid)
        deltax = np.min(xgrid)
        
    ## loop through the grid instead of each object
        for xi in range(xgrid.size):
            
            xmin = xgrid[xi] - deltax
            xmax = xgrid[xi] + deltax
            ymin = ygrid[xi] - deltay
            ymax = ygrid[xi] + deltay

            ii = (x > xmin) & (x <= xmax) & (y > ymin) & (y <= ymax)
            
            psfArray[ii] = PSFsOverField[xi]
                
            # use this to average over nearby pixels.
            
                
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
                print("The nearby PSFs are = ", PSFsOverField[mini])
                print("Before = ", PSFsOverField[mini][0])
            
            psfArray[xi] = PSFsOverField[mini][0]


    return psfArray
        
        
   
