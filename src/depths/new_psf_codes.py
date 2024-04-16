import numpy as np
import os
from astropy.table import Table,Column, hstack
from pathlib import Path
from typing import Union, List
from astropy.io import fits
import sep
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.backends.backend_pdf
from new_depth_codes import read_image_lis
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit



def psfex(image_name: Path, filter_name: str, field_name: str, zeropoint: float, depth_dir: Path,
          ap_diametersAS: np.ndarray = np.array([0.2, 0.3, 0.5, 1.0, 1.8, 2.0, 3.0]),
          wht_name: str = 'NONE', wht_type: str = 'NONE', seg_name: str = 'NONE',
          output_dir: str = 'none', overwrite: bool = False, overwrite_PSF: bool = False,
          input_sex: Path = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'HSC_SSP_DR3' / 'config_files' / 'euclid_psfex.sex',
          input_PSFEx: Path = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'HSC_SSP_DR3' / 'config_files' / 'default_cube.psfex',
          bzk: bool = False, stars_only: bool = False, use_cat: str = 'NONE') -> None:

    """
    Run PSFEx (Point Spread Function Extractor) on an image.

    Args:
    -----
    - image_name (Path): The path to the input image.
    - filter_name (str): The name of the filter.
    - field_name (str): The name of the field.
    - zeropoint (float): The zeropoint of the image.
    - depth_dir (str): The directory containing depth information.
    - ap_diametersAS (np.ndarray, optional): Array of aperture diameters in arcseconds. Default is [0.2, 0.3, 0.5, 1.0, 1.8, 2.0, 3.0].
    - wht_name (str, optional): The name of the weight file. Default is 'NONE'.
    - wht_type (str, optional): The type of weight. Default is 'NONE'.
    - seg_name (str, optional): The name of the segmentation map file. Default is 'NONE'.
    - output_dir (str, optional): The directory to save output files. Default is 'none'.
    - overwrite (bool, optional): Whether to overwrite existing files. Default is False.
    - overwrite_PSF (bool, optional): Whether to overwrite existing PSF files. Default is False.
    - input_sex (Union[str, Path], optional): Path to the SExtractor configuration file. Default is 'video_psfex.sex'.
    - input_PSFEx (Union[str, Path], optional): Path to the PSFEx configuration file. Default is 'default_cube.psfex'.
    - bzk (bool, optional): Whether to use BZK stars. Default is False.
    - stars_only (bool, optional): Whether to use stars only. Default is False.
    - use_cat (str, optional): The catalog to use. Default is 'NONE'.

    Returns:
    --------
    - None

    PSFEx outputs PSF models, samples, and catalogues to the output_dir. Enclosed fluxes are saved to the enclosedflux directory.

    """

    os.environ['EXTRACTOR_DIR'] = str(Path.home().parent.parent / 'users' / 'videouser' / 'sextractor' / 'share' / 'sextractor')
    plot = True
    
    ##############################################
    psfSize = 75
    psfDegree = 5
    psfNSNAP = 2

    assocRad = 1.0/0.1 # Euclid pix scale

    rstepAS = 0.05 ## arcsec
    rmaxAS = 9.0 ## arcsec
    numRad = int(rmaxAS/rstepAS)
    rAS = np.arange(numRad)*rstepAS + rstepAS
    normRas = 5.0
    
    ##############################################
    
    # main code where all the action happens!
    # first define the output files
    if output_dir == 'none':
        # put everything here!
        # make a new directory
        output_dir = 'psf/'
        if os.path.isdir(output_dir) == False:
            os.system('mkdir ' + output_dir)

    # for the seg map
    #imageDir = output_dir + 'images/'
    #if os.path.isdir(imageDir) == False:
    #    os.system('mkdir ' + imageDir)

    plot_dir = output_dir + '/plots/'
    if os.path.isdir(plot_dir) == False:
        os.system('mkdir ' + plot_dir)

    catDir = output_dir + '/catalogues/'
    if os.path.isdir(catDir) == False:
        os.system('mkdir ' + catDir)

    efDir = output_dir + '/enclosedflux/'
    if os.path.isdir(efDir) == False:
        os.system('mkdir ' + efDir)

    resultsDir = output_dir + '/results/'
    if os.path.isdir(resultsDir) == False:
        os.system('mkdir ' + resultsDir)
        
    # make a sub-directory here?  no but useful
    parts = image_name.split('/')
    if len(parts) > 1:
        baseName = parts[-1]

    # also remove the cat bit
    pparts = baseName.split('.')
    baseName = pparts[0]
    print("The base name is ", baseName)
        
    # check all the necessary files exist
    imyes = os.path.isfile(image_name)
    if imyes == False:
        print("ERROR: the image file does not exist, check the path:\n" + image_name)
        exit()
        
    if wht_name != 'NONE':
        whtyes = os.path.isfile(wht_name)
        if whtyes == False:
            print("ERROR: the weight file does not exist, check the path:\n" + wht_name)            
            exit()

    
    ####################################################
    # Get a star catalogue
    if filter_name == 'NONE':
        filter_name = baseName

    # get the pixel scale
    hdulist = fits.open(image_name)
    imageHeader = hdulist[0].header
    cdone_o = -3600.0*imageHeader['CD1_1']
    pix_scale = round(cdone_o, 5)
    print("The pixel scale is {0:.4f}".format(pix_scale))
    naxis1 = imageHeader['NAXIS1']
    naxis2 = imageHeader['NAXIS2']
    
    rpix = rAS/pix_scale

    #####################################################
    # Define the point sources for the run
    depthCat = depth_dir + '/{0}/catalogues/d{1}.fits'.format(field_name, filter_name)

    # get the five sigma limit
    five_sigFile = depth_dir + '/{0}/results/{1}_300.txt'.format(field_name, filter_name)
    five_sig = Table.read(five_sigFile, format = 'ascii.commented_header')

    print('################' + filter_name + '#####################')

    if filter_name[0] != 'f':
        five_sigHere = five_sig['1.0as'][0]

    if filter_name[0] == 'f':
        five_sigHere = five_sig['0.3as'][0]

    
    if bzk:
        print('BZK stars...')
    else:
        # use the depth catalogue
        print('Using catalogue', depthCat)
        plot_name = plot_dir + '/{0}_sizemag.pdf'.format(filter_name)
        starName = sizemag_stars(depthCat, five_sigHere, filter_name, field_name, pix_scale, plot_dir = plot_dir, star_dir = catDir, plot_name = plot_name, size_type = 'FLUX_RADIUS')

    if stars_only:
        return

    if use_cat != 'NONE':
        starName = starName.replace(filter_name, use_cat)
        print('Using star cat name ', starName)
        
    #############################################################
    # Run a catalogue in the required format for PSFEx
    seCatalogue = catDir  + filter_name + '.fits'
    
    # FITS_LDAC
    keywords = ' -CATALOG_TYPE FITS_LDAC -CATALOG_NAME ' + seCatalogue + \
               ' -MAG_ZEROPOINT '+ str(zeropoint)
    
    keywordsWeight = ' -WEIGHT_TYPE ' + wht_type + \
                     ' -WEIGHT_IMAGE ' +  wht_name
    
    #starCatalogue = 'sky.list'
    associateKeywords = ' -ASSOC_DATA 1,2,3 -ASSOC_NAME ' + starName + \
        ' -ASSOC_PARAMS 4,5 -ASSOC_RADIUS ' + str(assocRad) + \
        ' -ASSOC_TYPE NEAREST' + \
        ' -ASSOCSELEC_TYPE MATCHED'
        
    extraKeywords = associateKeywords
          
    if os.path.isfile(seCatalogue) and overwrite == False:
        print("Catalogue has been created previous.  Skipping this step.")
    else:
        print('/mnt/users/videouser/sextractor/bin/sex '+  image_name +' -c ' + str(input_sex) + extraKeywords + keywords + keywordsWeight )
        os.system('/mnt/users/videouser/sextractor/bin/sex '+  image_name +' -c ' + str(input_sex) + extraKeywords + keywords + keywordsWeight) 
        print("Source extractor catalogue has been saved to ", seCatalogue)
        

        if (filter_name[0:2] == 'ch'): #field_name == 'COSMOS' or
            hdulist = fits.open(seCatalogue)
            datahere = hdulist[2].data
            datahere['FLUX_RADIUS'] = datahere['FLUX_RADIUS']/3.0
            
            hdulist.writeto(seCatalogue, clobber = True)
            print("Editing the se catalogue for my ch results.")
            
    #############################################################
    # Run PSFEx on this catalogue!
    
    ## keywords
    keywords = ' -PSF_DIR ' + resultsDir + ' -PSF_SIZE ' + str(psfSize) + \
               ',' + str(psfSize) + ' -PSFVAR_DEGREES ' + str(psfDegree) + \
               ' -PSFVAR_NSNAP ' + str(psfNSNAP) + ' -CHECKIMAGE_TYPE ' + \
               ' SNAPSHOTS,SAMPLES,RESIDUALS ' + ' -CHECKIMAGE_NAME ' + \
               '{0}snap.fits,{0}samples.fits,{0}resi.fits '.format(resultsDir) + \
               ' -SAMPLE_AUTOSELECT N -PSF_SAMPLING 1.0 -SAMPLE_FWHMRANGE 1,10' #+ \
                #! CHECK PLOT KEYWORDS BELOW HERE.
               #' -CHECKPLOT_TYPE FWHM -CHECKPLOT_DEV PSC' + \
               #' -CHECKPLOT_NAME {0}fwhm'.format(plot_dir) # -BASIS_TYPE GAUSS-LAGUERRE -SAMPLE_FWHMRANGE 2.0,8.0 -VERBOSE_TYPE FULL'.format(plot_dir, filter_name) #,{0}ellipticity,{0}residuals,{0}chis'.format(resultsDir) + \
               #'-SAMPLE_FWHMRANGE 2.0,8.0 ' -PSF_ACCURACY 0.1'#+ ' -CHECKIMAGE_TYPE ' #,ELLIPTICITY,RESIDUALS,CHI2 ' + \
               #' PROTOTYPES,SNAPSHOTS,RESIDUALS ' + ' -CHECKIMAGE_NAME ' + \
               #'proto.fits,snap.fits,resi.fits'
    
    # Define the important psf file
    psfOutFile = resultsDir + 'snap_{0}.fits'.format(filter_name)
    
    # Can uncomment this line if we want to look at some snaps.
    #psfOutFile = resultsDir + 'snap_{0}_rohan.fits'.format(filter_name)
    
    print('psfex '+ seCatalogue +' -c ' + str(input_PSFEx) + keywords)
    overwrite_PSF = True
    
    if os.path.isfile(psfOutFile) and overwrite_PSF == False:
        print("PSFEx has already been run, do not overwrite.")
    else:
        print("Running psfex", seCatalogue)
        os.system('/mnt/zfsusers/varadaraj/psfex/bin/psfex '+ seCatalogue +' -c ' + str(input_PSFEx) + keywords)    
        
    #############################################################
    # Extract the PSFs over the full FOV
    print("Now extracting the PSF snapshots over the field. ")
    hdulist = fits.open(psfOutFile)
    hdulist.info()
    fullPSFdata = hdulist[0].data
    fullPSFheader = hdulist[0].header
    shapePSF = fullPSFdata.shape
    print("The PSF dimensions are ", shapePSF)
    
    psf_samp = fullPSFheader['PSF_SAMP']
    
    halfwidth = (shapePSF[2]-1)/2
    print("The central pixel is (", halfwidth, ',', halfwidth, ')')
    
    stepsizex = naxis1/psfNSNAP
    stepsizey = naxis2/psfNSNAP
    
    print("There are ", psfNSNAP*psfNSNAP, " psfs to extract.")
    
    ## Make an output array
    enclosedFluxArray = np.zeros((len(ap_diametersAS), psfNSNAP*psfNSNAP))
    x = np.zeros(psfNSNAP*psfNSNAP)
    y = np.zeros(psfNSNAP*psfNSNAP)

    for xi in range(psfNSNAP):
        for yi in range(psfNSNAP):
            
            posx = stepsizex*(xi+0.5)
            posy = stepsizey*(yi+0.5)
            index = yi + xi*psfNSNAP
            
            ## This is correct!!  Checked with UltarVISTA - looking at the corners
            ## not covered by HSC
            x[index] = naxis1-posx
            y[index] = naxis2 -posy

            data = fullPSFdata[xi, yi, :, :]
            psfHereFull = data.byteswap().newbyteorder()    
            psfHere = psfHereFull/np.sum(psfHereFull)

            ## plot the results
            if plot:
                if (index > 5) & (index < 7):
                    cogplot = resultsDir + '{0}_COG.pdf'.format(filter_name)
                    fig, ax1 = plt.subplots()    
            
            # sep expects array-like inputs.
            xArray = [halfwidth]
            yArray = [halfwidth]
            xArray = np.array(xArray)
            yArray = np.array(yArray)
            
            ##############################################
            ## Run my enclosed flux code on this cut-out
            
            ## IMPORTANT - I need to scale by PSF_SAMP
            rpixHere = rpix/psf_samp
            flux, fluxerr, flag = sep.sum_circle(psfHere, xArray, yArray, rpixHere, subpix = 5)

            # normalise this
            jj = (rAS == normRas)
            fluxnorm = flux/flux[jj]
            
            ## Plot the results
            if plot:
                if (index > 5) & (index < 7):
                    line1, = plt.plot(rAS, flux, label = 'Global bg', linestyle = '--', color = 'deepskyblue', )


                    plt.plot(rAS, fluxnorm, linestyle = ':', color = 'red')
                    
                    ax1.set_ylim([0.4, 1.1])
                    plt.plot([0, 10], [1.0, 1.0], 'k:', linewidth = 2)
                    
                    plt.xlabel('Aperture Radius (arsec)')
                    plt.title('Enclosed flux for filter ' + filter_name)
                    
                    ## Add an image
                    left, bottom, width, height = [0.55, 0.3, 0.3, 0.3]
                    ax2 = fig.add_axes([left, bottom, width, height])
                    ax2.imshow(psfHere, cmap = 'GnBu_r', origin = 'lower', interpolation = 'none')
                
            ## Find the values of interest
            for ri, diameterRequired in enumerate(ap_diametersAS):

                radiusRequired = diameterRequired/2.0
                cond = (rAS > radiusRequired-0.0001) & (rAS < radiusRequired+0.0001)
                fluxEnclosed = np.extract(cond, flux)
                fluxEnclosedNorm = np.extract(cond, fluxnorm)
                
                #print  "Aperture = ", radiusRequired, " AS, fenc = ", fluxEnclosed
                enclosedFluxArray[ri, index] = fluxEnclosed[0]
                if plot:
                    if (index > 5) & (index < 7):
                        ax1.text(2.0, 0.5+ri*0.05, '{0:.2f} ({2:.2f}), {1:.1f}as'.format(fluxEnclosed[0], radiusRequired, fluxEnclosedNorm[0]))
                
            if plot:
                if (index > 5) & (index < 7):
                    fig.savefig(cogplot)
                    print(cogplot)
                               
    ############################################################
    ## Tabulate these results, so they can be extracted
    ## at will e.g. in my flux_errors code
    localTable = Table([x, y], names = ['x', 'y'], dtype = ['f4', 'f4'])
    
    ##############################################################
    ## Make plots of the enclosed flux and check these make sense
    ## Create pdf with multiple pages
    plot_name = plot_dir + filter_name + '.pdf'
    pdf = matplotlib.backends.backend_pdf.PdfPages(plot_name)
    
    outputTable = efDir + filter_name + '.fits'

    # make a txt file of the median ef values
    medianFile = efDir + filter_name + '_peak.txt'
    efmedian = np.zeros(ap_diametersAS.size)
    
    for ri, diameterRequired in enumerate(ap_diametersAS):
        
    ## plot the results
        fig = plt.figure(figsize=(6,8))
        gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
        plt.subplot(gs[0])
        plt.axis('equal')
    ## get colours
        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(x, y, s = 30, c = enclosedFluxArray[ri, :], cmap = cm, \
                             linewidth = 0.0) #, vmin = magmin, vmax = magmax)
        plt.colorbar(sc)
        apDstring = "{0:.1f}".format(diameterRequired)
        plt.title('Enclosed flux for filter = ' +filter_name + " apD = " + apDstring)
        
        ## now make a histogram to go underneath
        ax = plt.subplot(gs[1])
        binwidth = 0.01
        efhere = enclosedFluxArray[ri,:]
        median = np.median(efhere)
        mad = 1.58*np.median(np.abs(efhere-median))

        minef = median - 3.0*mad
        maxef = median + 3.0*mad

        if (maxef - minef) < 0.05:
            minef = minef-0.02
            maxef = maxef+ 0.02

        print('Min = {0}, max = {1}'.format(minef, maxef))
        bins = np.arange(minef+binwidth/2.0, maxef, binwidth)
        histy, histx, _ = plt.hist(enclosedFluxArray[ri, :], facecolor= 'k', alpha = 0.8, bins = bins,  density = True)
        
        ## calculate the median
        smallbins = bins[:-1]
        print(histy)
        mm = (histy == np.max(histy))

        enclosedFlux = smallbins[mm]+binwidth/2.0

        efmedian[ri] = enclosedFlux[0]
        
        plt.plot([enclosedFlux, enclosedFlux], ax.get_ylim())
        
        pdf.savefig(fig)
        
        ## add the enclosed flux values to the table
        newcolumn = Column(name = 'ef_' + apDstring, data = enclosedFluxArray[ri, :], dtype = 'f4')
        localTable.add_column(newcolumn)     
        
    pdf.close()
    plt.close()
    print("Plot saved to ", plot_name)

    # Write out the median
    newtb = Table([ap_diametersAS, efmedian], names = ['apD', 'ef'])
    newtb['ef'].format = '7.2f'
    newtb.write(medianFile, overwrite = True, format = 'ascii.commented_header')
    print(medianFile)
    
    # Write out table
    localTable.write(outputTable, overwrite=True)
    print("Table saved to ", outputTable)

    return

def get_psf(field_name: str, req_filters: List[str], queue: str = 'none', ap_diametersAS: np.ndarray = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), 
            data_dir: Path = Path.home() / 'euclid', output_dir: str = 'none', overwrite: bool = False, stars_only: bool = False, use_cat: str = 'NONE') -> None:
    """
    Run PSF extraction on images for a given field and set of filters.

    Args:
    - field_name (str): The name of the field.
    - req_filters (List[str]): List of required filters.
    - queue (str, optional): The queue to use for spawning. Default is 'none'.
    - ap_diametersAS (np.ndarray, optional): Array of aperture diameters in arcseconds. Default is [1.8, 2.0, 3.0, 4.0, 5.0].
    - data_dir (str or Path, optional): The directory containing image data. Default is Path.home() / 'euclid'.
    - output_dir (str, optional): The directory to save output files. Default is 'none'.
    - overwrite (bool, optional): Whether to overwrite existing files. Default is False.
    - stars_only (bool, optional): Whether to process stars only. Default is False.
    - use_cat (str, optional): The catalog to use. Default is 'NONE'.

    Returns:
    - None
    This function runs PSFEx on the queue.
    """

    #! loop through each Euclid filter
    for fi, filter_name in enumerate(req_filters):

        # Read in the images file
        dirHere = data_dir / filter_name / field_name
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
            imageDir = imagedata['directory'][j]

            depth_dir = str(Path.cwd().parent.parent / 'data' / 'depths')

            
            if imageDir == 'here':
                imageDir = data_dir / filter_name / field_name
            
            # Now spawn the depths!
            if queue == 'none':
                print("Running here ")
                
                psfex(imageDir+'/'+image_name, tile_name, field_name, zeropoint, depth_dir, wht_name = imageDir+'/'+wht_name, wht_type = wht_type, output_dir = output_dir, overwrite = overwrite, stars_only = stars_only, overwrite_PSF = False, use_cat= use_cat)

            else:
                print("Spawning in the queue...", queue)
                # make shell script
                tmpName = "tmp_{1}_{0}.sh".format(tile_name, field_name)
                f = open(tmpName, 'w')
                f.write('#!/bin/bash\n')
                f.write('python3 stupid_psf.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}'.format(imageDir+'/'+image_name, imageDir+'/'+wht_name, wht_type, zeropoint, output_dir, tile_name, overwrite, stars_only, field_name, depth_dir))
                f.close()
                
                # now execute this
                command = "addqueue -c 'tmp_{0}' -m 8 -q {0} -d ./{1}".format(queue, tmpName)
                #print(command)
                os.system('chmod +x {0}'.format(tmpName))
                os.system(command)
        
    return

def sizemag_stars(input_cat: str, five_sig: float, filter_name: str, field: str, pix_scale: float, 
                  plot_dir: str = 'plots/', star_dir: str = '', size_type: str = 'FWHM_IMAGE', extra: str = '', plot_name: str = 'sizemag.pdf') -> str:
    """
    Process a catalog of stars to extract information about their size and magnitude.

    Args:
    - input_cat (str): Path to the input catalog.
    - five_sig (float): Five sigma detection limit.
    - filter_name (str): Name of the filter.
    - field (str): Name of the field.
    - pix_scale (float): Pixel scale in arcseconds/pixel.
    - plot_dir (str, optional): Directory to save plots. Default is 'plots/'.
    - star_dir (str, optional): Directory to save extracted star data. Default is ''.
    - size_type (str, optional): Type of size measurement. Default is 'FWHM_IMAGE'.
    - extra (str, optional): Extra identifier for output filenames. Default is ''.
    - plot_name (str, optional): Name of the output plot. Default is 'sizemag.pdf'.

    Returns:
    - str: Path to the output file containing extracted star data.
    """
    

    # read in table
    tbdata = Table.read(input_cat)

    classstarcut = 0.95

    ps = (tbdata['CLASS_STAR'] > classstarcut)
    pssize = np.array(tbdata[size_type][ps])
    expectedSize = np.median(pssize)
    print('Expected size is ', expectedSize*pix_scale)

    #minSize = 0.1*expectedSize*pix_scale     #! IF FILTER IS VIS, SET LOWER MIN SIZE
    minSize = 0.5*expectedSize*pix_scale 
    maxSize = 1.5*expectedSize*pix_scale # arcseconds
    
    apMag = tbdata['MAG_APER'][:,0]
    apMag = tbdata['MAG_AUTO']
    faintMag = five_sig - 2.0 # 1.5
    brightMag = faintMag -7.5
    
    bright = (apMag < faintMag) & (apMag > brightMag) & (tbdata[size_type]*pix_scale < maxSize) & (tbdata[size_type]*pix_scale > minSize) # cut at 20 sigma, and a reasonable size /minsize
    newtb = tbdata[bright]
    print('Cutting catalogue at five sig = {0:.1f}, {1:.2f} < m < {2:.2f}, {3:.2f} < fwhm < {4:.2f}, {5}'.format(five_sig, faintMag, brightMag, minSize, maxSize, size_type))
    print("There are {0}/{1} positions left after inital cuts.".format(len(newtb), len(tbdata)))
    
    ##############################################################
    # Now make a sanity plot
    pdf = matplotlib.backends.backend_pdf.PdfPages(plot_name)
    
    apMagR = apMag[bright]
    size = newtb[size_type]*pix_scale
    
    fig = plt.figure()
    plt.hexbin(apMagR, pix_scale*newtb[size_type], bins = 'log')
    
    plt.xlabel('mAB/mag')
    plt.ylabel('{0} (/arcsec)'.format(size_type))
    
    pdf.savefig()
    
    fig = plt.figure()
    plt.scatter(apMagR, pix_scale*newtb[size_type],s=0.2, rasterized = True)
    
    plt.xlabel('mAB/mag')
    plt.ylabel('{0} (/arcsec)'.format(size_type))
    
    
    #############################################
    # Try and fit iteratively
    # with minimal intervention.
    
    # read in the cut from a nice file
    cuts = Table.read('./stars/star_param_{0}.txt'.format(field), format = 'ascii.commented_header')
    ff = (cuts['name'] == filter_name)
    brMag = cuts['bright'][ff]
    faMag = cuts['faint'][ff]
    maxsize = cuts['upperfwhm'][ff]
    satLevel = cuts['saturate'][ff]
    nosigma = 2.0
    nosigmaup = 1.5
    faintLevel = cuts['faintest'][ff]
    
    # first get a rough estimate from the brighter stars
    br = (apMagR < faMag) & (apMagR > brMag)  & (size < maxsize)
    firstmedian = np.median(size[br])
    mad = np.median(np.abs(size[br]-firstmedian))
    firstsigma = 1.58*mad

    # do a sigma cut here
    sc = (apMagR < faMag) & (apMagR > brMag)  & (size < firstmedian + 2.0*firstsigma) & (size > firstmedian - 2.0*firstsigma)
    
    # and remeasure
    estmedian = np.median(size[sc])
    mad = np.median(np.abs(size[sc]-firstmedian))
    estsigma = 1.58*mad
    
    #fwhmcut = estlocus + 1.58*mad
    xline = np.arange(brightMag, faintMag, 0.1)
    yline = xline*0.0 + estmedian
    plt.plot(xline, yline, color= 'k', label = 'median of bright stars')
    
    plt.plot(xline, yline+estsigma, color= 'green', label = 'one sigma')
    plt.plot(xline, yline-estsigma, color= 'green')
    plt.legend()
    
    #################################################
    # now fit these points to get a line
    keep = (apMagR < faintLevel) & (apMagR > satLevel) & (size < estmedian+nosigmaup*estsigma) & (size > estmedian-nosigma*estsigma)
    
    x = apMagR[keep]
    y = pix_scale*newtb[size_type]
    y = y[keep]
    
    A,B = curve_fit(line, x, y)[0]
    
    plt.scatter(x, y, color = 'green', s= 0.2, rasterized = True)
    
    # plot the result!
    yline = xline*A + B
    plt.plot(xline, yline, color= 'red') #, label = 'fit to bright stars')
    
    # cut according to this fittted line
    midpoint = apMagR*A + B
    
    tofit = (apMagR < faintLevel) & (apMagR > satLevel) & (size < midpoint+nosigmaup*estsigma) & (size > midpoint-nosigma*estsigma) & (size < maxsize)

    # this is the final data
    
    locusdata = newtb[tofit]
    plt.scatter(apMagR[tofit], size[tofit], color = 'red', s= 0.2)
    
    pdf.savefig()
    
    # look at the cross section
    fig = plt.figure()
    ax = plt.subplot()
    plt.hist(size[br], 50, label = 'Bright')
    plt.hist(size[keep], 50, label = 'Keep')
    bottom, top = ax.get_ylim()
    plt.xlim([0.9*estmedian, 1.1*estmedian])
    plt.plot([estmedian, estmedian], [0, top], color = 'k')
    plt.plot([estmedian+estsigma, estmedian+estsigma], [0, top], color = 'green')
    plt.legend()
    pdf.savefig()
    
    pdf.close()
    
    # great now return the coordinates or whatever it needs
    tablename = star_dir + '{0}_stars{1}.fits'.format(filter_name, extra)
    locusdata.rename_column('NUMBER', 'ID')
    locusdata.rename_column('ALPHA_J2000', 'RA')
    locusdata.rename_column('DELTA_J2000', 'DEC')
    
    locusdata.write(tablename, overwrite=True) #, format = 'ascii')
    print("Table saved to ", tablename, " with {0} objects.".format(len(locusdata)))
    
    smalldata = locusdata[['ID', 'RA', 'DEC', 'X_IMAGE', 'Y_IMAGE', 'CLASS_STAR', 'FLUX_RADIUS']]
    starName = star_dir + '{0}_stars{1}.ascii'.format(filter_name, extra)
    smalldata.write(starName, overwrite=True, format = 'ascii')
    print("Stars written to ", starName)

    print("Inspect ", plot_name, " and adjust parameters.")
    
    return starName



def line(x, A, B):
    return A*x + B