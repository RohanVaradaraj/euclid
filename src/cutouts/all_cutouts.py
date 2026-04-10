"""
cutout_codes.py

Functions and classes used to make cutouts

Created: Thursday 21st March 2024.
"""

from shapely.geometry import Point, Polygon
import numpy as np
from pathlib import Path
from astropy.io import ascii
from typing import Optional
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
import glob
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from astropy.table import Table
import astropy.units as u
from scipy import ndimage
from astropy.nddata.utils import NoOverlapError

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

# plt.rcParams['axes.linewidth'] = 2.5
# plt.rcParams.update({'font.size': 15})
# plt.rcParams['figure.dpi'] = 100

ground_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data'
euclid_dir = Path.home() / 'euclid'
cweb_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'CWEB'
primer_dir = Path.home() / 'JWST'
plot_dir = Path.cwd().parent.parent / 'plots' / 'cutouts'
refcat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'


def readAllCOSMOSGalaxies() -> Table:

    """
    Read in all crossmatched objects found from NED from the text file to an astropy table.
    """

    # Read in the file
    cosmos_file = refcat_dir / 'all_COSMOS_highz.txt'
    cosmos = np.loadtxt(cosmos_file, dtype=str, delimiter='|', skiprows=21)

    # Convert to astropy table
    cosmos_colnames = cosmos[0]
    cosmos_data = cosmos[1:]
    cosmos_table = Table(cosmos_data, names=cosmos_colnames)

    cosmos_table.write(refcat_dir / 'all_COSMOS_highz.fits', format='fits', overwrite=True)




def findPlotLimits(data: np.ndarray) -> tuple:

    mean = np.mean(data)
    std_dev = np.std(data)

    # Sigma clip the data
    data = data[(data < mean + 3 * std_dev)]

    # Recalculate mean and std_dev
    mean = np.mean(data)
    std_dev = np.std(data)

    lower = mean - 2 * std_dev
    upper = mean + 5 * std_dev

    return lower, upper


def AllCutout(ra: float, dec:float, contained_in: Optional[np.array] = None, size: float = 10.0, 
           save_cutout: bool = True, save_dir: Path = Path.cwd().parent.parent / 'data' / 'cutouts',
           plot_title: Optional[str] = None,
           add_centre_lines: Optional[bool] = False) -> None:

    """
    Create cutouts from Euclid, ground-based and JWST imaging.

    Parameters
    ----------

    ra : float
        Right ascension of the cutout center in degrees.
    dec : float
        Declination of the cutout center in degrees.
    contained_in : np.array (str, str, str)
        Array with three strings. The strings are the values returned by isCoordInSurveyFootprints.
        The first string is the Euclid tile label (0 if absent), the second is the CWEB tile label (0 if absent) and the third is '1' if in PRIMER.
    size : float, optional
        Size of the cutout in arcseconds. The default is 10.0
    save_cutout : bool, optional
        Whether to save the cutout. The default is False.
    save_dir : Path, optional
        Directory to save the cutout. The default is Path.cwd().parent.parent / 'data' / 'cutouts'.
    plot_title: str, optional
        Can pass a title to the stamp plot. The default is None.
    add_centre_lines: bool, optional
        Whether to add lines to the plot showing the centre of the cutout. The default is False.

    
    Returns
    -------
    None

    """
    
    # Convert RA, DEC to skycoords.
    c = SkyCoord(ra, dec, unit='deg')

    #! First get the ground-based cutouts
    vista_dir = ground_dir / 'COSMOS'
    euclid_dir = ground_dir / 'euclid' / 'images'

    ### G ####
    with fits.open(vista_dir / 'HSC-G_DR3.fits') as hdu_G:

        data_G = hdu_G[0].data
        hdr_G = hdu_G[0].header
        pix_scale = np.abs(hdr_G['CD1_1']) * 3600
        wcs_G = WCS(hdr_G)

        cutout_grG = Cutout2D(data_G, c, size=size/pix_scale, wcs=wcs_G)

    ### R ###
    with fits.open(vista_dir / 'HSC-R_DR3.fits') as hdu_R:

        data_R = hdu_R[0].data
        hdr_R = hdu_R[0].header
        pix_scale = np.abs(hdr_R['CD1_1']) * 3600
        wcs_R = WCS(hdr_R)

        cutout_grR = Cutout2D(data_R, c, size=size/pix_scale, wcs=wcs_R)


    #### I ####
    with fits.open(vista_dir / 'HSC-I_DR3.fits') as hdu_I:
            
            data_I = hdu_I[0].data
            hdr_I = hdu_I[0].header
            pix_scale = np.abs(hdr_I['CD1_1']) * 3600
            wcs_I = WCS(hdr_I)
    
            cutout_grI = Cutout2D(data_I, c, size=size/pix_scale, wcs=wcs_I)

    #### NB0816 ####
    with fits.open(vista_dir / 'HSC-NB0816_DR3.fits') as hdu_NB0816:

        data_NB0816 = hdu_NB0816[0].data
        hdr_NB0816 = hdu_NB0816[0].header
        pix_scale = np.abs(hdr_NB0816['CD1_1']) * 3600
        wcs_NB0816 = WCS(hdr_NB0816)

        cutout_grNB0816 = Cutout2D(data_NB0816, c, size=size/pix_scale, wcs=wcs_NB0816)

    #### Z ####
    with fits.open(vista_dir / 'HSC-Z_DR3.fits') as hdu_Z:

        data_Z = hdu_Z[0].data
        hdr_Z = hdu_Z[0].header
        pix_scale = np.abs(hdr_Z['CD1_1']) * 3600
        wcs_Z = WCS(hdr_Z)

        cutout_grZ = Cutout2D(data_Z, c, size=size/pix_scale, wcs=wcs_Z)

    #### NB0921 ####
    with fits.open(vista_dir / 'HSC-NB0921_DR3.fits') as hdu_NB0921:

        data_NB0921 = hdu_NB0921[0].data
        hdr_NB0921 = hdu_NB0921[0].header
        pix_scale = np.abs(hdr_NB0921['CD1_1']) * 3600
        wcs_NB0921 = WCS(hdr_NB0921)

        cutout_grNB0921 = Cutout2D(data_NB0921, c, size=size/pix_scale, wcs=wcs_NB0921)


    #### y ####
    with fits.open(vista_dir / 'HSC-Y_DR3.fits') as hdu_y:

        data_y = hdu_y[0].data
        hdr_y = hdu_y[0].header
        pix_scale = np.abs(hdr_y['CD1_1']) * 3600
        wcs_y = WCS(hdr_y)

        cutout_gry = Cutout2D(data_y, c, size=size/pix_scale, wcs=wcs_y)  

    #### Y ####
    with fits.open(vista_dir / 'UVISTA_Y_DR6.fits') as hdu_Y:
    #with fits.open(vista_dir / 'UVISTA_Y_dr5_rc1.fits') as hdu_Y:

        data_Y = hdu_Y[0].data
        hdr_Y = hdu_Y[0].header
        pix_scale = np.abs(hdr_Y['CD1_1']) * 3600
        wcs_Y = WCS(hdr_Y)

        cutout_grY = Cutout2D(data_Y, c, size=size/pix_scale, wcs=wcs_Y)

    #### J ####
    with fits.open(vista_dir / 'UVISTA_J_DR6.fits') as hdu_J:
    #with fits.open(vista_dir / 'UVISTA_J_dr5_rc1.fits') as hdu_J:

        data_J = hdu_J[0].data
        hdr_J = hdu_J[0].header
        pix_scale = np.abs(hdr_J['CD1_1']) * 3600
        wcs_J = WCS(hdr_J)

        cutout_grJ = Cutout2D(data_J, c, size=size/pix_scale, wcs=wcs_J)

    #### H ####
    with fits.open(vista_dir / 'UVISTA_H_DR6.fits') as hdu_H:
    #with fits.open(vista_dir / 'UVISTA_H_dr5_rc1.fits') as hdu_H:
            
        data_H = hdu_H[0].data
        hdr_H = hdu_H[0].header
        pix_scale = np.abs(hdr_H['CD1_1']) * 3600
        wcs_H = WCS(hdr_H)

        cutout_grH = Cutout2D(data_H, c, size=size/pix_scale, wcs=wcs_H)

    #### Ks ####
    with fits.open(vista_dir / 'UVISTA_K_DR6.fits') as hdu_K:
            
        data_K = hdu_K[0].data
        hdr_K = hdu_K[0].header
        pix_scale = np.abs(hdr_K['CD1_1']) * 3600
        wcs_K = WCS(hdr_K)

        cutout_grK = Cutout2D(data_K, c, size=size/pix_scale, wcs=wcs_K)

    ### VIS ###
    with fits.open(euclid_dir / 'COSMOS_VIS_MOSAIC.fits') as hdu_VIS:
                
        data_VIS = hdu_VIS[0].data
        hdr_VIS = hdu_VIS[0].header
        pix_scale = np.abs(hdr_VIS['CD1_1']) * 3600
        wcs_VIS = WCS(hdr_VIS)

        # Get cutout, but catch no overlap errors
        try:
            cutout_VIS = Cutout2D(data_VIS, c, size=size/pix_scale, wcs=wcs_VIS)

        except NoOverlapError:
            print('No overlap with VIS data')
            # Create empty data array as placeholder for when there is no overlap with euclid image.
            empty_data = np.zeros((int(size/pix_scale), int(size/pix_scale)))
            cutout_VIS = Cutout2D(empty_data, (empty_data.shape[0]//2, empty_data.shape[1]//2), 
                                size=size/pix_scale, wcs=None) 
    
    ### Ye ###
    with fits.open(euclid_dir / 'COSMOS_Y_MOSAIC.fits') as hdu_Ye:

        data_Ye = hdu_Ye[0].data
        hdr_Ye = hdu_Ye[0].header
        pix_scale = np.abs(hdr_Ye['CD1_1']) * 3600
        wcs_Ye = WCS(hdr_Ye)

        try:
            cutout_Ye = Cutout2D(data_Ye, c, size=size/pix_scale, wcs=wcs_Ye)
        except NoOverlapError:
            # Create empty data array as placeholder for when there is no overlap with euclid image.
            empty_data = np.zeros((int(size/pix_scale), int(size/pix_scale)))
            cutout_Ye = Cutout2D(empty_data, (empty_data.shape[0]//2, empty_data.shape[1]//2), 
                                size=size/pix_scale, wcs=None) 
    


    ### Je ###
    with fits.open(euclid_dir / 'COSMOS_J_MOSAIC.fits') as hdu_Je:

        data_Je = hdu_Je[0].data
        hdr_Je = hdu_Je[0].header
        pix_scale = np.abs(hdr_Je['CD1_1']) * 3600
        wcs_Je = WCS(hdr_Je)

        try:
            cutout_Je = Cutout2D(data_Je, c, size=size/pix_scale, wcs=wcs_Je)
        except NoOverlapError:
            # Create empty data array as placeholder for when there is no overlap with euclid image.
            empty_data = np.zeros((int(size/pix_scale), int(size/pix_scale)))
            cutout_Je = Cutout2D(empty_data, (empty_data.shape[0]//2, empty_data.shape[1]//2), 
                                size=size/pix_scale, wcs=None)


    ### He ###
    with fits.open(euclid_dir / 'COSMOS_H_MOSAIC.fits') as hdu_He:

        data_He = hdu_He[0].data
        hdr_He = hdu_He[0].header
        pix_scale = np.abs(hdr_He['CD1_1']) * 3600
        wcs_He = WCS(hdr_He)

        try:
            cutout_He = Cutout2D(data_He, c, size=size/pix_scale, wcs=wcs_He)
        except NoOverlapError:
            # Create empty data array as placeholder for when there is no overlap with euclid image.
            empty_data = np.zeros((int(size/pix_scale), int(size/pix_scale)))
            cutout_He = Cutout2D(empty_data, (empty_data.shape[0]//2, empty_data.shape[1]//2), 
                                size=size/pix_scale, wcs=None)

    #### IRAC ch1 ####
    with fits.open(vista_dir / 'COSMOS_ch1_COMPLETE_microJy.fits') as hdu_ch1:
            
        data_ch1 = hdu_ch1[0].data
        hdr_ch1 = hdu_ch1[0].header
        pix_scale = np.abs(hdr_ch1['CD1_1']) * 3600
        wcs_ch1 = WCS(hdr_ch1)

        cutout_ch1 = Cutout2D(data_ch1, c, size=size/pix_scale, wcs=wcs_ch1)

    #### IRAC ch2 ####
    with fits.open(vista_dir / 'COSMOS_ch2_COMPLETE_microJy.fits') as hdu_ch2:

        data_ch2 = hdu_ch2[0].data
        hdr_ch2 = hdu_ch2[0].header
        pix_scale = np.abs(hdr_ch2['CD1_1']) * 3600
        wcs_ch2 = WCS(hdr_ch2)

        cutout_ch2 = Cutout2D(data_ch2, c, size=size/pix_scale, wcs=wcs_ch2)

    # Make a GRI stack
    optical_stack = np.sum([cutout_grG.data, cutout_grR.data, cutout_grI.data], axis=0) / 3.
    optical_stack = Cutout2D(optical_stack, c, size=size/pix_scale, wcs=cutout_grG.wcs)

    # Smooth the stack
    smooth_optical_stack = ndimage.gaussian_filter(optical_stack.data, sigma=2)
    smooth_optical_stack = Cutout2D(smooth_optical_stack, c, size=size/pix_scale, wcs=cutout_grG.wcs)

    # Make a YJHK stack
    NIR_stack = np.sum([cutout_grY.data, cutout_grJ.data, cutout_grH.data, cutout_grK.data], axis=0) / 4.
    NIR_stack = Cutout2D(NIR_stack, c, size=size/pix_scale, wcs=cutout_grY.wcs)

    # Make a spitzer stack
    spitzer_stack = np.sum([cutout_ch1.data, cutout_ch2.data], axis=0) / 2.
    spitzer_stack = Cutout2D(spitzer_stack, c, size=size/pix_scale, wcs=cutout_ch1.wcs)

    #! Plot
    cutout_titles = [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$', 'ch1', 'ch2', 'NB0816', 'NB0921', 'HSC-GRI', 'Y', 'J', 'H', 'K', r'GRI$_{\sigma}$', 'HSC-Z', 'HSC-Y']
    cutouts = [cutout_VIS, cutout_Ye, cutout_Je, cutout_He, cutout_ch1, cutout_ch2, cutout_grNB0816, cutout_grNB0921, optical_stack, cutout_grY, cutout_grJ, cutout_grH, cutout_grK, smooth_optical_stack, cutout_grZ, cutout_gry]

    #! saving fits files for John
    # cutouts = [cutout_VIS, cutout_Ye, cutout_Je, cutout_He]
    # cutout_name_appends = ['VIS', 'Y', 'J', 'H']
    # cutout_dir = Path.cwd().parents[1] / 'data' / 'cutouts' / 'bowler_sample'
    # for i, cutout in enumerate(cutouts):

    #     # Save the cutout
    #     hdu_cutout = fits.PrimaryHDU(data=cutout.data)
    #     hdu_cutout.header.update(cutout.wcs.to_header())
    #     hdu_cutout.writeto(cutout_dir / f'{plot_title}_{cutout_name_appends[i]}.fits', overwrite=True)


    #! PLOTTING: UNCOMMENT
    fig, ax = plt.subplots(2, 8, figsize=(20, 12))

    # Maximise the size of the subplots
    plt.subplots_adjust(wspace=0.01, hspace=0.01)

    for i, cutout in enumerate(cutouts):
            
        vmin, vmax = findPlotLimits(cutout.data)
        ax[i//8, i%8].imshow(cutout.data, origin='lower', cmap='gist_yarg', vmin=vmin, vmax=vmax)
        #ax[i//8, i%8].set_title(cutout.wcs.wcs.ctype[0])
        ax[i//8, i%8].set_xticks([])
        ax[i//8, i%8].set_yticks([])
        ax[i//8, i%8].set_title(cutout_titles[i])

        # If the title has 'ch', draw a red circle of diameter 2.8 arcsec at the centre
        # if 'ch' in cutout_titles[i]:
        #     circle = plt.Circle((cutout.data.shape[0]//2, cutout.data.shape[1]//2), 1.4/pix_scale, color='red', fill=False, lw=2)
        #     ax[i//8, i%8].add_artist(circle)

        # Remove axis labels
        ax[i//8, i%8].set_xticks([])
        ax[i//8, i%8].set_yticks([])

    plt.tight_layout()
    plot_dir = Path.cwd().parent.parent / 'plots' / 'cutouts'
    #plt.savefig(plot_dir / f'LAE_stamps_10arcsec.pdf')
    plt.show()

    return fig, ax
    # plt.show()
    # plt.close()
    # return None



if __name__ == '__main__':
    
    #! PRIMER STARS
    t = Table.read(Path.cwd().parents[3] / 'data' / 'psf' / 'COSMOS' / 'catalogues' / 'f410m_stars.fits')

    ra = t['ALPHA_J2000']
    dec = t['DELTA_J2000']
    ID = t['NUMBER']

    #! REBELS sources
    # t = ascii.read(Path.cwd().parent.parent / 'data' / 'mosaic' / 'REBELS.csv', format='csv')
    # t = t[t['RA'] > 148]

    # # # Skip first few
    # # #t = t[:]

    # ra = t['RA']
    # dec = t['Dec']
    # z = t['Redshift (z)']
    # ID = t['Object Name']
    # ID = [name.split('>')[1].split('<')[0] for name in ID]

    #! Strong lens
    #ra = [150.00280406167596]
    #dec = [2.2002751856804608]

    #! Nathan's z=3 sources
    # nathan_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'HSC_SSP_DR3' / 'ref_catalogues' / 'nathan'
    # cat = Table.read(nathan_dir / 'Z3_FinalSample.fits')
    # cat.sort('MUV', reverse=False)
    # cat = cat[cat['RA'] > 148]
    # ra = cat['RA']
    # dec = cat['DEC']
    # muv = cat['MUV']
    # z = cat['ZPhot']
    # ID = cat['UID']

    #! DEVILS sources
    # devils_dir = Path.home() / 'DEVILS' / 'dr1cats' / 'data' / 'catalogues'
    # cat = Table.read(devils_dir / 'D10VisualMorphology.csv', format='csv')
    # print(cat.colnames)
    # cat = cat[(cat['RAcen'] > 148) & (cat['DECcen'] > 1.8)]
    # # Drop 'NA' strings from zBest column
    # cat = cat[cat['zBest'] != 'NA']
    # # Convert zBest to float
    # cat['zBest'] = cat['zBest'].astype(float)
    # # Remove zBest values of -99
    # cat = cat[cat['zBest'] > 0.]

    # cat = cat[cat['FIRST_CLASS'] != 'NA']
    # cat.sort('zBest')
    # ra = cat['RAcen']
    # dec = cat['DECcen']

    #! Harikane z=12-16 sources
    # hd1 = '10:01:51.31 02:32:50.0'
    # hd2 = '02:18:52.44 -05:08:36.1'

    # # Use skycoord to convert these coordinates to degrees
    # c1 = SkyCoord(hd1, unit=(u.hourangle, u.deg))
    # c2 = SkyCoord(hd2, unit=(u.hourangle, u.deg))

    # ra = [c1.ra.deg, c2.ra.deg]
    # dec = [c1.dec.deg, c2.dec.deg]

    # ID = ['HD1', 'HD2']

    #! Rebecca's z>8.5 sources
    # uvista_1212 = '10:02:31.81 02:31:17.10'
    # uvista_237 = '10:00:31.88 01:57:50.04'

    # c1 = SkyCoord(uvista_1212, unit=(u.hourangle, u.deg))
    # c2 = SkyCoord(uvista_237, unit=(u.hourangle, u.deg))

    # ra = [c1.ra.deg, c2.ra.deg]
    # dec = [c1.dec.deg, c2.dec.deg]

    # ID = ['uvista_1212', 'uvista_237']

    #! Big three dragons
    # b1 = '10:01:40.69 01:54:52.42'
    # b1 = SkyCoord(b1, unit=(u.hourangle, u.deg))
    # ra = [b1.ra.deg]
    # dec = [b1.dec.deg]

    #! Bad astrometry sources
    # stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
    # stars = ascii.read(stars_dir / f'H_outside_pixscale_vista_euclid_coords.ascii')
    # ra = stars['RA_euclid']
    # dec = stars['DEC_euclid']

    #! Normal stars
    # stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
    # stars = ascii.read(stars_dir / f'Y_vista_euclid_coords.ascii')
    # ra = stars['RA_euclid']
    # dec = stars['DEC_euclid']

    #! Couple bad JWST astrometry sources
    # stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
    # stars = ascii.read(stars_dir / f'H_outside_pixscale_jwst_euclid_coords.ascii')
    # ra = stars['RA_euclid']
    # dec = stars['DEC_euclid'] 

    #! Spectra from Vincent!
    #ra = [150.21552279971348, 149.98832200622317, 150.26771198215147, 150.24988047861532]
    #dec = [2.02579740958649, 2.2070238197331493, 2.6213922569606107, 2.4663986583043647]

    #ID = [1502155319020257969, 1499883241022070278, 1502677115026213801, 1502498770024664099]

    #! Casey CWEB z>10 sources

    # Super bright 10<z<12 sources
    # ra = ['10:01:26.00', '09:58:55.21', '09:59:59.91', '09:59:49.04']
    # dec = ['01:55:59.70', '02:07:16.77', '02:06:59.90', '01:53:26.19']
    # z = [10.27, 12.54, 11.92, 12.03]
    # Muv = [-21.53, -22.19, -21.89, -21.58]
    # ID = ['COS-z10-1', 'COS-z12-1', 'COS-z12-2', 'COS-z12-3']

    # Bright 10<z<12 sources
    # ra = ['09:59:51.77', '09:59:57.50', '10:00:37.96', '09:59:52.53', '10:01:34.80']
    # dec = ['2:07:15.02', '02:06:20.06', '01:49:32.43', '02:00:23.53', '02:05:41.48']
    # z = [10.06, 10.17, 11.07, 11.70, 11.50]
    # Muv = [-20.62, -20.97, -20.85, -21.13, -20.77]
    # ID = ['COS-z10-2', 'COS-z10-3', 'COS-z11-1', 'COS-z11-2', 'COS-z11-3']

    # # z>13
    # ra = ['09:59:05.75', '10:00:04.24', '10:01:31.17']
    # dec = ['02:04:04.39', '02:02:11.19', '01:58:45.00']
    # z = [13.10, 12.50, 14.7]
    # Muv = [-21.27, -21.03, -20.75]
    # ID = ['COS-z13-1', 'COS-z13-2', 'COS-z14-1']

    # Probable contaminants
    # ra = ['09:59:30.49', '09:59:31.30', '10:00:20.38']
    # dec = ['02:14:44.10', '02:08:33.85', '01:49:58.33']
    # z = [12.63, 13.8, 14.7]
    # Muv = [-21.90, -20.97, -21.32]
    # ID = ['COS-z12-4', 'COS-z13-3', 'COS-z14-2']

    # Convert ra, dec to degrees
    # ra = [SkyCoord(r, d, unit=(u.hourangle, u.deg)).ra.deg for r, d in zip(ra, dec)]
    # dec = [SkyCoord(r, d, unit=(u.hourangle, u.deg)).dec.deg for r, d in zip(ra, dec)]

    #! CR7
    # cr7 = '10:00:58.005 01:48:15.251'
    # cr7 = SkyCoord(cr7, unit=(u.hourangle, u.deg))

    # ra = [cr7.ra.deg]
    # dec = [cr7.dec.deg]
    # ID = ['CR7']

    #! All COSMOS crossmatched sources
    # t = Table.read(refcat_dir / 'all_COSMOS_highz.fits')

    # # Sort by decreasing redshift
    # t.sort('Redshift')

    # # Flip
    # t = t[::-1]


    # ra = t['RA']
    # dec = t['DEC']

    # z = t['Redshift']
    # ID = t['Object Name']

    #! Initial det_YJH candidates
    # t = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'XMATCH_COSMOS_5sig_Ye_2sig_VISTA_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits')
    # ra = t['RA_1']
    # dec = t['DEC_1']
    # ID = t['Object Name']
    # z = t['Redshift']

    # # Sort by flux in Y, reversed
    # t.sort('flux_Ye', reverse=True)
    # print(t['flux_Ye'])

    #! My LAE candidate
    # ra = [150.11833152095758]
    # dec = [2.2522416552619746]
    # ID = [178396]

    #! Stars
    # t = Table.read(Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues' / 'df444w_locus_stars.fits')

    # # Sort by class_star
    # t.sort('CLASS_STAR', reverse=True)
    # print(len(t))

    # t = t[t['CLASS_STAR'] < 0.99]

    # STARS = (t['FLAGS'] < 2) & (t['ELONGATION'] < 1.5) & (t['FWHM_IMAGE'] < 6) & (t['FWHM_IMAGE'] > 5) & (t['CLASS_STAR'] < 0.99) & (t['CLASS_STAR'] > 0.8)
    # t = t[STARS]

    # print(len(t))

    # #t = t[(t['FWHM_IMAGE'] > 5.5) & (t['FWHM_IMAGE'] < 7)]

    # ra = t['RA']
    # dec = t['DEC']
    # class_star = t['CLASS_STAR']
    # flag = t['FLAGS']
    # elong = t['ELONGATION']
    # fwhm = t['FWHM_IMAGE']

    #! Rebecca's sample
    # t = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates' / 'UVISTA_REBELS_sizes.fits')

    # IDs = t['ID_prev']
    # IDs = [ID.replace(' ', '_') for ID in IDs]
    # print(IDs)

    # ra = t['RA']
    # dec = t['DEC']

    #! Brown dwarfs - any proper motion?
    t = Table.read('/mnt/vardy/vardygroupshare/rohan/euclid/data/catalogues/candidates/COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_really_good_BDs_INTERLOPERS_2025_08_22_with_euclid.fits')

    t = t[t['ID'] == 829044]

    ra = t['RA']
    dec = t['DEC']
    ID = t['ID']


    #! Stavrox z=6.8
    ra = [150.1200752011961]
    dec = [2.0898737362171156]

    #! Charlotte z=7 candidate
    ra = [150.6525914016524]
    dec = [2.7949620730349847]


    #! Nathan's z=5 candidates in COSMOS
    # ra = [ 150.0646487736565, 150.2970171829096, 150.0015577303966, 149.93031878361805, 149.55790768581574, 149.81976533364988, 150.67691341136822, 150.37025906264196, 149.9179524404775, 150.41129219687292, 149.6785624899251, 150.02803178397048, 150.4029315800968, 150.31595510191195]
    # dec = [2.1240832246701857, 2.488110301538341, 2.6952885516554907, 1.7686795080791051, 1.790581039771246, 1.9773914886177946, 2.22463751680825, 2.3805154948508487, 2.392064687286876, 2.4150960438461673, 2.5690418786852987, 2.6098005600763656, 2.760439230110506, 2.3366573607029992]


    ############! GET CUTOUTS ############
    for i in range(len(ra)):

        print(f'Object number {i+1} of {len(ra)}')

    #     print(f'Object number {i+1} of {len(ra)}')

    #     if isCoordInCWEB(ra[i], dec[i])[0] == '0':
    #         print('Not in CWEB')
    #         continue

        #print(cat[i]['zBest'])
        #print(cat[i]['FIRST_CLASS'])
        #print(ID[i])
        #print(muv[i])

        #Cutout(ra[i], dec[i], size=10., plot_title=str(ID[i]) + ', z=' + str(z[i]), save_cutout=False)
        #print(class_star[i], flag[i], elong[i], fwhm[i])
        #Cutout(ra[i], dec[i], size=6., save_cutout=False)
        AllCutout(ra[i], dec[i], size=10.) #, add_centre_lines=True)
        #AllCutout(ra[i], dec[i], size=10., plot_title=ID[i], add_centre_lines=True)
        #Cutout(ra[i], dec[i], size=10., plot_title='Big Three Dragons')   
        #Cutout(ra[i], dec[i], size=4., plot_title=ID[i] + ', z=' + str(z[i]) + ', Muv=' + str(Muv[i]))


    
