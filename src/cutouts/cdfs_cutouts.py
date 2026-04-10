"""
cutout_codes.py

Functions and classes used to make cutouts, for the CDFS field

Created: Thursday 21st March 2024.
Adapted from xmm_cutouts.py: Wednesday 8th April 2026.
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

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

# plt.rcParams['axes.linewidth'] = 2.5
# plt.rcParams.update({'font.size': 15})
# plt.rcParams['figure.dpi'] = 100

ground_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data'
plot_dir = Path.cwd().parent.parent / 'plots' / 'cutouts' / 'z4'
refcat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'
video_dir = Path.home().parents[1] / 'hoy' / 'VIDEO_FINAL'
euclid_base_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'euclid'


def find_VIDEO_tile(ra, dec):
    """Given an RA, Dec, find the VIDEO tile the object lies in."""
    point = Point(ra, dec)

    tiles = {
        'CDFS1': ground_dir / 'CDFS1' / 'HSC-G.fits',
        'CDFS2': ground_dir / 'CDFS2' / 'HSC-G.fits',
        'CDFS3': ground_dir / 'CDFS3' / 'HSC-G.fits',
    }

    for name, path in tiles.items():
        with fits.open(path) as hdu:
            wcs = WCS(hdu[0].header)
            footprint = wcs.calc_footprint()
            polygon = Polygon(footprint)
            if polygon.contains(point):
                return name  # Found the matching tile

    return None  # No match found



def findPlotLimits(data: np.ndarray) -> tuple:

    mean = np.mean(data)
    std_dev = np.std(data)

    # Sigma clip the data
    data = data[(data < mean + 3 * std_dev)]

    # Recalculate mean and std_dev
    mean = np.mean(data)
    std_dev = np.std(data)

    # # Sigma clip the data
    # data = data[(data < mean + 3 * std_dev)]

    # # Recalculate mean and std_dev
    # mean = np.mean(data)
    # std_dev = np.std(data)

    lower = mean - 2 * std_dev
    upper = mean + 5 * std_dev

    return lower, upper



def CDFSCutout(ra: float, dec:float, contained_in: Optional[np.array] = None, size: float = 10.0, 
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

    # Find the VIDEO tile
    tile_name = find_VIDEO_tile(ra, dec)

    #! First get the ground-based cutouts
    vista_dir = ground_dir / tile_name

    ### u ###
    with fits.open(vista_dir / f'{tile_name}_VOICE_u.fits') as hdu_u:
        data_u = hdu_u[0].data
        hdr_u = hdu_u[0].header
        pix_scale = np.abs(hdr_u['CD1_1']) * 3600
        wcs_u = WCS(hdr_u)

        cutout_gru = Cutout2D(data_u, c, size=size/pix_scale, wcs=wcs_u)

    ### g ###
    with fits.open(vista_dir / f'{tile_name}_VOICE_g.fits') as hdu_g:
        data_g = hdu_g[0].data
        hdr_g = hdu_g[0].header
        pix_scale = np.abs(hdr_g['CD1_1']) * 3600
        wcs_g = WCS(hdr_g)

        cutout_grg = Cutout2D(data_g, c, size=size/pix_scale, wcs=wcs_g)
    
    ### r ####
    with fits.open(vista_dir / f'{tile_name}_VOICE_r.fits') as hdu_r:
        data_r = hdu_r[0].data
        hdr_r = hdu_r[0].header
        pix_scale = np.abs(hdr_r['CD1_1']) * 3600
        wcs_r = WCS(hdr_r)

        cutout_grr = Cutout2D(data_r, c, size=size/pix_scale, wcs=wcs_r)

    ### i ####
    with fits.open(vista_dir / f'{tile_name}_VOICE_i.fits') as hdu_i:
        data_i = hdu_i[0].data
        hdr_i = hdu_i[0].header
        pix_scale = np.abs(hdr_i['CD1_1']) * 3600
        wcs_i = WCS(hdr_i)

        cutout_gri = Cutout2D(data_i, c, size=size/pix_scale, wcs=wcs_i)

    ### G ####
    with fits.open(vista_dir / 'HSC-G.fits') as hdu_G:

        data_G = hdu_G[0].data
        hdr_G = hdu_G[0].header
        pix_scale = np.abs(hdr_G['CD1_1']) * 3600
        wcs_G = WCS(hdr_G)

        cutout_grG = Cutout2D(data_G, c, size=size/pix_scale, wcs=wcs_G)

    ### R ###
    with fits.open(vista_dir / 'HSC-R.fits') as hdu_R:

        data_R = hdu_R[0].data
        hdr_R = hdu_R[0].header
        pix_scale = np.abs(hdr_R['CD1_1']) * 3600
        wcs_R = WCS(hdr_R)

        cutout_grR = Cutout2D(data_R, c, size=size/pix_scale, wcs=wcs_R)


    #### I ####
    with fits.open(vista_dir / 'HSC-I.fits') as hdu_I:
            
            data_I = hdu_I[0].data
            hdr_I = hdu_I[0].header
            pix_scale = np.abs(hdr_I['CD1_1']) * 3600
            wcs_I = WCS(hdr_I)
    
            cutout_grI = Cutout2D(data_I, c, size=size/pix_scale, wcs=wcs_I)

    #### Z ####
    with fits.open(vista_dir / 'HSC-Z.fits') as hdu_Z:

        data_Z = hdu_Z[0].data
        hdr_Z = hdu_Z[0].header
        pix_scale = np.abs(hdr_Z['CD1_1']) * 3600
        wcs_Z = WCS(hdr_Z)

        cutout_grZ = Cutout2D(data_Z, c, size=size/pix_scale, wcs=wcs_Z)


    #### Y ####
    with fits.open(video_dir / f'{tile_name.lower()}_Y.fits') as hdu_Y:

        data_Y = hdu_Y[0].data
        hdr_Y = hdu_Y[0].header
        pix_scale = np.abs(hdr_Y['CD1_1']) * 3600
        wcs_Y = WCS(hdr_Y)

        cutout_grY = Cutout2D(data_Y, c, size=size/pix_scale, wcs=wcs_Y)

    #### J ####
    with fits.open(video_dir / f'{tile_name.lower()}_J.fits') as hdu_J:
    #with fits.open(vista_dir / 'UVISTA_J_dr5_rc1.fits') as hdu_J:

        data_J = hdu_J[0].data
        hdr_J = hdu_J[0].header
        pix_scale = np.abs(hdr_J['CD1_1']) * 3600
        wcs_J = WCS(hdr_J)

        cutout_grJ = Cutout2D(data_J, c, size=size/pix_scale, wcs=wcs_J)

    #### H ####
    with fits.open(video_dir / f'{tile_name.lower()}_H.fits') as hdu_H:
    #with fits.open(vista_dir / 'UVISTA_H_dr5_rc1.fits') as hdu_H:
            
        data_H = hdu_H[0].data
        hdr_H = hdu_H[0].header
        pix_scale = np.abs(hdr_H['CD1_1']) * 3600
        wcs_H = WCS(hdr_H)

        cutout_grH = Cutout2D(data_H, c, size=size/pix_scale, wcs=wcs_H)

    #### Ks ####
    with fits.open(video_dir / f'{tile_name.lower()}_Ks.fits') as hdu_K:
            
        data_K = hdu_K[0].data
        hdr_K = hdu_K[0].header
        pix_scale = np.abs(hdr_K['CD1_1']) * 3600
        wcs_K = WCS(hdr_K)

        cutout_grK = Cutout2D(data_K, c, size=size/pix_scale, wcs=wcs_K)

    #! Euclid Q1
    euclid_dir = euclid_base_dir / tile_name

    #### VIS_Q1 ####
    with fits.open(euclid_dir / f'{tile_name}_VIS_MOSAIC.fits') as hdu_vis:

        data_vis = hdu_vis[0].data
        hdr_vis = hdu_vis[0].header
        pix_scale = np.abs(hdr_vis['CD1_1']) * 3600
        wcs_vis = WCS(hdr_vis)

        cutout_vis = Cutout2D(data_vis, c, size=size/pix_scale, wcs=wcs_vis)
    
    ### YE_Q1 ####
    with fits.open(euclid_dir / f'{tile_name}_Y_MOSAIC.fits') as hdu_ye:

        data_ye = hdu_ye[0].data
        hdr_ye = hdu_ye[0].header
        pix_scale = np.abs(hdr_ye['CD1_1']) * 3600
        wcs_ye = WCS(hdr_ye)

        cutout_ye = Cutout2D(data_ye, c, size=size/pix_scale, wcs=wcs_ye)
    
    ### JE_Q1 ####
    with fits.open(euclid_dir / f'{tile_name}_J_MOSAIC.fits') as hdu_je:

        data_je = hdu_je[0].data
        hdr_je = hdu_je[0].header
        pix_scale = np.abs(hdr_je['CD1_1']) * 3600
        wcs_je = WCS(hdr_je)

        cutout_je = Cutout2D(data_je, c, size=size/pix_scale, wcs=wcs_je)

    ### HE_Q1 ####
    with fits.open(euclid_dir / f'{tile_name}_H_MOSAIC.fits') as hdu_he:

        data_he = hdu_he[0].data
        hdr_he = hdu_he[0].header
        pix_scale = np.abs(hdr_he['CD1_1']) * 3600
        wcs_he = WCS(hdr_he)

        cutout_he = Cutout2D(data_he, c, size=size/pix_scale, wcs=wcs_he)

    # # Make a GR stack
    # optical_stack = np.sum([cutout_grG.data, cutout_grr.data], axis=0) / 2.
    # optical_stack = Cutout2D(optical_stack, c, size=size/pix_scale, wcs=cutout_grG.wcs)

    # # Smooth the stack
    # smooth_optical_stack = ndimage.gaussian_filter(optical_stack.data, sigma=2)
    # smooth_optical_stack = Cutout2D(smooth_optical_stack, c, size=size/pix_scale, wcs=cutout_grG.wcs)

    # # Make a yYJHK stack
    # NIR_stack = np.sum([cutout_gry.data, cutout_grY.data, cutout_grJ.data, cutout_grH.data, cutout_grK.data], axis=0) / 4.
    # NIR_stack = Cutout2D(NIR_stack, c, size=size/pix_scale, wcs=cutout_grY.wcs)

    # # Make a spitzer stack
    # spitzer_stack = np.sum([cutout_ch1.data, cutout_ch2.data], axis=0) / 2.
    # spitzer_stack = Cutout2D(spitzer_stack, c, size=size/pix_scale, wcs=cutout_ch1.wcs)

    #! Plot
    # cutouts = [optical_stack, smooth_optical_stack, cutout_grI, cutout_grNB0816, cutout_grZ, cutout_grNB0921, NIR_stack, spitzer_stack, cutout_gry, cutout_grY, cutout_grJ, cutout_grH, cutout_grK, cutout_ch1, cutout_ch2]
    # cutout_titles = ['GR', r'$GR_{\sigma}$', 'HSC-I', 'NB0816', 'HSC-Z', 'NB0921', 'yYJHK', 'ch1+ch2', 'HSC-Y', 'Y', 'J', 'H', 'K', 'ch1', 'ch2']

    # cutouts are u, g, HSC-G, r, HSC-R, i, HSC-I, HSC-Z, VIS_Q1, Y, YE_Q1, J, JE_Q1, H, HE_Q1, Ks
    cutouts = [cutout_gru, cutout_grg, cutout_grG, cutout_grr, cutout_grR, cutout_gri, cutout_grI, cutout_grZ, cutout_vis, cutout_grY, cutout_ye, cutout_grJ, cutout_je, cutout_grH, cutout_he, cutout_grK]
    cutout_titles = ['u', 'g', 'HSC-G', 'r', 'HSC-R', 'i', 'HSC-I', 'HSC-Z', 'IE', 'Y', 'YE', 'J', 'JE', 'H', 'HE', 'Ks']

    fig, ax = plt.subplots(2, 8, figsize=(20, 12))

    # Maximise the size of the subplots
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    for i, cutout in enumerate(cutouts):
            
        vmin, vmax = findPlotLimits(cutout.data)
        ax[i//8, i%8].imshow(cutout.data, origin='lower', cmap='gist_yarg', vmin=vmin, vmax=vmax)
        ax[i//8, i%8].set_title(cutout.wcs.wcs.ctype[0])
        ax[i//8, i%8].set_xticks([])
        ax[i//8, i%8].set_yticks([])
        ax[i//8, i%8].set_title(cutout_titles[i])

        # Remove axis labels
        ax[i//8, i%8].set_xticks([])
        ax[i//8, i%8].set_yticks([])

    # Remove unused subplot (16th subplot)
    fig.delaxes(ax[1, 7])

    plt.tight_layout()

    # Save fig
    if save_cutout:
        plt.savefig(save_dir / f'XMM_ID_{plot_title}.pdf', dpi=300)

    return fig, ax
    #plt.show()
    #plt.close()
    #return None



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

    #! Nathan z=5 XMM sample
    # radio_dir = Path.home() / 'highz_rgs' / 'data' / 'catalogues'
    # t = Table.read(radio_dir / 'combined_z5_matches_nospecz.fits')

    #! Nathan z=4 XMM sample
    radio_dir = Path.home() / 'highz_rgs' / 'data' / 'catalogues'
    t = Table.read(radio_dir / 'Nathan_z4_XMM_Xmatch.fits')

    #t = t[t['RA'] < 140]
    ra = t['RA']
    dec = t['DEC']
    ID = t['UID']


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
        #Cutout(ra[i], dec[i], size=6., add_centre_lines=True)
        XMMCutout(ra[i], dec[i], size=10., plot_title=ID[i], save_cutout=True, save_dir=plot_dir)
        #Cutout(ra[i], dec[i], size=10., plot_title='Big Three Dragons')   
        #Cutout(ra[i], dec[i], size=4., plot_title=ID[i] + ', z=' + str(z[i]) + ', Muv=' + str(Muv[i]))


    
