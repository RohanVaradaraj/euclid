"""
Plot the Euclid VIS sizes of z=6 candidates in COSMOS and CDFS.
Want to see if BDs in CDFS can easily be identified from their VIS morphology.

Created: Friday 10th April 2026.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
import numpy as np
import glob
from astropy.wcs import WCS
import sep
from astropy.nddata import Cutout2D
from matplotlib.patches import Ellipse
from astropy.modeling import models, fitting
from matplotlib.ticker import FixedLocator, FixedFormatter
from scipy.ndimage import gaussian_filter
from scipy.stats import ks_2samp

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['figure.dpi'] = 100

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 3, 'ytick.major.width': 3,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 2, 'ytick.minor.width': 2,
})

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

def fit_gaussian2d(stamp: np.ndarray, box: int = 21):
    """
    Fit a rotated 2D Gaussian to the brightest source in a postage stamp.

    Returns
    -------
    fitted : Gaussian2D
    subdata : 2D array
    xx, yy : meshgrid for the subdata
    x0, y0 : fitted center in subdata coordinates
    """
    data = stamp.astype(float)

    # Simple local background estimate
    bkg = np.nanmedian(data)
    data = data - bkg

    # Find brightest pixel as initial guess
    y_peak, x_peak = np.unravel_index(np.nanargmax(data), data.shape)

    # Define fitting box around the peak
    half = box // 2
    y1 = max(0, y_peak - half)
    y2 = min(data.shape[0], y_peak + half + 1)
    x1 = max(0, x_peak - half)
    x2 = min(data.shape[1], x_peak + half + 1)

    sub = data[y1:y2, x1:x2]

    yy, xx = np.mgrid[y1:y2, x1:x2]

    amp0 = np.nanmax(sub)
    g0 = models.Gaussian2D(
        amplitude=amp0,
        x_mean=x_peak,
        y_mean=y_peak,
        x_stddev=2.0,
        y_stddev=2.0,
        theta=0.0,
    )

    # Keep the fit sane
    g0.amplitude.min = 0
    g0.x_stddev.min = 0.3
    g0.y_stddev.min = 0.3

    fitter = fitting.LevMarLSQFitter()

    try:
        fitted = fitter(g0, xx, yy, sub)
    except Exception as e:
        print(f"Gaussian fit failed: {e}")
        return None, sub, xx, yy, x_peak, y_peak

    return fitted, sub, xx, yy, x_peak, y_peak

def density_contours(ax, x, y, bins=40, levels=(0.68, 0.95), smooth=1.0,
                     color='k', linewidths=1.5, label=None):
    x = np.asarray(x)
    y = np.asarray(y)

    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]

    H, xedges, yedges = np.histogram2d(x, y, bins=bins)
    H = gaussian_filter(H, smooth)

    # Convert histogram counts to cumulative probability levels
    Hflat = H.flatten()
    idx = np.argsort(Hflat)[::-1]
    Hsorted = Hflat[idx]
    cumsum = np.cumsum(Hsorted)
    cumsum /= cumsum[-1]

    # Find contour thresholds corresponding to enclosed probabilities
    thresholds = []
    for lev in levels:
        thresholds.append(Hsorted[np.searchsorted(cumsum, lev)])
    thresholds = sorted(thresholds)

    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    cs = ax.contour(X, Y, H.T, levels=thresholds, colors=color,
                    linewidths=linewidths)
    if label is not None:
        cs.collections[0].set_label(label)
    return cs

######################################! Save cdfs VIS footprints? ######################################
save_footprints = False

#? measure_sizes?
measure_sizes = False

#? If measure_sizes, which fields?
cdfs = True
cosmos = True

#? Whilst measuring sizes, plot cutouts with Gaussian fit overlay?
plot_cutouts = False

######################################! MAKE PLOTS? ######################################
#* RUNS IF MEASURE_FOOTPRINTS = False, ASSUMES THE FWHM MEASUREMENTS ARE ALREADY DONE AND SAVED IN THE TABLES.

#? PLOT 1: Histogram of each FWHM measurement for COSMOS and CDFS candidates
plot_1 = True

#? Plot 2: Scatter plot of FWHM vs Muv for COSMOS and CDFS candidates
plot_2 = False

#? Plot 3: FWHM vs photometric redshift for CDFS and COSMOS
plot_3 = False

#! Plotting dir
plot_dir = Path.cwd().parents[1] / 'plots' / 'brown_dwarfs'

if measure_sizes:
    #! Paths to the data
    cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'

    cosmos_cat = Table.read(cat_dir / 'COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06.fits')
    cdfs_cat = Table.read(cat_dir / 'CDFS_5sig_HSC_Z_nonDet_HSC_G_nonDet_r_candidates_2026_04_08_with_euclid.fits')

    # If there are not already columns for the VIS FWHM, add them to the tables
    if 'VIS_FWHM_PIX' not in cosmos_cat.colnames:
        cosmos_cat['VIS_FWHM_PIX_MOM'] = np.full(len(cosmos_cat), np.nan)
        cosmos_cat['VIS_FWHM_AS_MOM'] = np.full(len(cosmos_cat), np.nan)
        cosmos_cat['VIS_FWHM_PIX_GAUSS'] = np.full(len(cosmos_cat), np.nan)
        cosmos_cat['VIS_FWHM_AS_GAUSS'] = np.full(len(cosmos_cat), np.nan)
    if 'VIS_FWHM_PIX' not in cdfs_cat.colnames:
        cdfs_cat['VIS_FWHM_PIX_MOM'] = np.full(len(cdfs_cat), np.nan)
        cdfs_cat['VIS_FWHM_AS_MOM'] = np.full(len(cdfs_cat), np.nan)
        cdfs_cat['VIS_FWHM_PIX_GAUSS'] = np.full(len(cdfs_cat), np.nan)
        cdfs_cat['VIS_FWHM_AS_GAUSS'] = np.full(len(cdfs_cat), np.nan)

    #! Location of original COSMOS VIS image
    cosmos_vis_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'euclid' / 'DR1' / 'COSMOS'
    cosmos_vis = 'COSMOS_VIS_DR1.fits'
    cosmis_vis_wht = 'COSMOS_VIS_DR1_WHT.fits'

    #! Location of original CDFS VIS images
    cdfs_vis_dir = Path.home() / 'euclid' / 'CDFS'

    if save_footprints:
        #! Loop through all *MOSAIC-VIS-TILE* fits files in the CDFS dir and save the footprints to a .npy file.
        #! Then will check for the CDFS candidates, which VIS tile it lies in, open this, and measure size.
        footprints = []
        for i, fits_file in enumerate(sorted(glob.glob(str(cdfs_vis_dir / '*MOSAIC-VIS_TILE*.fits')))):
            print('File {}/{}: {}'.format(i+1, len(glob.glob(str(cdfs_vis_dir / '*MOSAIC-VIS_TILE*.fits'))), fits_file))
            with fits.open(fits_file) as hdul:
                header = hdul[0].header
                wcs = WCS(header)
                footprints.append(wcs.calc_footprint())
        footprints = np.array(footprints)
        print(footprints)
        np.save('CDFS_VIS_tile_footprints.npy', footprints)
    else:
        footprints = np.load('CDFS_VIS_tile_footprints.npy')

    #! Loop through CDFS cands, find which VIS tile they lie in, show 5x5" cutout
    cdfs_moment_fwhm_pix = []
    cdfs_moment_fwhm_as = []

    if cdfs:
        for i, row in enumerate(cdfs_cat):
            print('Candidate {}/{}'.format(i+1, len(cdfs_cat)))
            ra = row['RA']
            dec = row['DEC']
            for footprint in footprints:
                if (footprint[:,0].min() < ra < footprint[:,0].max()) and (footprint[:,1].min() < dec < footprint[:,1].max()):

                    #! Open the corresponding VIS tile
                    file_name = sorted(glob.glob(str(cdfs_vis_dir / '*MOSAIC-VIS_TILE*.fits')))[footprints.tolist().index(footprint.tolist())]
                    print(file_name)
                    hdu = fits.open(str(cdfs_vis_dir / file_name))[0]
                    wcs = WCS(hdu.header)
                    pix_scale = np.abs(hdu.header['CD1_1'] * 3600)  # arcsec/pixel

                    #! Get pixel coordinates of the candidate
                    x, y = wcs.world_to_pixel_values(ra, dec)

                    #! Make a cutout
                    cutout_size = 3.5 / pix_scale
                    cutout = Cutout2D(hdu.data, (x, y), cutout_size, wcs=wcs)

                    #! Get corresponding weight map
                    tile_id = file_name.split('TILE')[-1].split('-')[0]
                    wht_name = glob.glob(str(cdfs_vis_dir / '*VIS-MAP_WEIGHT_TILE{}*.fits'.format(tile_id)))[0]

                    hdu_wht = fits.open(wht_name)[0]
                    cutout_wht = Cutout2D(hdu_wht.data, (x, y), cutout_size, wcs=wcs)

                    # Cutout weight is inverse variance, convert to rms
                    cutout_wht.data = np.sqrt(1 / cutout_wht.data)

                    #! Source extract with sep
                    # Make cutout data sep friendly
                    cutout.data = cutout.data.astype(np.float32)
                    objects = sep.extract(cutout.data, 2, err=cutout_wht.data)

                    #! Fit a Source Extractor like Gaussian to the candidate and get the FWHM.
                    #? (A) moment-based  

                    # Pick the object closest to the candidate position
                    yc, xc = np.array(cutout.data.shape) / 2.0

                    if len(objects) == 0:
                        print("No objects detected")
                        continue

                    dist2 = (objects['x'] - xc)**2 + (objects['y'] - yc)**2
                    j = np.argmin(dist2)
                    obj = objects[j]

                    # SExtractor-like FWHM in pixels
                    fwhm_pix = 2.0 * np.sqrt(np.log(2.0) * (obj['a']**2 + obj['b']**2))

                    # Convert to arcsec
                    fwhm_arcsec = fwhm_pix * pix_scale

                    print(f"FWHM ~ {fwhm_pix:.2f} pix = {fwhm_arcsec:.3f} arcsec")       

                    #? (B) 2D Gaussian fit
                    fitted, sub, xx, yy, x_peak, y_peak = fit_gaussian2d(cutout.data, box=21)

                    if plot_cutouts:
                        fig, ax = plt.subplots(figsize=(5, 5))
                        vmin, vmax = findPlotLimits(cutout.data)
                        ax.imshow(cutout.data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)

                        # SEP ellipse for the matched object
                        e = Ellipse(
                            xy=(objects['x'][j], objects['y'][j]),
                            width=6 * objects['a'][j],
                            height=6 * objects['b'][j],
                            angle=objects['theta'][j] * 180.0 / np.pi,
                            facecolor='none',
                            edgecolor='red',
                            lw=1.5,
                        )
                        ax.add_patch(e)

                    # Gaussian fit overlay
                    if fitted is not None:
                        model_img = fitted(xx, yy)
                        if plot_cutouts:
                            ax.contour(xx, yy, model_img, levels=5, colors='cyan', linewidths=1)

                        sig_x = float(fitted.x_stddev.value)
                        sig_y = float(fitted.y_stddev.value)
                        theta = float(fitted.theta.value)

                        fwhm_x = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sig_x
                        fwhm_y = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sig_y
                        fwhm_circ = 2.0 * np.sqrt(2.0 * np.log(2.0)) * np.sqrt(sig_x * sig_y)

                        # print(f"Gaussian fit center: ({fitted.x_mean.value:.2f}, {fitted.y_mean.value:.2f})")
                        # print(f"Gaussian sigmas: {sig_x:.2f}, {sig_y:.2f} pix")
                        # print(f"Gaussian FWHM major/minor: {max(fwhm_x, fwhm_y):.2f} / {min(fwhm_x, fwhm_y):.2f} pix")
                        print(f"Gaussian circularized FWHM: {fwhm_circ:.2f} pix = {fwhm_circ * pix_scale:.3f} arcsec")

                        # Print ratio/difference of moment-based and Gaussian FWHMs
                        print(f"FWHM ratio (Gaussian / SExtractor): {fwhm_circ / fwhm_pix:.2f}")

                        #! Save the Gaussian fit FWHM in pixels and arcsec to the list
                        cdfs_moment_fwhm_pix.append(fwhm_circ)
                        cdfs_moment_fwhm_as.append(fwhm_circ * pix_scale)

                        if plot_cutouts:
                            gauss_ellipse = Ellipse(
                                xy=(fitted.x_mean.value, fitted.y_mean.value),
                                width=fwhm_x,
                                height=fwhm_y,
                                angle=np.degrees(theta),
                                facecolor='none',
                                edgecolor='cyan',
                                lw=1.5,
                                linestyle='--',
                            )
                            ax.add_patch(gauss_ellipse)

                    if plot_cutouts:
                        ax.set_title(f"Candidate {i+1}")
                        plt.show()
                        plt.close()

                    #! Save the FWHM measurements to the table
                    cdfs_cat['VIS_FWHM_PIX_MOM'][i] = fwhm_pix
                    cdfs_cat['VIS_FWHM_AS_MOM'][i] = fwhm_arcsec
                    if fitted is not None:
                        cdfs_cat['VIS_FWHM_PIX_GAUSS'][i] = fwhm_circ
                        cdfs_cat['VIS_FWHM_AS_GAUSS'][i] = fwhm_circ * pix_scale

                    #! Break if there are multiple footprints, just want to check the first one for now
                    break

    #! Write the updated tables with the FWHM measurements to new .fits files
    if cdfs:
        cdfs_cat.write('CDFS_5sig_HSC_Z_nonDet_HSC_G_nonDet_r_candidates_2026_04_08_with_euclid_VIS_sizes.fits', overwrite=True)


    #! Loop through COSMOS cands and do the same thing
    cosmos_moment_fwhm_pix = []
    cosmos_moment_fwhm_as = []

    with fits.open(cosmos_vis_dir / cosmos_vis) as hdu:
        wcs = WCS(hdu[0].header)
        pix_scale = np.abs(hdu[0].header['CD1_1'] * 3600)  # arcsec/pixel
        image = hdu[0].data
        footprint = wcs.calc_footprint()

    with fits.open(cosmos_vis_dir / cosmis_vis_wht) as hdu_wht:
        wcs_wht = WCS(hdu_wht[0].header)
        weight_map = hdu_wht[0].data

    for i, row in enumerate(cosmos_cat):
        print('Candidate {}/{}'.format(i+1, len(cosmos_cat)))
        ra = row['RA']
        dec = row['DEC']


        #! Get pixel coordinates of the candidate
        x, y = wcs.world_to_pixel_values(ra, dec)

        #! Check if the candidate lies within the VIS footprint
        if not (footprint[:,0].min() < ra < footprint[:,0].max()) or not (footprint[:,1].min() < dec < footprint[:,1].max()):
            print("Candidate is outside the VIS footprint")
            continue

        #! Make a cutout
        cutout_size = 3.5 / pix_scale
        cutout = Cutout2D(image, (x, y), cutout_size, wcs=wcs)
        cutout_wht = Cutout2D(weight_map, (x, y), cutout_size, wcs=wcs_wht)


        #! Source extract with sep
        # Make cutout data sep friendly
        cutout.data = np.ascontiguousarray(cutout.data, dtype=np.float32)
        cutout_wht.data = np.ascontiguousarray(cutout_wht.data, dtype=np.float32)
        # Cutout weight is inverse variance, convert to rms
        cutout_wht.data = np.sqrt(1 / cutout_wht.data)

        objects = sep.extract(cutout.data, 2, err=cutout_wht.data)

        #! Fit a Source Extractor like Gaussian to the candidate and get the FWHM.
        #? (A) moment-based  

        # Pick the object closest to the candidate position
        yc, xc = np.array(cutout.data.shape) / 2.0

        if len(objects) == 0:
            print("No objects detected")
            continue

        dist2 = (objects['x'] - xc)**2 + (objects['y'] - yc)**2
        j = np.argmin(dist2)
        obj = objects[j]

        # SExtractor-like FWHM in pixels
        fwhm_pix = 2.0 * np.sqrt(np.log(2.0) * (obj['a']**2 + obj['b']**2))

        # Convert to arcsec
        fwhm_arcsec = fwhm_pix * pix_scale

        print(f"FWHM ~ {fwhm_pix:.2f} pix = {fwhm_arcsec:.3f} arcsec")       

        #? (B) 2D Gaussian fit
        fitted, sub, xx, yy, x_peak, y_peak = fit_gaussian2d(cutout.data, box=21)

        if plot_cutouts:
            fig, ax = plt.subplots(figsize=(5, 5))
            vmin, vmax = findPlotLimits(cutout.data)
            ax.imshow(cutout.data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)

            # SEP ellipse for the matched object
            e = Ellipse(
                xy=(objects['x'][j], objects['y'][j]),
                width=6 * objects['a'][j],
                height=6 * objects['b'][j],
                angle=objects['theta'][j] * 180.0 / np.pi,
                facecolor='none',
                edgecolor='red',
                lw=1.5,
            )
            ax.add_patch(e)

        # Gaussian fit overlay
        if fitted is not None:
            model_img = fitted(xx, yy)
            if plot_cutouts:
                ax.contour(xx, yy, model_img, levels=5, colors='cyan', linewidths=1)

            sig_x = float(fitted.x_stddev.value)
            sig_y = float(fitted.y_stddev.value)
            theta = float(fitted.theta.value)

            fwhm_x = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sig_x
            fwhm_y = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sig_y
            fwhm_circ = 2.0 * np.sqrt(2.0 * np.log(2.0)) * np.sqrt(sig_x * sig_y)

            # print(f"Gaussian fit center: ({fitted.x_mean.value:.2f}, {fitted.y_mean.value:.2f})")
            # print(f"Gaussian sigmas: {sig_x:.2f}, {sig_y:.2f} pix")
            # print(f"Gaussian FWHM major/minor: {max(fwhm_x, fwhm_y):.2f} / {min(fwhm_x, fwhm_y):.2f} pix")
            print(f"Gaussian circularized FWHM: {fwhm_circ:.2f} pix = {fwhm_circ * pix_scale:.3f} arcsec")

            # Print ratio/difference of moment-based and Gaussian FWHMs
            print(f"FWHM ratio (Gaussian / SExtractor): {fwhm_circ / fwhm_pix:.2f}")

            #! Save the Gaussian fit FWHM in pixels and arcsec to the list
            cdfs_moment_fwhm_pix.append(fwhm_circ)
            cdfs_moment_fwhm_as.append(fwhm_circ * pix_scale)

            #! Plot
            if plot_cutouts:
                gauss_ellipse = Ellipse(
                    xy=(fitted.x_mean.value, fitted.y_mean.value),
                    width=fwhm_x,
                    height=fwhm_y,
                    angle=np.degrees(theta),
                    facecolor='none',
                    edgecolor='cyan',
                    lw=1.5,
                    linestyle='--',
                )
                ax.add_patch(gauss_ellipse)

        if plot_cutouts:
            ax.set_title(f"Candidate {i+1}")
            plt.show()
            plt.close()

        #! Save the FWHM measurements to the table
        cosmos_cat['VIS_FWHM_PIX_MOM'][i] = fwhm_pix
        cosmos_cat['VIS_FWHM_AS_MOM'][i] = fwhm_arcsec
        if fitted is not None:
            cosmos_cat['VIS_FWHM_PIX_GAUSS'][i] = fwhm_circ
            cosmos_cat['VIS_FWHM_AS_GAUSS'][i] = fwhm_circ * pix_scale

    #! Write the updated tables with the FWHM measurements to new .fits files
    cosmos_cat.write('COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06_with_euclid_VIS_sizes.fits', overwrite=True)

if measure_sizes == False:

    t_cosmos = Table.read('COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06_with_euclid_VIS_sizes.fits')
    t_cdfs = Table.read('CDFS_5sig_HSC_Z_nonDet_HSC_G_nonDet_r_candidates_2026_04_08_with_euclid_VIS_sizes.fits')

    # Restrict to things with a ratio of Gaussian FWHM to moment FWHM between 0.5 and 2, to remove bad fits
    mask_cosmos = (t_cosmos['VIS_FWHM_PIX_GAUSS'] / t_cosmos['VIS_FWHM_PIX_MOM'] > 0.5) & (t_cosmos['VIS_FWHM_PIX_GAUSS'] / t_cosmos['VIS_FWHM_PIX_MOM'] < 2)
    mask_cdfs = (t_cdfs['VIS_FWHM_PIX_GAUSS'] / t_cdfs['VIS_FWHM_PIX_MOM'] > 0.5) & (t_cdfs['VIS_FWHM_PIX_GAUSS'] / t_cdfs['VIS_FWHM_PIX_MOM'] < 2)

    t_cosmos = t_cosmos[mask_cosmos]
    t_cdfs = t_cdfs[mask_cdfs]

    #? Restrict t_cdfs to Zphot > 5.5
    t_cdfs = t_cdfs[t_cdfs['Zphot'] > 5.5]

    #? Make a copy of the cosmos table restricted to z=5.7-6.3
    t_cosmos_restricted = t_cosmos[(t_cosmos['Zphot'] > 5.7) & (t_cosmos['Zphot'] < 6.3)]

    #? And the corresponding sources outside this range
    t_cosmos_outside = t_cosmos[(t_cosmos['Zphot'] <= 5.7) | (t_cosmos['Zphot'] >= 6.3)]

    if plot_1:
        #! PLOT 1: Histogram of each FWHM measurement for COSMOS and CDFS candidates
        
        #? Plot CDFS vs COSMOS, or COSMOS restricted to z=5.7-6.3 vs outside this range?
        restrict_cosmos = False

        plt.figure(figsize=(10, 6))
        if restrict_cosmos:
            #? restricted cosmos
            plt.hist(t_cosmos_restricted['VIS_FWHM_AS_GAUSS'], bins=np.arange(0.0, 0.7, 0.01), alpha=0.8, label='COSMOS (z=5.7-6.3)', histtype='step', color='tab:blue', linewidth=2, zorder=1)

            #? Cosmos outside z=5.7-6.3
            plt.hist(t_cosmos_outside['VIS_FWHM_AS_GAUSS'], bins=np.arange(0.0, 0.7, 0.01), alpha=0.8, label='COSMOS (z<5.7, z>6.3)', histtype='step', color='cyan', linewidth=2, zorder=1, linestyle='--')

            #? Parent sample
            plt.hist(t_cosmos['VIS_FWHM_AS_GAUSS'], bins=np.arange(0.0, 0.7, 0.01), alpha=0.8, label='COSMOS full', histtype='step', color='gray', linewidth=1.5, zorder=0, linestyle=':')

        else:
            plt.hist(t_cosmos['VIS_FWHM_AS_GAUSS'], bins=np.arange(0.0, 0.7, 0.01), alpha=0.8, label='COSMOS (z=5.5-6.5)', histtype='step', color='tab:blue', linewidth=2, zorder=1)
            # plt.hist(t_cosmos['VIS_FWHM_AS_MOM'], bins=np.arange(0.01, 1, 0.01), alpha=0.8, label='COSMOS Moment FWHM', histtype='step', color='cyan', linewidth=2)

            plt.hist(t_cdfs['VIS_FWHM_AS_GAUSS'], bins=np.arange(0.0, 0.7, 0.01), alpha=0.8, label='CDFS', histtype='step', color='tab:red', linewidth=2, zorder=2)
            # plt.hist(t_cdfs['VIS_FWHM_AS_MOM'], bins=np.arange(0.01, 1, 0.01), alpha=0.8, label='CDFS Moment FWHM', histtype='step', color='orange', linewidth=2)



        plt.xlabel('FWHM (arcsec)')
        plt.ylabel('Number of candidates')

        # Vertical line at VIS PSF FWHM of ~0.2"
        plt.axvline(0.208, color='black', linestyle='--', label='VIS PSF FWHM', linewidth=2, alpha=0.7, zorder=0)

        plt.xlim(0., 1)

        plt.legend(loc='upper right')

        ax = plt.gca()

        # ax.set_xscale('log')
        ax.set_xlim(0.1, 0.7)

        ticks = [0.1, 0.2, 0.5, 1.0]
        ticks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FixedFormatter(['0.1', '0.2', '0.5', '1']))
        ax.xaxis.set_major_formatter(FixedFormatter(['{:.1f}'.format(tick) for tick in ticks]))

        # plt.title('Euclid VIS FWHM measurements for z=6 candidates')

        # Optional: remove minor ticks (they clutter log plots)
        ax.xaxis.set_minor_locator(plt.NullLocator())

        plt.tight_layout()

        # Remove NaNs just in case
        cosmos_fwhm = t_cosmos['VIS_FWHM_AS_GAUSS']
        cdfs_fwhm = t_cdfs['VIS_FWHM_AS_GAUSS']
        cosmos_restricted_fwhm = t_cosmos_restricted['VIS_FWHM_AS_GAUSS']
        cosmos_outside_fwhm = t_cosmos_outside['VIS_FWHM_AS_GAUSS']

        mask = np.isfinite(cosmos_fwhm)
        cosmos_fwhm = cosmos_fwhm[mask]

        mask = np.isfinite(cdfs_fwhm)
        cdfs_fwhm = cdfs_fwhm[mask]

        mask = np.isfinite(cosmos_restricted_fwhm)
        cosmos_restricted_fwhm = cosmos_restricted_fwhm[mask]

        mask = np.isfinite(cosmos_outside_fwhm)
        cosmos_outside_fwhm = cosmos_outside_fwhm[mask]

        # KS test
        if not restrict_cosmos:
            ks_stat, p_value = ks_2samp(cosmos_fwhm, cdfs_fwhm)

        if restrict_cosmos:
            # KS test for cosmos vs restricted cosmos
            ks_stat, p_value = ks_2samp(cosmos_restricted_fwhm, cosmos_outside_fwhm)

        print(f"KS statistic: {ks_stat:.4f}")
        print(f"p-value: {p_value:.4e}")
        
        # Save the plot
        if not restrict_cosmos:
            plt.savefig(plot_dir / 'candidates_z6_vis_fwhm_histogram.pdf')
        if restrict_cosmos:
            plt.savefig(plot_dir / 'candidates_cos_z6_zrange_restricted_vis_fwhm_histogram.pdf', bbox_inches='tight')
        
        plt.show()

    if plot_2:

        #! Plot 2: Scatter plot of FWHM vs Muv for COSMOS and CDFS candidates
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.scatter(t_cosmos['Muv'], t_cosmos['VIS_FWHM_AS_GAUSS'],
                label='COSMOS', color='tab:blue', alpha=0.5, zorder=1)

        ax.scatter(t_cdfs['Muv'], t_cdfs['VIS_FWHM_AS_GAUSS'],
                label='CDFS', color='tab:red', edgecolor='black', alpha=0.8, zorder=2)

        density_contours(ax,
                        t_cosmos['Muv'],
                        t_cosmos['VIS_FWHM_AS_GAUSS'],
                        bins=30, smooth=1.5,
                        color='tab:blue', linewidths=2)

        # density_contours(ax,
        #                 t_cdfs['Muv'],
        #                 t_cdfs['VIS_FWHM_AS_GAUSS'],
        #                 bins=35, smooth=1.0,
        #                 color='tab:red', linewidths=2)

        ax.set_xlabel(r'M$_{\mathrm{UV}}$')
        ax.set_ylabel('FWHM (arcsec)')

        ax.axhline(0.208, color='black', linestyle='--',
                label='VIS PSF FWHM', linewidth=2, alpha=0.7, zorder=0)

        ax.set_xlim(-24.5, -19)
        ax.set_ylim(0., 2)

        ax.invert_xaxis()
        ax.legend(loc='upper left')

        plt.tight_layout()
        plt.savefig(plot_dir / 'candidates_z6_vis_fwhm_vs_muv.pdf')
        plt.show()

    if plot_3:
        #! Plot 3: FWHM vs photometric redshift for CDFS and COSMOS
        plt.figure(figsize=(10, 6))

        plt.scatter(t_cdfs['Zphot'], t_cdfs['VIS_FWHM_AS_GAUSS'],
                    label='CDFS', color='tab:red', edgecolor='black', alpha=0.8, zorder=2)
        
        plt.scatter(t_cosmos['Zphot'], t_cosmos['VIS_FWHM_AS_GAUSS'],
                    label='COSMOS', color='tab:blue', alpha=0.4, zorder=1)

        # Add line for VIS PSF FWHM
        plt.axhline(0.208, color='black', linestyle='--',
                    label='VIS PSF FWHM', linewidth=2, alpha=0.7, zorder=0)

        # Add lines at z=5.7 and z=6.3, cuts used in literature to minimise BD contamination
        plt.axvline(5.7, color='gray', linestyle='--', linewidth=2, alpha=0.7, zorder=0)
        plt.axvline(6.3, color='gray', linestyle='--', linewidth=2, alpha=0.7, zorder=0)

        plt.xlabel(r'$z_{\mathrm{phot}}$')
        plt.ylabel('FWHM (arcsec)')

        plt.xlim(5.4, 6.6)
        plt.ylim(0., 2)

        plt.legend(loc='upper left')

        plt.tight_layout()

        plt.savefig(plot_dir / 'candidates_z6_vis_fwhm_vs_zphot.pdf')

        plt.show()
        

        