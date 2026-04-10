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
from astropy.nddata.utils import Cutout2D, NoOverlapError
from astropy.wcs import WCS
import glob
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from astropy.table import Table
import astropy.units as u

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

# plt.rcParams['axes.linewidth'] = 2.5
# plt.rcParams.update({'font.size': 15})
# plt.rcParams['figure.dpi'] = 100

ground_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data'
euclid_dir = Path.home() / 'euclid'
cweb_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'CWEB'
primer_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS'
plot_dir = Path.cwd().parent.parent / 'plots' / 'cutouts' / 'z4'
# Check if exists, if not then make


refcat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'

def save_cutout_fits(cutout, filename):
    hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
    hdu.writeto(filename, overwrite=True)

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
    
    # Sigma clip the data
    data = data[(data < mean + 3 * std_dev)]

    # Recalculate mean and std_dev
    mean = np.mean(data)
    std_dev = np.std(data)

    lower = mean - 2 * std_dev
    upper = mean + 5 * std_dev

    return lower, upper

def truncate(value, decimals):
    factor = 10 ** decimals
    return int(value * factor) / factor


def make_euclid_name(ra_deg, dec_deg):
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    ra_h = c.ra.hms
    dec_d = c.dec.dms

    # Truncate RA to two d.p.
    ra_sec_trunc = truncate(ra_h.s, 2)
    ra_str = f"{int(ra_h.h):02d}{int(ra_h.m):02d}{int(ra_sec_trunc):02d}.{int((ra_sec_trunc % 1) * 100):02d}"

    # Truncate DEC to one d.p.
    sign = '+' if dec_d.d >= 0 else '-'
    dec_d_abs = abs(dec_d.d)
    dec_m_abs = abs(dec_d.m)
    dec_s_trunc = truncate(abs(dec_d.s), 1)
    dec_str = f"{sign}{int(dec_d_abs):02d}{int(dec_m_abs):02d}{int(dec_s_trunc):02d}.{int((dec_s_trunc % 1) * 10):1d}"

    return f"J{ra_str}${sign}${dec_str[1:]}"  

def isCoordInSurveyFootprints(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """
    Check if a coordinate lies in the footprints of:
        1) Euclid On-the-fly COSMOS release
        2) COSMOS-Web
        3 )PRIMER

    Parameters
    ----------
    ra : np.array
        Right ascension of the coordinate in degrees.
    dec : np.array
        Declination of the coordinate in degrees.


    Returns
    -------
    np.ndarray
        Array with rows corresponding to the input coordinates and columns corresponding to the footprints.
        The value is '0' if it is not in the footprint.
        If it is in the footprint then: 
            for Euclid it takes a numeric string value of the tile label I have assigned (1-22).
            for CWEB it takes a numeric string value of the tile label Nathan has assigned (4A, 5A, 5B, 6A, etc.)
            for PRIMER it takes the value '1'.

    """

    # Check if ra and dec are arrays
    ra = np.array(ra)
    dec = np.array(dec)

    # Reshape ra and dec to make sure they are at least one-dimensional arrays
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    assert ra.shape == dec.shape, "ra and dec arrays must have the same shape"

    points = [Point(x, y) for x, y in zip(ra, dec)]

    # Load the footprints
    euclid_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_footprint.npy')
    primer_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'primer_footprint.npy')
    cweb_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_footprint.npy')
    hubble_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'hubble_footprint.npy')

    # Load the labels for CWEB and Euclid
    euclid_labels = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_labels.npy')
    cweb_labels = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_labels.npy')

    # Create shapely polygons
    euclid_polygons = [Polygon(footprint) for footprint in euclid_footprint]
    cweb_polygons = [Polygon(footprint) for footprint in cweb_footprint]
    primer_polygons = Polygon(primer_footprint)
    hubble_polygons = Polygon(hubble_footprint)

    points_in_euclid = []
    points_in_cweb = []
    points_in_primer = []
    points_in_hubble = []

    for point in points:

        # 0 if not in footprint, else the tile label (1 for HST and PRIMER since it is one tile each)
        in_euclid = '0'
        in_cweb = '0'
        in_primer = '0'
        in_hubble = '0'

        #! Check if the point lies in any of the Euclid polygons
        for i, polygon in enumerate(euclid_polygons):
            if polygon.contains(point):
                in_euclid = euclid_labels[i]
                break

        #! Check if the point lies in any of the COSMOS-Web polygons
        for i, polygon in enumerate(cweb_polygons):
            if polygon.contains(point):
                in_cweb = cweb_labels[i]
                break

        #! Check if point lies in PRIMER footprint
        if primer_polygons.contains(point):
            in_primer = '1'

        #! Check if point lies in Hubble footprint
        if hubble_polygons.contains(point):
            in_hubble = '1'

        # Append results of search to lists
        points_in_euclid.append(in_euclid)
        points_in_cweb.append(in_cweb)
        points_in_primer.append(in_primer)
        points_in_hubble.append(in_hubble)
    
    # Convert lists to numpy arrays
    euclid_mask = np.array(points_in_euclid)
    cweb_mask = np.array(points_in_cweb)
    primer_mask = np.array(points_in_primer)
    hubble_mask = np.array(points_in_hubble)

    return np.column_stack((euclid_mask, cweb_mask, primer_mask, hubble_mask))


def isCoordInCWEB(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """
    Check if a coordinate lies in the footprints of CWEB.

    If it doesn't, return zero for the coord, else return the mosaic string.

    Parameters
    ----------
    ra : np.array
        Right ascension of the coordinate in degrees.
    dec : np.array
        Declination of the coordinate in degrees.
    
    Returns
    -------
    np.ndarray
        Array with rows corresponding to the input coordinates and columns corresponding to the footprints.
        The value is '0' if it is not in the footprint.
        If it is in the footprint then it takes a numeric string value of the tile label Nathan has assigned (4A, 5A, 5B, 6A, etc.)

    """

    # Check if ra and dec are arrays
    ra = np.array(ra)
    dec = np.array(dec)

    # Reshape ra and dec to make sure they are at least one-dimensional arrays
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    assert ra.shape == dec.shape, "ra and dec arrays must have the same shape"

    points = [Point(x, y) for x, y in zip(ra, dec)]

    mosaics = ['0A', '0B', '1A', '1B', '2A', '2B', '3A', '3B', '4A', '4B', '5A', '5B', '6A', '6B', '7A', '7B']

    polygons = []

    for mosaic in mosaics:
        footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / f'CWEB_footprint_{mosaic}.npy')

        # Create shapely polygons
        cweb_polygon = Polygon(footprint)
        polygons.append(cweb_polygon)

    points_in_cweb = []
    for point in points:

        # 0 if not in footprint, else the tile label
        in_cweb = '0'

        #! Check if the point lies in any of the COSMOS-Web polygons
        #! If so, find all of them
        for i, polygon in enumerate(polygons):
            if polygon.contains(point):
                in_cweb = mosaics[i]
                break

        # Append results of search to lists
        points_in_cweb.append(in_cweb)

    # Convert lists to numpy arrays
    points_in_cweb = np.array(points_in_cweb)

    return points_in_cweb




def Cutout(ra: float, dec: float, ID=None, contained_in: Optional[np.array] = None, size: float = 10.0,
           save_cutout: bool = True, save_dir: Path = Path.cwd().parent.parent / 'data' / 'cutouts',
           plot_title: Optional[str] = None,
           add_centre_lines: Optional[bool] = False) -> None:
    """
    Create and plot *only CWEB* cutouts (VIS, Y, J, H) around (ra, dec).

    Notes
    -----
    - This function *only* looks for CWEB tiles. It reads the second element of the
      contained_in array (CWEB tile label). If that's '0' the function returns None.
    - It attempts robust globbing for files; if a band tile is missing it makes a
      zero-filled cutout with the correct WCS (so the layout stays consistent).
    - The function plots the 2x2 stamp grid (VIS, Y, J, H). If no CWEB tile -> None.
    """

    import glob
    import matplotlib.pyplot as plt

    # ensure contained_in known
    if contained_in is None:
        contained_in = isCoordInSurveyFootprints(ra, dec)

    # second element -> CWEB tile label
    # contained_in shape assumed similar to original: contained_in[0][1] is CWEB tile
    try:
        cweb_tile = contained_in[0][1]
    except Exception:
        # incompatible input format
        print("contained_in has unexpected format; expected contained_in[0][1] to be CWEB tile label.")
        return None

    if cweb_tile == '0' or cweb_tile is None:
        # not in CWEB footprint
        print("Coordinate not in CWEB footprint; no cutout made.")
        return None

    # prepare skycoord
    c = SkyCoord(ra, dec, unit='deg')

    # directory where CWEB tiles live; try sensible location(s)
    # adapt these paths to your environment if different
    cweb_dir_candidates = []
    try:
        cweb_dir_candidates.append(euclid_dir / 'CWEB' / 'COSMOS')
    except NameError:
        pass
    try:
        cweb_dir_candidates.append(Path.cwd().parent.parent / 'data' / 'CWEB' / 'COSMOS')
    except Exception:
        pass
    # final fallback: euclid_dir / 'COSMOS' / 'CWEB'
    try:
        cweb_dir_candidates.append(euclid_dir / 'COSMOS' / 'CWEB')
    except Exception:
        pass

    # bands and preferred filename patterns (ordered)
    bands = {
        'VIS': ['*VIS*', 'CWEB_VIS', 'VIS'],
        'Y':   ['*Y*', 'CWEB_Y', 'Y'],
        'J':   ['*J*', 'CWEB_J', 'J'],
        'H':   ['*H*', 'CWEB_H', 'H'],
    }

    cutouts = {}
    wcs_ref = None  # keep a reference WCS if available

    for band, patterns in bands.items():
        tile_path = None
        for base in cweb_dir_candidates:
            if base is None:
                continue
            for patt in patterns:
                glob_pattern = str(base / f'*{patt}*_{cweb_tile}.fits')
                found = glob.glob(glob_pattern)
                if not found:
                    # try without the tile suffix (some files may not include tile id)
                    found = glob.glob(str(base / f'*{patt}*.fits'))
                if found:
                    tile_path = found[0]
                    break
            if tile_path:
                break

        if tile_path is None:
            # fallback: try any fits under candidate dirs containing band name
            for base in cweb_dir_candidates:
                if base is None:
                    continue
                found = glob.glob(str(base / f'*{band}*.fits'))
                if found:
                    tile_path = found[0]
                    break

        if tile_path is None:
            # no file available -> create zero-filled fake cutout using last-known WCS if present
            if wcs_ref is None:
                # try to find any fits in the candidate dirs to get a WCS
                for base in cweb_dir_candidates:
                    if base is None:
                        continue
                    any_fits = glob.glob(str(base / '*.fits'))
                    if any_fits:
                        with fits.open(any_fits[0]) as hdu_any:
                            wcs_ref = WCS(hdu_any[0].header)
                        break

            if wcs_ref is None:
                # cannot determine WCS; skip band
                cutouts[band] = None
                continue

            if np.isscalar(size):
                shape = (int(size / (np.abs(wcs_ref.wcs.cdelt[0]) * 3600)),
                         int(size / (np.abs(wcs_ref.wcs.cdelt[0]) * 3600)))
            else:
                # size provided as (x_arcsec, y_arcsec)
                shape = (int(size[1] / (np.abs(wcs_ref.wcs.cdelt[0]) * 3600)),
                         int(size[0] / (np.abs(wcs_ref.wcs.cdelt[0]) * 3600)))

            empty_data = np.zeros(shape, dtype=np.float32)
            cutouts[band] = Cutout2D(empty_data,
                                     position=(shape[1] // 2, shape[0] // 2),
                                     size=shape,
                                     wcs=wcs_ref)
            continue

        # file exists: open and cutout
        with fits.open(tile_path) as hdu:
            data = hdu[0].data
            hdr = hdu[0].header
            try:
                pix_scale = np.abs(hdr['CD1_1']) * 3600
            except Exception:
                # try cdelt fallback
                try:
                    pix_scale = np.abs(hdr['CDELT1']) * 3600
                except Exception:
                    pix_scale = 1.0  # best-effort fallback
            wcs_band = WCS(hdr)
            wcs_ref = wcs_ref or wcs_band

            try:
                cut = Cutout2D(data, c, size=size / pix_scale, wcs=wcs_band)
            except NoOverlapError:
                # produce zero-filled cutout with same WCS
                if np.isscalar(size):
                    shape = (int(size / pix_scale), int(size / pix_scale))
                else:
                    shape = (int(size[1] / pix_scale), int(size[0] / pix_scale))
                empty_data = np.zeros(shape, dtype=data.dtype)
                cut = Cutout2D(empty_data,
                               position=(shape[1] // 2, shape[0] // 2),
                               size=shape,
                               wcs=wcs_band)
            cutouts[band] = cut

    # If no valid cutouts at all -> return None
    if all((cutouts.get(b) is None or cutouts[b].data.size == 0) for b in bands.keys()):
        return None

    # Plot 2x2 grid: VIS | Y
    #               J   | H
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    axes = axes.flatten()
    vmax_vmin = {}
    # compute percentiles per-band for sensible stretching
    for i, band in enumerate(['VIS', 'Y', 'J', 'H']):
        ax = axes[i]
        cut = cutouts.get(band)
        if cut is None:
            ax.axis('off')
            ax.set_title(f"{band} (missing)")
            continue
        arr = cut.data
        # handle masked or NaN data
        arr = np.array(arr, dtype=np.float32)
        if arr.size == 0 or np.all(np.isnan(arr)):
            ax.axis('off')
            ax.set_title(f"{band} (empty)")
            continue
        # robust vmin/vmax
        try:
            vmin, vmax = np.nanpercentile(arr, [2, 98])
            if vmin == vmax:
                vmin = np.nanmin(arr)
                vmax = np.nanmax(arr) if np.nanmax(arr) != vmin else vmin + 1.0
        except Exception:
            vmin, vmax = np.nanmin(arr), np.nanmax(arr)
            if vmin == vmax:
                vmax = vmin + 1.0
        vmax_vmin[band] = (vmin, vmax)
        ax.imshow(arr, origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(band)
        ax.set_xticks([])
        ax.set_yticks([])

        if add_centre_lines:
            ny, nx = arr.shape
            ax.axvline(nx / 2, linestyle='--', linewidth=0.6)
            ax.axhline(ny / 2, linestyle='--', linewidth=0.6)

    if plot_title:
        fig.suptitle(plot_title)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # create save directory if requested
    if save_cutout:
        save_dir.mkdir(parents=True, exist_ok=True)
        # save figure
        fname_fig = save_dir / f'ID_{ID}_CWEB_cutouts_{ra:.5f}_{dec:.5f}.png'
        fig.savefig(fname_fig, dpi=150, bbox_inches='tight')
        print(f'Saved cutout figure to {save_dir / fname_fig}')
        # also save individual fits cutouts (if present)
        for band, cut in cutouts.items():
            if cut is None:
                continue
            out_f = save_dir / f'ID_{ID}_CWEB_{band}_{ra:.5f}_{dec:.5f}.fits'
            save_cutout_fits(cut, out_f)

    plt.show()
    #plt.close(fig)

    return None




if __name__ == '__main__':
    

    #! HzRSs
    # 149.9179524404775,2.392064687286876
    # 150.00206688878058,2.3821626758186336
    # 150.1424249418437,1.9593914588078478

    ra = [149.9179524404775, 150.00206688878058, 150.1424249418437]
    dec = [2.392064687286876, 2.3821626758186336, 1.9593914588078478]
    IDs = ['MGTJ09594+02233', 'MGTJ10000+02225', 'MGTJ10003+01573']



    ############! GET CUTOUTS ############
    for i in range(len(ra)):

        print(f'Object number {i+1} of {len(ra)}')
        print(IDs[i])
        # continue

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
        #Cutout(ra[i], dec[i], size=6., add_centre_lines=True, save_cutout=False, plot_title=ID[i])
        Cutout(ra[i], dec[i], size=10, ID=IDs[i], save_cutout=True, plot_title=IDs[i]) #, add_centre_lines=True, save_cutout=False)
        #Cutout(ra[i], dec[i], size=600, plot_title=names[i], save_cutout=True)
        #Cutout(ra[i], dec[i], size=6., plot_title=ID[i])
        #Cutout(ra[i], dec[i], size=10., plot_title='Big Three Dragons')   #
        #Cutout(ra[i], dec[i], size=4., plot_title=ID[i] + ', z=' + str(z[i]) + ', Muv=' + str(Muv[i]))


    
