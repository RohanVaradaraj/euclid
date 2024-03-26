"""
collate_tile_depths.py

The euclid image is too big to run depths in one go, so we run on the individual mosaic tiles and then bring them together.
Here we make the maps and plots.

Created: Tueaday 26th March 2024.
"""


from astropy.table import Table, vstack
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.wcs import WCS
import glob
from shapely.geometry import Polygon, Point
from scipy.interpolate import griddata
from scipy import stats

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

plot_dir = Path.cwd().parent.parent / 'plots'

def collate_tile_depths(filter_name: str, field_name: str, ap_diametersAS: float, mask=False) -> Table:

    '''
    Opens the local depth tables of the individual Euclid tiles and collates them into a single table.

    Parameters
    ----------

    filter_name : str
        The filter we want to collate.
    field_name : str
        The field name.
    ap_diametersAS : float
        The aperture diameter used to measure depths in arcseconds.
    mask : bool
        Whether to mask the depths to the Euclid region.
    
    Returns
    -------

    local_depths: astropy.table.Table
        A large table of local depths at positions x, y (and RA, DEC) across the full masked Euclid mosaic.
    '''

    # Set up directories
    table_dir = Path.cwd().parent.parent / 'data' / 'depths' / field_name / 'phot'
    image_dir = Path.home() / 'euclid' / filter_name / field_name
    mask_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'masks'
    grid_dir = Path.cwd().parent.parent / 'data' / 'grids'

    # Get all files with the filter name and aperture diameter.
    files = glob.glob(str(table_dir / f'*{filter_name}*{ap_diametersAS}as_gridDepths_300_200.fits*'))
    print(len(files))

    # Empty list to store each table
    tables = []
    tiles = []

    for filename in files:

        # Open the corresponding image header
        tile_index = filename.split('/')[-1].split('_')[1]
        if filter_name != 'VIS':
            image_name = glob.glob(str(image_dir/f'EUC_MER_BGSUB-MOSAIC-NIR-{filter_name}_TILE_*_{tile_index}.fits'))
        else:
            image_name = glob.glob(str(image_dir/f'EUC_MER_BGSUB-MOSAIC-VIS_TILE_*_{tile_index}.fits'))

        if len(image_name) == 0:
            continue

        image_name = image_name[0]

        tiles.append(tile_index)

        # Get wcs
        hdu = fits.open(image_name)
        header = hdu[0].header
        wcs = WCS(header)

        # Open the table
        t = Table.read(filename)

        # Add the RA and DEC columns
        t['RA'], t['DEC'] = wcs.all_pix2world(t['x'], t['y'], 0)

        # Append to the list
        tables.append(t)

    # Stack the tables
    local_depths = vstack(tables)

    # Mask the depths
    if mask:

        with open(mask_dir / f'Euclid_square_COSMOS.reg', 'r') as file:
            
            for line in file:
                if line.startswith("polygon"):
                    # Split the line by commas and remove the "polygon(" part
                    coords = line.split('(')[1].rstrip(')\n').split(',')
                    coords = [float(coord) for coord in coords]

        coords = [(coords[i], coords[i+1]) for i in range(0, len(coords), 2)]

        # Set up shapely polgyon and check if coordinates lie within it.
        polygon = Polygon(coords)
        points_inside_polygon = [Point(ra, dec).within(polygon) for ra, dec in zip(local_depths['RA'], local_depths['DEC'])]
        local_depths = local_depths[points_inside_polygon]


    local_depths.write(grid_dir / f'{field_name}_{filter_name}_{ap_diametersAS}as_gridDepths.fits', overwrite=True)
    return local_depths



def plotDepthGrid(t: Table, filter_name: str, ap_diametersAS: float) -> None:
    '''
    Plot the depth grid in RA and DEC.

    Parameters
    ----------
    t: astropy.table.Table
        The table of depths.
    filter_name: str
        The filter name.
    ap_diametersAS: float
        The aperture diameter used to calculate depth, in arcseconds.

    Returns
    -------
    None
    '''

    # Determine the bounding box of the data
    min_RA, max_RA = np.min(t['RA']), np.max(t['DEC'])
    min_DEC, max_DEC = np.min(t['RA']), np.max(t['DEC'])

    # Create a meshgrid covering the square region
    grid_RA, grid_DEC = np.meshgrid(np.linspace(min_RA, max_RA, 100), np.linspace(min_DEC, max_DEC, 100))


    # Interpolate depth values onto the meshgrid
    interpolated_depths = griddata((t['RA'], t['DEC']), t['depths'], (grid_RA, grid_DEC), method='linear')

    qlow, qhigh = np.percentile(t['depths'], [1, 99])

    plt.figure(figsize=(10,10))
    plt.scatter(t['RA'], t['DEC'], c=t['depths'], s=30, vmin=qlow, vmax=qhigh, alpha=0.8)
    plt.xlabel('RA (degrees)')
    plt.ylabel('Dec (degrees)')
    plt.gca().invert_xaxis()
    plt.colorbar(label='Depth')
    plt.gca().set_aspect('equal', adjustable='box')
    print('Plotting...')
    plot_name = plot_dir / f'{filter_name}_{ap_diametersAS}as_depth_grid.png'
    plt.savefig(plot_name)



def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth



def plotDepthHist(t: Table, filter_name: str, ap_diametersAS: float) -> None:

    '''
    Plot a histogram of the depths.

    Parameters
    ----------
    t: astropy.table.Table
        The table of depths.
    filter_name: str
        The filter name.
    ap_diametersAS: float
        The aperture diameter used to calculate depth, in arcseconds.
    
    Returns
    -------
    None
    '''

    mags = t['depths']

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

    print(bins[modeMag[0]])

    # Mode
    modeMag = round(bins[modeMag[0]], 2)

    # Clear previous plot so that saving figures works!
    plt.cla()
    plt.clf()

    # Histogram
    n, binsMode = np.histogram(mags, bins=bins)
    plt.step(binsMode[0:-1], n, color='orange', alpha=0.6, lw=2)

    # Smooth the mags for calculation.
    smoothMags = smooth(n, 20)

    # Split the peaks for UltraVISTA
    split = int(0.99*len(smoothMags))

    # Histogram and plots of smoothed mags.
    vals, valBins = np.histogram(smoothMags[0:split], bins=bins)
    plt.plot(binsMode[0:-1], smoothMags, color='black', linestyle='None', marker='.', linewidth=1, markersize=4)

    # Mode
    modeMag = round(binsMode[np.argmax(smoothMags)], 2)

    qlow, qhigh = np.percentile(mags, [2, 99])

    plt.axvline(medianMag, color='r', linestyle='dashed', label='median={0}'.format(round(medianMag, 2)), lw=3)
    plt.axvline(meanMag, color='green', linestyle='dashed', label='mean={0}'.format(round(meanMag, 2)), lw=3)
    plt.axvline(modeMag, color='orange', linestyle='dashed', label='mode={0}'.format(modeMag), lw=3)

    plt.title('{0} {1}'.format(filter_name, str(ap_diametersAS)+ ' as'))
    plt.xlim(left=qlow)
    plt.xlim(right=qhigh)
    plt.xlabel('mag')
    plt.ylabel('Count')
    plt.legend()

    plot_name = plot_dir / f'{filter_name}_{ap_diametersAS}as_depth_hist.png'
    plt.savefig(plot_name)



if __name__ == "__main__":

    for filter_name in ['H']:

        t = collate_tile_depths(filter_name, 'COSMOS', 1.0, mask=True)

        t = t[t['mask']==0.0]
        t = t[np.isfinite(t['depths'])]

        plotDepthGrid(t, filter_name, 1.0)
        plotDepthHist(t, filter_name, 1.0)



    