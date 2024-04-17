"""
astrometric_accuracy.py

Plot the delta between RA and DEC in various images/catalogues.

Takes inputs from crossmatch_gaia_stars.py

Created: Friday 12th April 2024.

"""

from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
plot_dir = Path.cwd().parent.parent / 'plots' / 'astrometry'

# Run Gaia crossmatched stars or VISTA crossmatched stars? If False, run VISTA
gaia = False

psf_fwhm = [0.48, 0.51, 0.55] # Taken from pipeline PSFs, psf_fwhm_plots.py 0.21, 

filter_names = ['Y', 'J', 'H']

for i, filter_name in enumerate(filter_names):

    print(f'Processing filter {filter_name}')

    # Read the commented header stars.ascii files
    if gaia:
        stars_file = stars_dir / f'{filter_name}_euclid_gaia_coords.ascii' #! Stars used for PSFEx 
        #stars_file = stars_dir / f'{filter_name}_brightEuclid_gaia_coords.ascii' #! Brighter stars (inc. saturated)
    else:
        stars_file = stars_dir / f'{filter_name}_vista_euclid_coords.ascii'
    stars = ascii.read(stars_file)

    # Calculate the difference between the RA and DEC
    if gaia:
        delta_ra = stars['RA_euclid'] - stars['RA_Gaia']
        delta_dec = stars['DEC_euclid'] - stars['DEC_Gaia']
    else:
        delta_ra = stars['RA_euclid'] - stars['RA_vista']
        delta_dec = stars['DEC_euclid'] - stars['DEC_vista']

    # Convert delta coords to arcsec
    delta_ra *= 3600
    delta_dec *= 3600
    
    # Plot the difference
    plt.figure(figsize=(8, 8))
    plt.scatter(delta_ra, delta_dec, s=1, marker='x', color='black', alpha=0.6)
    plt.xlabel(r'$\Delta \mathrm{RA \ (as)}$')
    plt.ylabel(r'$\Delta \mathrm{Dec \ (as)}$')
    if gaia:
        plt.title(f'{filter_name} - Euclid vs Gaia')
    else:
        plt.title(f'{filter_name} - Euclid vs VISTA')

    # Plot a circle corresponding to the Euclid pixel scale
    circle1 = plt.Circle((0, 0), 0.1, color='r', fill=False, linewidth=2.5)
    plt.gca().add_artist(circle1)

    # And plot a Circle corresponding to the VIS FWHM peak
    circle2 = plt.Circle((0, 0), psf_fwhm[filter_names.index(filter_name)], color='black', fill=False, linewidth=2.5)
    plt.gca().add_artist(circle2)

    # Find the subset of stars that lie outside the Euclid pixel scale and the FWHM
    outside = np.where((np.abs(delta_ra) > 0.1) | (np.abs(delta_dec) > 0.1) | (np.sqrt(delta_ra**2 + delta_dec**2) > psf_fwhm[filter_names.index(filter_name)]))[0]
    
    # Save these to a file
    if gaia:
        outside_file = stars_dir / f'{filter_name}_outside_pixscale_euclid_gaia_coords.ascii'
    else:
        outside_file = stars_dir / f'{filter_name}_outside_pixscale_vista_euclid_coords.ascii'

    # Dummy plots for the legend
    plt.plot([-99, -98], [-99, -98], color='r', lw=2.5, label='Euclid pixel scale')
    plt.plot([-99, -98], [-99, -98], color='black', lw=2.5, label=f'{filter_name} PSF FWHM')

    plt.axis('equal')

    plt.xlim(-0.6, 0.6)
    plt.ylim(-0.6, 0.6)

    plt.tight_layout()
    plt.legend()

    #if gaia:
    #    plt.savefig(plot_dir / f'{filter_name}_bright_delta_ra_dec.png')
    #else:
    #    plt.savefig(plot_dir / f'VISTA_{filter_name}_delta_ra_dec.png')
    #plt.close()

    #plt.show()
    #plt.close()