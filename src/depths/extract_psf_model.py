#!/usr/bin/env python3

'''Get out PSF model from the .psf file made by PSFEx'''

from astropy.io import fits
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import glob

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

filter_names = ['VIS', 'Y', 'J', 'H']

colors = ['blue', 'green', 'red', 'purple']

plt.figure(figsize=(10, 6))

for filter_name in filter_names:

    psfDir = Path.home() / 'euclid' / filter_name / 'COSMOS'

    outputDir = Path.cwd().parent.parent / 'data' / 'psf' / 'COSMOS' / 'results'

    psf_files = glob.glob(str(psfDir / '*PSF*.fits'))

    all_fwhms = []
    all_ra = []
    all_dec = []

    for psf_file in psf_files:
        print(psf_file.split('/')[-1])
        hdu = fits.open(psf_file)

        header = hdu[2].header
        data = hdu[2].data

        fwhms = [t[-1] for t in data]
        Decs = [t[-2] for t in data]
        RAs = [t[-3] for t in data]
        all_fwhms.extend(fwhms)
        all_ra.extend(RAs)
        all_dec.extend(Decs)

    #plt.hist(fwhms, bins=np.arange(0.1, 0.6, 0.001), color=colors[filter_names.index(filter_name)], alpha=0.8, label=filter_name, density=True)
    plt.scatter(all_ra, all_dec, c=all_fwhms, cmap='viridis', s=10)


    plt.legend()
    plt.show()

#plt.ylabel('Normalised frequency')
#plt.xlabel('FWHM (arcsec or pixels???)')

