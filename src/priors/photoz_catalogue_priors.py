"""
Use Natalia's photo-z catalogues to determine a prior on the P(z) of my high-redshift sources.

Created: Friday 23rd May 2025.
"""

from astropy.table import Table
import numpy as np
import os
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

field_name = 'COSMOS'

#! Define prior functional form
# def prior_func(z, C, gamma, z_0):
#     """
#     Define the functional form of the prior.
#     """
#     return C * np.power(z, gamma) * np.exp(-np.power(z / z_0, gamma))
def prior_func(z, C, gamma, z_0, beta):
    return C * z**gamma / (1 + (z / z_0)**beta)


if __name__ == '__main__':

    #! Natalia's catalogues
    cat_dir = Path.cwd().parents[3] / 'natalia' / 'photo-z_catalogues'

    field_name_dict = {
        'XMM': 'XMM-LSS',
        'CDFS': 'E-CDFS',
        'COSMOS': 'COSMOS_DR6'
    }

    field_cat_dict = {
        'XMM': 'XMM_DR3_MASKED_Ks-stellar-cut_peters-LP-photo-z_LPz-AGN_chi2_LPclass_ModBest_LPmasses_LPsfr_LPssfr_flag_20-05-2025.fits',
        'CDFS': 'ECDFS_MASKED_Ks-stellar-cut_LP-photo-z_LPz-AGN_chi2_LPclass_ModBest_LPmasses_LPsfr_LPssfr_flag_Spec-z_20-05-2025.fits',
        'COSMOS': 'COSMOS_DR6_HSC-DR3_MASKED_Ks-stellar-cut_peters-LP-photo-z_LPz-AGN_chi2_LPclass_ModBest_LPmasses_LPsfr_LPssfr_flag_25-05-2025.fits'
    }

    # Load the catalogue
    cat_file = cat_dir / field_name_dict[field_name] / field_cat_dict[field_name]

    t = Table.read(cat_file)
    print(t.colnames)

    t = t[t['MAG_APER_1.8_HSC-Z'] < 99.]

    # Histogram of all the redshifts
    # plt.hist(t['Z_BEST_peak'], bins=np.arange(0, 10, 0.1), histtype='step', lw=2, color='k')
    # plt.yscale('log')
    # plt.show()
    # exit()

    # plt.hist(t['MAG_APER_1.8_HSC-Z'], bins=np.arange(15, 30, 0.1), histtype='step', lw=2, color='k', label='HSC-Z')
    # plt.xlabel(r'$m_{Z}$')
    # plt.ylabel('Count')
    # plt.title('XMM-LSS')
    # plt.show()

    # Bin the table in terms of MAG_APER_1.8_HSC-Z
    mag_bins = np.arange(18, 28, 1)

    cmap = cm.viridis
    norm = mcolors.Normalize(vmin=min(mag_bins), vmax=max(mag_bins))

    with open('prior_fits.txt', 'a') as f:
        f.write(f'# Field name \t mag \t C \t gamma \t z_0 \t beta \n') 

    for mag_bin in mag_bins:
        # Restrict table to these magnitudes
        t_restricted = t[(t['MAG_APER_1.8_HSC-Z'] > mag_bin) & (t['MAG_APER_1.8_HSC-Z'] < mag_bin + 1)]
        print(f'Number of sources in bin {mag_bin}: {len(t_restricted)}')

        # Plot the redshift distribution of these sources
        #plt.hist(t_restricted['Z_BEST_peak'], bins=np.arange(0, 4, 0.1), histtype='step', lw=2, color=cmap(norm(mag_bin)), alpha=0.5, density=False)

        # Histogram
        hist, bin_edges = np.histogram(t_restricted['Z_BEST_peak'], bins=np.arange(0.05, 10, 0.1), density=True)
        z_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

        mask = (hist > 0) & (z_centers > 0.05)
        z_fit = z_centers[mask]
        hist_fit = hist[mask]

        # Initial guess and bounds
        p0 = [1.0, 2.0, 0.5, 4.0]
        bounds = ([0, 0.5, 0.1, 2], [100, 8, 5, 15]) # C, gamma, z_0, beta

        popt, _ = curve_fit(prior_func, z_fit, hist_fit, p0=p0, bounds=bounds)


        try:
            popt, _ = curve_fit(prior_func, z_fit, hist_fit, p0=p0, bounds=bounds, maxfev=10000)
            fit_vals = prior_func(z_centers, *popt)

            # Normalise the fit_vals curve such that the integral of the curve is 1
            integral = np.trapz(fit_vals, z_centers)
            fit_vals /= integral

            # Normalise the histogram
            hist /= np.trapz(hist, z_centers)

            plt.plot(z_centers, hist, drawstyle='steps-mid', lw=2, alpha=0.3, color=cmap(norm(mag_bin)))
            plt.plot(z_centers, fit_vals, lw=3, alpha=0.8, color=cmap(norm(mag_bin)))
            print(f"Fit: bin {mag_bin}-{mag_bin+2} → C={popt[0]:.2f}, γ={popt[1]:.2f}, z₀={popt[2]:.2f}, β={popt[3]:.2f}")
        except RuntimeError as e:
            print(f"⚠️ Fit failed in bin {mag_bin}-{mag_bin+2}: {e}")

        # Save the fitted values to a file
        with open('prior_fits.txt', 'a') as f:
            f.write(f'{field_name} {mag_bin} {popt[0]} {popt[1]} {popt[2]} {popt[3]}\n')

    plt.xlabel(r'$z$')
    plt.ylabel(r'$p(z|m_{\mathrm{HSC-Z}})$')
    plt.tight_layout()
    plt.text(2, 2, r'$18 < m_{\mathrm{HSC-Z}} < 28$')
    #plt.yscale('log')

    plt.show()

