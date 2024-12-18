"""
plot some things related to our LAE

Created: Tuesday 17th December 2024.
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import simps
import matplotlib.ticker as mticker

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

def load_obj(ID):

    """ Load objects from the xmm or cdfs catalogue in the format bagpipes needs."""

    # Open catalogue
    cat = Table.read(catDir / catName)

    # Extract object with ID
    #cat['ID'].format

    # BAGPIPES ID is just the row number
    ID = 178396
    cat = cat[:][cat['ID'] == ID]

    # # But also print our ID for this object.
    # print('########## OBJECT ID: ' + str(cat['ID']) + '##############')

    # Empty arrays for getting names of columns
    fluxes = []
    fluxerr = []

    # Loop through filters
    for i, filterName in enumerate(filters):

        # Append flux name
        fluxes = fluxes + ['flux_{0}'.format(filterName)]

        # Append flux error name
        fluxerr = fluxerr + ['err_{0}'.format(filterName)]

    # Now get the flux
    flux = cat[fluxes][0][:]
    flux = np.array(flux)
    error = cat[fluxerr][0][:]
    error = np.array(error)

    # print the ch1 and ch2 magnitudes
    mag_ch1 = -2.5*np.log10(cat['flux_ch1cds'][0]) - 48.6
    mag_ch2 = -2.5*np.log10(cat['flux_ch2cds'][0]) - 48.6
    print('ch1 mag: ', mag_ch1)
    print('ch2 mag: ', mag_ch2)

    # Convert to microJanskys
    flux = flux * 1e29 #* 3.631 * 10**28.44
    error = error * 1e29 #* 3.631 * 10**28.44

    # Turn these into a 2D array
    photometry = np.c_[flux, error]

    return flux, error

# Load the bagpipes fits
table_dir = Path.cwd().parent / 'sed_fitting' / 'pipes'
sed = Table.read(table_dir / '178396_microJy.fits')
model_phot = Table.read(table_dir / '178396_microJy_synthpoints.fits')

# Load actual phot
catDir = Path.cwd().parents[1] / 'data' / 'catalogues'
catName = 'LAE.fits'
filters = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'VIS', 'Ye', 'Je', 'He', 'f115w', 'f150w', 'f277w', 'f444w'] #, 'ch1cds', 'ch2cds'] 
#filters = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'VIS', 'Ye', 'Je', 'He', 'ch1cds', 'ch2cds']
phot, err = load_obj(178396)

# Replace fluxes for ch1 and ch2 with tractor values from COSMOS2020
phot[-1] = 0.39101
err[-1] = 0.00924
phot[-2] = 0.2084
err[-2] = 0.00958

phot *= 1e-29
err *= 1e-29

fig, ax1 = plt.subplots(figsize=(10, 6))


# Main plot: Flux vs Wavelength
ax1.plot(sed['wave'], sed['flux'] * 1e-29, c='deepskyblue', alpha=1, lw=2.5, zorder=-1)
ax1.set_yscale('log')
ax1.set_ylim(1e-31, 1e-28)
ax1.set_xlim(3000, 50000)

# Synthetic photometry (open gray circles)
ax1.scatter(model_phot['wave'], model_phot['flux'] * 1e-29, edgecolor='gray', s=120, facecolors='none', label="Model photometry", lw=3)

# Actual photometry with error bars
for i in range(len(phot)):
    if phot[i] / err[i] < 2:
        ax1.errorbar(model_phot['wave'][i], phot[i] + 2 * err[i], yerr=0, fmt='v', c='black', markersize=7)
    else:
        ax1.errorbar(model_phot['wave'][i], phot[i], yerr=err[i], fmt='o', c='black', markersize=7, elinewidth=2)

# Labels for the primary axis
ax1.set_xlabel(r'$\lambda (\mu m)$')
ax1.set_ylabel(r'$f_{\nu}$ (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')

# Secondary y-axis for magnitude
def flux_to_mag(flux, zp=48.6):
    """Convert flux to AB magnitude."""
    return -2.5 * np.log10(flux) - zp

def mag_to_flux(mag, zp=48.6):
    """Convert AB magnitude to flux."""
    return 10 ** (-0.4 * (mag + zp))

# Add magnitude axis
secax = ax1.secondary_yaxis('right', functions=(flux_to_mag, mag_to_flux))

mags = [30, 29, 28, 27, 26, 25, 24, 23, 22]
secax.set_yticks(mags)

secax.yaxis.set_major_formatter(mticker.ScalarFormatter())

ax1.tick_params(which='major', length=5, width=3)
ax1.tick_params(axis='both', which='minor', length=3, width=2)

secax.tick_params(which='major', length=5, width=3)
secax.tick_params(axis='both', which='minor', length=5, width=2)

secax.set_ylabel(r'$\mathrm{m_{AB}}$')

# Convert the x axis into microns by dividing by 10000
ax1.set_xticks([5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000])
ax1.set_xticklabels([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])

# Add legend and show the plot
ax1.legend(loc='upper left', fontsize=15)
plt.tight_layout()
plot_dir = Path.cwd().parents[1] / 'plots' / 'seds'
plt.savefig(plot_dir / '178396_BAGPIPES_CWEB_fit.pdf')
plt.show()
exit()

# Load wavelength and flux data (replace with your SED data)
wave = sed['wave']  # Wavelengths in Angstroms
flux = sed['flux']*1e-29  # Flux in erg/s/cm^2/Hz

# Define the wavelength range for the [OIII] and H-beta lines
line_region = (wave > 39400) & (wave < 41300)

# Define continuum regions on either side of the lines
cont_region_left = (wave > 37500) & (wave < 39400)
cont_region_right = (wave > 41300) & (wave < 47500)

# Estimate the continuum level using a linear fit
cont_wave = np.concatenate([wave[cont_region_left], wave[cont_region_right]])
cont_flux = np.concatenate([flux[cont_region_left], flux[cont_region_right]])

# Fit a linear continuum
p = np.polyfit(cont_wave, cont_flux, 1)
continuum_fit = np.polyval(p, wave[line_region])

# Subtract continuum from flux in the line region
line_flux = flux[line_region] - continuum_fit

# Integrate the line flux using the trapezoidal rule
EW = simps(line_flux / continuum_fit, wave[line_region])

# Print the equivalent width
print(f"Observed Equivalent Width of [OIII] + H-beta: {EW:.2f} Angstroms")

# Correct for redshift
lae_z = 7.19309
EW_0 = EW / (1 + lae_z)
print(f"Rest-frame Equivalent Width of [OIII] + H-beta: {EW_0:.2f} Angstroms")

# Plot the result
plt.plot(wave, flux, label='SED')
plt.plot(wave[line_region], continuum_fit, label='Fitted Continuum', linestyle='--')
plt.fill_between(wave[line_region], continuum_fit, flux[line_region], color='gray', alpha=0.5, label='Line Flux')
plt.ylim(1e-31, 1e-28)
plt.xlim(3000, 50000)
plt.xlabel("Wavelength (Angstroms)")
plt.ylabel("Flux")
plt.yscale('log')
plt.legend()
plt.show()
