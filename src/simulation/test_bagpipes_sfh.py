#!/usr/bin/env python3

"""
Playing around with SFHs to produce emission lines

Created: Thursday 11th July 2024
"""

from brown_dwarf_colours import makeLBG
import matplotlib.pyplot as plt
import bagpipes as pipes
import numpy as np
from pathlib import Path

filter_list = np.loadtxt(str(Path.cwd() / 'filter_list.txt'), dtype="str")

sfh = {}
sfh["massformed"] = 9.  # Log_10 total stellar mass formed: M_Solar
sfh["metallicity"] = 0.1 # Metallicity: Z_sol = 0.02

delayed = {}                         # Delayed Tau model t*e^-(t/tau)
delayed["age"] = 0.5           # Time since SF began: Gyr
delayed["tau"] = 0.1           # Timescale of decrease: Gyr

burst = {}                           # Delta function burst
burst["age"] = 0.01            # Time since burst: Gyr

nebular = {}
nebular["logU"] = -1.          # Log_10 of the ionization parameter

dust = {}
dust["type"] = 'Calzetti'   # Attenuation law: "Calzetti", "Cardelli", "CF00" or "Salim"
dust["Av"] = 0.1   # Absolute attenuation in the V band: magnitudes

model_components = {}
model_components["redshift"] = 7 # Observed redshift
model_components["t_bc"] = 0.01           # Max age of birth clouds: Gyr

model_components["sfh"] = sfh   # Dict containing SFH info
model_components["delayed"] = delayed     # Dict containing delayed SFH info
model_components["burst"] = burst         # Dict containing burst SFH info
model_components["nebular"] = nebular     # Dict containing nebular info

model_components["dust"] = dust # Dict containing dust attenuation info

model = pipes.model_galaxy(model_components, filt_list=filter_list)

wlen = model.wavelengths
flux = model.spectrum_full

# Redshift the wavelength
wlen *= (1 + 7)
# Convert flux to correct units
flux = flux * (wlen**2) / (10**-29 * 2.9979 * 10**18) * 1e-19

plt.plot(wlen, flux)
plt.yscale('log')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Flux')
plt.legend()
plt.xlim(5000, 50000)
plt.yscale('log')
plt.ylim(3e-32, 1e-27)
plt.tight_layout()
plt.savefig('galaxy.png')
#plt.show()
