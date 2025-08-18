import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

r = np.linspace(0, 5, 500)

def sersic_n1(r, Reff):
    b_n = 1.678  # for n = 1
    Rs = Reff / b_n
    return np.exp(-r / Rs)

def gaussian(r, I0, sigma):
    return I0 * np.exp(-r**2 / (2 * sigma**2))

Reff = 1.0  
profile = sersic_n1(r, Reff)

# Fit only the core (e.g., r < 1.5)
fit_mask = r < 4
r_fit = r[fit_mask]
profile_fit = profile[fit_mask]

# Fit Gaussian to the core
popt, _ = curve_fit(gaussian, r_fit, profile_fit, p0=[1.0, 0.5])
I0_fit, sigma_fit = popt

# Compute FWHM of fitted Gaussian
fwhm_gaussian = 2 * np.sqrt(2 * np.log(2)) * sigma_fit
fwhm_ratio = fwhm_gaussian / Reff

# Plot
plt.figure(figsize=(8, 5))
plt.plot(r, profile, label="Sérsic n=1 (Reff=1)", lw=2)
plt.plot(r, gaussian(r, *popt), "--", label=f"Gaussian Fit\nFWHM ≈ {fwhm_gaussian:.2f} ({fwhm_ratio:.2f}×Reff)")
plt.axvline(fwhm_gaussian/2, color='gray', ls=':', label='FWHM/2')
plt.xlabel("Radius (r)")
plt.ylabel("Normalized Intensity")
plt.legend()
plt.title("Gaussian Fit to Core of Sérsic n=1 Profile")
plt.grid(True)
plt.tight_layout()
plt.show()
