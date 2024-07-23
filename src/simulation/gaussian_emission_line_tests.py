from astropy import units as u
from synphot import SourceSpectrum
from synphot.models import GaussianFlux1D, ConstFlux1D
import matplotlib.pyplot as plt
from synphot import etau_madau


# Create a flat source
sp = SourceSpectrum(ConstFlux1D, amplitude=1e-29)



# Apply extinction for a given redshift
z = 2
wave = range(2400, 4200)  # Angstrom
extcurve = etau_madau(wave, z)
sp *= extcurve

# Add emission line
emission = SourceSpectrum(GaussianFlux1D, mean=3750*u.angstrom, fwhm=10*u.angstrom,
                    total_flux= (1e-29*u.erg/(u.cm**2 * u.s * u.angstrom) * 10 * u.angstrom) - 1e-29*u.erg/(u.cm**2 * u.s))

sp += emission

plt.plot(wave, sp(wave), 'b--')
plt.show()

area = sp.integrate()

# Get equivalent width = area / continuum flux
