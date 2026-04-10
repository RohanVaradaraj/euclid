
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.io import ascii

fontsize = 25

# Plotting configuration
plt.rcParams.update({
    'axes.linewidth': 2.5,
    'font.size': 15,
    'figure.dpi': 100
})

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    # 'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 2, 'ytick.major.width': 2,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 1.5, 'ytick.minor.width': 1.5,
})

# --- Functions to convert Muv + EW to Lya luminosity ---

def Muv_to_L_lambda(Muv, lam=1500.0):
    """
    Convert absolute UV magnitude Muv (AB) to UV continuum luminosity density L_lambda [erg/s/Å]
    lam in Angstrom
    """
    # L_nu at 10 pc (erg/s/Hz)
    L_nu = 10**(-0.4*(Muv + 48.6)) * 4 * np.pi * (10*3.0857e18)**2  # convert 10 pc -> cm

    # L_lambda = L_nu * c / lambda^2
    L_lambda = L_nu * 2.9979e18 / lam**2  # lam in Å, L_lambda in erg/s/Å

    return L_lambda

def Lya_from_EW(Muv, EW0):
    """
    Compute Lya luminosity from Muv and EW in erg/s
    """
    L_lambda = Muv_to_L_lambda(Muv)
    return EW0 * L_lambda  # Lya luminosity in erg/s


# --- Read SHELLQs data ---
t = Table.read('SHELLQs.dat.txt', format='ascii')
shellqs_Muv = t['col4']
shellqs_EW = t['col8']

# Read SILVERRUSH data
colnames = ['HSC', 'HSC_id', 'RAh', 'RAm', 'RAs', 'DE_sign', 'DEd', 'DEm', 'DEs',
            'zsp', 'gmag', 'rmag', 'imag', 'zmag', 'ymag', 'NBmag', 'NBSet', 'Feature', 'r_zsp']

col_starts = [0, 4, 19, 22, 25, 31, 32, 35, 38,
              43, 49, 54, 59, 64, 69, 74, 79, 85, 98]

col_ends   = [3, 18, 21, 24, 30, 32, 34, 37, 42,
              48, 53, 58, 63, 68, 73, 78, 84, 97, 126]  # exclusive

t_nb = Table.read('table4.dat.txt',
               format='ascii.fixed_width',
               col_starts=col_starts,
               col_ends=col_ends,
               names=colnames)

# Convert NB mag into Lya luminosity, assuming z=7.3
# Restrict to z=6.6 sources, so NB
t_nb_921 = t_nb[t_nb['NBSet'] == 'NB921']
t_nb_973 = t_nb[t_nb['NBSet'] == 'NB973']

mag_921 = t_nb_921['NBmag']
mag_973 = t_nb_973['NBmag']

# Convert to flux
flux_921 = 10**(-0.4*(mag_921 + 48.6))
flux_973 = 10**(-0.4*(mag_973 + 48.6))

# z=6.6 and z=7.3
z_921 = 6.6
z_973 = 7.3
d_L_921 = cosmo.luminosity_distance(z_921).to(u.cm).value
d_L_973 = cosmo.luminosity_distance(z_973).to(u.cm).value

# Now compute Lya luminosity assuming all NB flux is Lya
L_921 = flux_921 * 4 * np.pi * d_L_921**2 * (10**8)  # convert Hz to Å
L_973 = flux_973 * 4 * np.pi * d_L_973**2 * (10**8)  # convert Hz to Å

# For each, estimate Muv from zmag and ymag.
zmag_921 = t_nb_921['zmag']
ymag_973 = t_nb_973['ymag']

d_L_921 = cosmo.luminosity_distance(z_921).to(u.pc).value
d_L_973 = cosmo.luminosity_distance(z_973).to(u.pc).value
# Compute Muv
M_z = zmag_921 - 5*np.log10(d_L_921/10) + 2.5*np.log10(1+z_921)
M_y = ymag_973 - 5*np.log10(d_L_973/10) + 2.5*np.log10(1+z_973)
print(M_z)

# Read chorus csv file
#t_nb = ascii.read('860693.csv')
# t_nb = ascii.read('860694.csv')


# nb_mag = t_nb['n921_undeblended_apertureflux_20_mag']
# # Convert to flux
# nb_flux = 10**(-0.4*(nb_mag + 48.6))
# print(nb_flux)

# # Assume z=6.6, compute luminosity distance
# z_nb = 6.6
# d_L = cosmo.luminosity_distance(z_nb).to(u.cm).value

# # Now compute Lya luminosity assuming all NB flux is Lya
# nb_L = nb_flux * 4 * np.pi * d_L**2 * (10**8)  # convert Hz to Å

# print(nb_L)

# # Now estimate Muv from z flux
# z_mag = t_nb['z_undeblended_apertureflux_20_mag']  
# z_flux = 10**(-0.4*(z_mag + 48.6))
# d_L = cosmo.luminosity_distance(z_nb).to(u.pc).value
# # Compute Muv
# M_z = z_mag - 5*np.log10(d_L/10) + 2.5*np.log10(1+z_nb)
# print(M_z)

# Plot as contours in Muv, Lya


# plt.hist(np.log10(nb_L), bins=100)
# plt.show()
# exit()



# Compute Lya luminosities
shellqs_Lya = Lya_from_EW(shellqs_Muv, shellqs_EW)

# --- Endsley data ---
endsley_Muv = np.array([-21.8, -21.6, -22.7, -22.6, -21.6, -21.8, -21.9, -22.7])
endsley_EW = np.array([14.7, 3.8, 14.6, 3.7, 20, 5.4, 5.6, 10])
endsley_Lya = Lya_from_EW(endsley_Muv, endsley_EW)

# --- Plot ---
plt.figure(figsize=(10, 8))

# Plot SHELLQs and Endsley
plt.scatter(shellqs_Muv, shellqs_Lya, color='gray', label=r'Matsuoka+22 ($5.8<z<7$)')
plt.scatter(endsley_Muv, endsley_Lya, color='tab:blue', label=r'Endsley+22 ($z\simeq7$)', s=100, marker='s')

# Plot my sources
Muv_LAE1 = -20.8
EW_LAE1 = 240
Lya_LAE1 = 10**(43.7)
plt.scatter(Muv_LAE1, Lya_LAE1, color='tab:red', s=600, label='EUCL LAE 1', edgecolor='black', linewidth=2, marker='*')

# Vertical lines for Muv = -21.3, -21.4 since we have no Lya for LBG1 and LBG2 yet
plt.axvline(-21.35, color='purple', linestyle='--', label='EUCL LBG 1/2', linewidth=4)
#plt.axvline(-21.4, color='brown', linestyle='--', label='LBG2 Muv')

# Add two green hexagons for CR7 and Himiko
Muvs = [-21.8, -21.8]
Lyas = [10**43.6, 10**43.85]
plt.scatter(Muvs, Lyas, color='green', s=300, label='CR7 / Himiko', edgecolor='black', linewidth=2, marker='h')

# Plot SILVERRUSH sources
# plt.scatter(M_z, nb_L, color='green', label='SILVERRUSH Sources', alpha=0.5)
# 2D histogra

# plt.scatter(M_z, L_921, color='green', label='SILVERRUSH z=6.6 Sources', alpha=0.5)
# plt.scatter(M_y, L_973, color='orange', label='SILVERRUSH z=7.3 Sources', alpha=0.5)

# mask = np.isfinite(M_z) & np.isfinite(nb_L) & (nb_L > 0)
# M_z_clean = M_z[mask]
# nb_L_clean = nb_L[mask]

# Find top 20 by luminosity
# top20_idx = np.argsort(nb_L_clean)[-20:]

# counts, xedges, yedges = np.histogram2d(M_z_clean , nb_L_clean, bins=100)
# X, Y = np.meshgrid(0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1]))

# plt.contour(X, Y, counts.T, levels=10, cmap='viridis')
# plt.scatter(M_z_clean[top20_idx], nb_L_clean[top20_idx], color='red', label='Top 20 brightest')



# Overlay lines of constant EW
Muvs = np.linspace(-26.2, -20, 100)
EWs_const = [10, 100, 1000]  # constant EW in Å
for ew_const in EWs_const:
    Ls = Lya_from_EW(Muvs, ew_const)
    plt.plot(Muvs, Ls, linestyle='--', color='black', linewidth=2)
    # Add text label, bold
    plt.text(-21.9, Lya_from_EW(-21.9, ew_const)*1.2, f"EW={ew_const} Å", rotation=32, color='black', fontsize=20, fontweight='bold')

plt.yscale('log')
plt.xlabel(r'$M_{\rm UV}$', fontsize=fontsize)
plt.ylabel(r'$\log_{10}(L_{\mathrm{Ly}\alpha} \, / \mathrm{erg \, s^{-1}})$', fontsize=fontsize)
plt.gca().invert_xaxis()
plt.ylim(2e42, 2e45)
plt.xlim(-20.5, -26.2)
plt.legend(loc='lower right', fontsize=18)
# Match axis tick labels to fontsize
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.tight_layout()
plt.savefig('Muv_Lya_luminosity.pdf', dpi=100)
plt.show()

