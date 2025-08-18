"""
Check the photometry of stars with a flat spectrum

Created: Tuesday 23rd July 2024.
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 3, 'ytick.major.width': 3,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 2, 'ytick.minor.width': 2,
})

sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import filter_widths

def flux_to_mag(flux):

    mag = -2.5*np.log10(flux) -48.6
    return mag

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 25})
plt.rcParams['figure.dpi'] = 100

cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'

# Read in the photometry
t = Table.read(cat_dir / 'det_YJH_flat_stars_VISTA_Y_with_JWST.fits')
#t = Table.read(cat_dir / 'det_YJH_stars_VISTA_Y.fits')
#t = Table.read(cat_dir / 'det_YJH_flat_stars_HSC_I.fits')

# Restrict JWST table to where error is > 0
t = t[t['err_f115w'] > 0]


filter_info = filter_widths()


Y = flux_to_mag(t['flux_Y'])
J = flux_to_mag(t['flux_J'])
H = flux_to_mag(t['flux_H'])
r = flux_to_mag(t['flux_HSC-R_DR3'])
i = flux_to_mag(t['flux_HSC-I_DR3'])
z = flux_to_mag(t['flux_HSC-Z_DR3'])
y = flux_to_mag(t['flux_HSC-Y_DR3'])

f115w = flux_to_mag(t['flux_f115w'])
f150w = flux_to_mag(t['flux_f150w'])
f277w = flux_to_mag(t['flux_f277w'])
f444w = flux_to_mag(t['flux_f444w'])


#! UNCOMMENT TO PLOT STELLAR LOCUS
#plt.scatter(r-i, i-z, color='black', s=5)
# plt.scatter(f115w-f150w, f150w-f277w, color='red', s=5)
# plt.show()
# exit()



#! UNCOMMENT TO PLOT PHOTOMETRY/SED
# for obj in t[0:10]:
#     plt.figure(figsize=(10, 6))
#     for filt in filter_info.keys():
#         # Check if name is vis
#         if filt == 'VIS':
#             plt.errorbar(filter_info[filt][0], flux_to_mag(obj[f'flux_{filt}']), yerr=np.abs((2.5 / np.log(10)) * (obj[f'err_{filt}'] / obj[f'flux_{filt}'])), fmt='o', color='red', markersize=5)
#         else:
#             plt.errorbar(filter_info[filt][0], flux_to_mag(obj[f'flux_{filt}']), yerr=np.abs((2.5 / np.log(10)) * (obj[f'err_{filt}'] / obj[f'flux_{filt}'])), fmt='o', color='black', markersize=5)
#     plt.yscale('log')
#     plt.gca().invert_yaxis()
#     plt.show()


# exit()

#! UNCOMMENT TO SAVE FLAT SED STARS TO A TABLE
#Select stars flat in i-z and z-y
# mask = (np.abs(r-i) < 0.2) & (np.abs(i-z) < 0.2)
# t = t[mask]
# t.write(cat_dir / 'det_YJH_flat_stars_HSC_I.fits', overwrite=True)
# exit()

Ye = flux_to_mag(t['flux_Ye'])
Je = flux_to_mag(t['flux_Je'])
He = flux_to_mag(t['flux_He'])
VIS = flux_to_mag(t['flux_VIS'])

# Errors
dY = (2.5 / np.log(10)) * (t['err_Y'] / t['flux_Y'])
dJ = (2.5 / np.log(10)) * (t['err_J'] / t['flux_J'])
dH = (2.5 / np.log(10)) * (t['err_H'] / t['flux_H'])
di = (2.5 / np.log(10)) * (t['err_HSC-I_DR3'] / t['flux_HSC-I_DR3'])

dYe = (2.5 / np.log(10)) * (t['err_Ye'] / t['flux_Ye'])
dJe = (2.5 / np.log(10)) * (t['err_Je'] / t['flux_Je'])
dHe = (2.5 / np.log(10)) * (t['err_He'] / t['flux_He'])
dVIS = (2.5 / np.log(10)) * (t['err_VIS'] / t['flux_VIS'])

df115w = (2.5 / np.log(10)) * np.abs((t['err_f115w'] / t['flux_f115w']))
df150w = (2.5 / np.log(10)) * (t['err_f150w'] / t['flux_f150w'])
df277w = (2.5 / np.log(10)) * (t['err_f277w'] / t['flux_f277w'])
df444w = (2.5 / np.log(10)) * (t['err_f444w'] / t['flux_f444w'])

# Propagate errors in subtraction
dYYe = np.sqrt(dY**2 + dYe**2)
dJJe = np.sqrt(dJ**2 + dJe**2)
dHHe = np.sqrt(dH**2 + dHe**2)
dVISi = np.sqrt(dVIS**2 + dVIS**2)

dYef115w = np.sqrt(dYe**2 + df115w**2)
dJef150w = np.sqrt(dJe**2 + df150w**2)



# print(dY)


plt.figure(figsize=(10, 6))
plt.errorbar(Y, Y-Ye, yerr=dYYe, xerr=dY, fmt='o', color='black', markersize=5, alpha=0.9)
plt.plot([min(Y), max(Y)], [-0.1, -0.1], color='red', linestyle='dotted')
plt.plot([min(Y), max(Y)], [0, 0], color='red', linestyle='--')
plt.plot([min(Y), max(Y)], [0.1, 0.1], color='red', linestyle='dotted')
#plt.text(17, 0.4, r'$| Y - J| < 0.05 \wedge |J - H| < 0.05$')
plt.xlabel(r'$Y$', fontsize=25)
plt.ylabel(r'$Y - Y_{\rm{E}}$', fontsize=25)
plt.ylim(-0.5, 0.5)
plt.xlim(17, 21.2)
#plt.tick_params(axis='both', which='major', width=2.5, length=5)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(Path.cwd().parents[1] / 'plots' / 'psf' / 'Y_Ye.pdf', bbox_inches='tight')
#plt.show()
plt.close()

plt.figure(figsize=(10, 6))
plt.errorbar(J, J-Je, yerr=dJJe, xerr=dJ, fmt='o', color='black', markersize=5, alpha=0.9)
plt.plot([min(J), max(J)], [-0.1, -0.1], color='red', linestyle='dotted')
plt.plot([min(J), max(J)], [0, 0], color='red', linestyle='--')
plt.plot([min(J), max(J)], [0.1, 0.1], color='red', linestyle='dotted')
plt.xlabel(r'$J$', fontsize=25)
plt.ylabel(r'$J - J_{\rm{E}}$', fontsize=25)
plt.ylim(-0.5, 0.5)
plt.xlim(17, 21.2)
#plt.tick_params(axis='both', which='major', width=2.5, length=5)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(Path.cwd().parents[1] / 'plots' / 'psf' / 'J_Je.pdf', bbox_inches='tight')
#plt.show()
plt.close()

plt.figure(figsize=(10, 6))

# Remove H-He > 0.1 or H-He < -0.3
mask = (H-He < 0.1) & (H-He > -0.3)
H = H[mask]
He = He[mask]
dH = dH[mask]
dHe = dHe[mask]
dHHe = dHHe[mask]

plt.errorbar(H, H-He, yerr=dHHe, xerr=dH, fmt='o', color='black', markersize=5)
plt.plot([min(H), max(H)], [-0.1, -0.1], color='red', linestyle='dotted')
plt.plot([min(H), max(H)], [0, 0], color='red', linestyle='--')
plt.plot([min(H), max(H)], [0.1, 0.1], color='red', linestyle='dotted')
#plt.title(r'$| Y - J| < 0.05 \wedge |J - H| < 0.05$')
plt.xlabel(r'$H$', fontsize=25)
plt.ylabel(r'$H - H_{\rm{E}}$', fontsize=25)
plt.ylim(-0.5, 0.5)
plt.xlim(17, 21.2)
#plt.tick_params(axis='both', which='major', width=2.5, length=5)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(Path.cwd().parents[1] / 'plots' / 'psf' / 'H_He.pdf', bbox_inches='tight')
#plt.show()
plt.close()

# plt.figure(figsize=(10, 6))
# plt.errorbar(i, i-VIS, yerr=dVISi, xerr=di, fmt='o', color='black', markersize=5)
# plt.plot([min(i), max(i)], [-0.1, -0.1], color='red', linestyle='dotted')
# plt.plot([min(i), max(i)], [0, 0], color='red', linestyle='--')
# plt.plot([min(i), max(i)], [0.1, 0.1], color='red', linestyle='dotted')
# plt.title(r'$| r - i| < 0.2 \wedge |i - z| < 0.2$')
# plt.xlabel(r'$i_{\mathrm{HSC}}$')
# plt.ylabel(r'$i_{\mathrm{HSC}} - i_{E}$')
# plt.ylim(-0.5, 0.5)
# plt.show()

# plt.figure(figsize=(10, 6))
# plt.errorbar(f115w, f115w-Y, yerr=dYef115w, xerr=np.abs(df115w), fmt='o', color='black', markersize=5)
# plt.plot([min(f115w), max(f115w)], [-0.1, -0.1], color='red', linestyle='dotted')
# plt.plot([min(f115w), max(f115w)], [0, 0], color='red', linestyle='--')
# plt.plot([min(f115w), max(f115w)], [0.1, 0.1], color='red', linestyle='dotted')
# plt.title(r'$| Y - J| < 0.05 \wedge |J - H| < 0.05$')
# plt.xlabel(r'$f115w$')
# plt.ylabel(r'$f115w - Y$')
# plt.ylim(-0.5, 0.5)
# plt.savefig(Path.cwd().parents[1] / 'plots' / 'psf' / 'f115w_Y.pdf')
# plt.show()
