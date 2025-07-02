"""
Experimenting with placing a luminosity function prior on the P(z) of an SED.

Just trying on one object or a handful of objects for now.

Created: Wednesday 14th May 2025.
"""
from pathlib import Path
import sys
sed_fitting_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_fitting_path))

from sed_fitting_codes import parse_spec_file
import matplotlib.pyplot as plt
import glob
import numpy as np
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from photoz_catalogue_priors import prior_func

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

#! SWITCH TO PLOT P(z) FOR ALL OBJECTS, TO TAKE A LOOK
plot_all_pdfs = False
plot = False

#! INTERESTING IDs
# IDs = [376795, 905849, 1053436]
# IDs = [460143]

#! Experimenting on an individual ID
# ID = IDs[0]
# #ID = 788054 # Good high-z
# ID = 376795 # Secondary low-z

# ID = 688921

sed_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'XMM' / 'best_fits' / 'det_HSC-Z_DR3_z7'
spec_files = glob.glob(str(sed_dir / '*.spec'))

#! Define a cosmological model
H = 70
omegaM = 0.3
omegaV = 0.7
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)

#! Bouwens+21 LF parameter evolution
def alpha_evolution(z):
    """
    Evolution of the faint end slope of the LF.
    """
    return -1.94 - 0.11 * (z - 6)

def phi_star_evolution(z):
    """
    Evolution of the characteristic value of the LF.
    """
    exponent = -0.33 * (z - 6) - 0.024 * (z - 6)**2
    return 0.4e-3 * 10 ** exponent

def Mstar_evolution(z):
    """
    Evolution of the characteristic magnitude of the LF.
    """
    z_t = 2.46
    return np.where(
        z < z_t,
        -20.89 - 1.09 * (z - z_t),
        -21.03 - 0.04 * (z - z_t)
    )


#! Schechter function from the parameters
def Schecter_function(M, alpha, phi_star, M_star):
    """
    Produces Schechter function given parameters: normalisation phi*, faint slope alpha, mag array M and char. mag M*. 
    """

    coeff = np.log(10) / 2.5

    faint = (10 ** (0.4 * (M_star - M))) ** (alpha+1)


    bright_exponent = -10 ** (0.4 * (M_star - M))
    bright = np.exp(bright_exponent)

    phi = coeff * phi_star * faint * bright

    return phi



#! Open the parent catalogue ID
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
parent_cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits'
parent_cat = Table.read(cat_dir / parent_cat_name)

parent_cat.sort('Muv', reverse=False)

# Array to store IDs which have low redshift when the prior is applied
low_z_ids = []

for i, ID in enumerate(parent_cat['ID']):
    print(f'Processing ID {i+1}/{len(parent_cat)}: {ID}')

    # Get Muv of this ID
    # if ID != 688921:
    #     continue
    ID_row = parent_cat['ID'] == ID
    Muv = parent_cat['Muv'][ID_row][0]
    m_Z = -2.5*np.log10(parent_cat['flux_HSC-Z_DR3'][ID_row][0])-48.6
    print('Apparent magnitude in HSC-Z:', m_Z)

    #! Get .spec file
    # Pad ID with zeros to left so that there are nine digits
    ID_str = str(ID).zfill(9)
    spec_file = f'{sed_dir}/Id{ID_str}.spec'
    file = parse_spec_file(spec_file)

    phot = file.get('phot')
    zpdf = file.get('zpdf')
    sed = file.get('sed')
    params = file.get('model')

    #! Open the prior file, prior_fits.txt
    prior_file = 'prior_fits.txt'
    prior_table = Table.read(prior_file, format='ascii.commented_header', names=['field_name', 'mag', 'C', 'gamma', 'z_0', 'beta'], comment='#')

    # Go through the table until m_Z < mag, then take the prevous one
    for i in range(len(prior_table)):
        if m_Z < prior_table['mag'][i]:
            break

    # Get the parameters for this bin
    C = prior_table['C'][i-1]
    gamma = prior_table['gamma'][i-1]
    z_0 = prior_table['z_0'][i-1]
    beta = prior_table['beta'][i-1]
    print(f'Parameters for m_Z < {m_Z}: C = {C}, gamma = {gamma}, z_0 = {z_0}, beta = {beta}')

    #! Compute the prior on the same redshift grid as the zpdf
    prior_vals = prior_func(zpdf['z'], C, gamma, z_0, beta)

    #! Normalise so that the integral is 1
    prior_vals /= np.trapz(prior_vals, zpdf['z'])

    #! Compute the comoving volume at each redshift in z_grid increment
    comoving_volume = cosmo.comoving_volume(zpdf['z']).value

    #! Compute LF * comoving volume as function of redshift
    M_values = np.arange(-23, -19, 0.1)
    cmap = cm.coolwarm
    norm = mcolors.Normalize(vmin=min(M_values), vmax=max(M_values))
    print('Original integral:', np.trapz(zpdf['P(z)'], zpdf['z']))

    pdf_with_prior = zpdf['P(z)'] * prior_vals

    # Normalize the PDF with prior, such that the integral is 1
    integral = np.trapz(pdf_with_prior, zpdf['z'])
    pdf_with_prior /= np.max(pdf_with_prior)

    # Find the peak of the PDF with prior
    peak_index = np.argmax(pdf_with_prior)
    peak_z = zpdf['z'][peak_index]
    print('Peak z:', zpdf['z'][peak_index])

    if peak_z < 3:
        low_z_ids.append(ID)

    # Find the 68th percentile of the data left and right of the peak
    left_mask = zpdf['z'] < peak_z
    right_mask = zpdf['z'] > peak_z

    left_z = zpdf['z'][left_mask]
    left_pdf = pdf_with_prior[left_mask]
    right_z = zpdf['z'][right_mask]
    right_pdf = pdf_with_prior[right_mask]

    # Get highest-posterior density interval (HPDI) for the left and right side
    left_cumulative = np.cumsum(left_pdf)
    right_cumulative = np.cumsum(right_pdf)
    left_cumulative /= np.max(left_cumulative)
    right_cumulative /= np.max(right_cumulative)

    # Find the 68th percentile
    left_68 = np.where(left_cumulative >= 0.68)[0][0]
    right_68 = np.where(right_cumulative >= 0.68)[0][0]
    left_68_z = left_z[left_68]
    right_68_z = right_z[right_68]

    # Draw vertical lines for the 68th percentile
    # plt.axvline(x=left_68_z, color='green', linestyle='--', label=r'$z_{\mathrm{68}}=$' + f'{left_68_z:.2f}', alpha=0.7, zorder=3)
    # plt.axvline(x=right_68_z, color='green', linestyle='--', alpha=0.7, zorder=3)

    # Draw a vertical line here


    # Find peak of PDF without prior
    peak_index_no_prior = np.argmax(zpdf['P(z)'])
    peak_z_no_prior = zpdf['z'][peak_index_no_prior]


    #! Plot the prior
    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(zpdf['z'], zpdf['P(z)'], color='black', linewidth=2, label='Flat prior', ls='--', alpha=1, zorder=2)
        plt.plot(zpdf['z'], prior_vals, color='gray', linewidth=3, label=r'$P(z | m_{\mathrm{AB}})$', alpha=0.7, zorder=0)
        plt.axvline(x=peak_z, color='red', linestyle='--', label=r'$z_{\mathrm{vol}}=$' + f'{peak_z:.2f}', alpha=0.7, zorder=3)
        plt.axvline(x=peak_z_no_prior, color='black', linestyle='--', label=r'$z_{\mathrm{flat}}=$' + f'{peak_z_no_prior:.2f}', alpha=0.7, zorder=3)
        plt.plot(zpdf['z'], pdf_with_prior, color='red', linewidth=3, label='Volume prior', alpha=1, zorder=1)
        plt.xlabel('z')
        plt.ylabel(r'$P(z | F) \propto e^{-\chi^2(z)} P(z | m_{\mathrm{AB}})$')
        plt.title(f'ID: {ID}, ' + r'$M_{\mathrm{UV}} = $' +f'{Muv:.2f}, ' + r'$m_{\mathrm{AB}} = $' + f'{m_Z:.2f}', pad=5)
        plt.legend()
        plt.savefig(f'ID_{ID}_prior.pdf', bbox_inches='tight')
        plt.show()

# Save low-z IDs to a file
output_file = Path.cwd() / f'low_z_with_prior_IDs.txt'
with open(output_file, 'w') as f:
    for id_val in low_z_ids:
        f.write(f'{id_val}\n')

print(f"Saved {len(low_z_ids)} IDs with peak_z < 3 to {output_file}")


