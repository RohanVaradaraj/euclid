"""
For every candidate galaxy, check the redshifts of all galaxies within a 1 arcminute radius.

We are looking for peaks in the redshift distribution at low-redshift (z~1).

Created: Friday 6th June 2025. First new script post-viva!
"""

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

#! Field name
field_name = 'COSMOS'

#! Dictionary to map field name into Natalia's directory structure
field_info = {
    'XMM': {
        'subdir': 'XMM-LSS',
        'filename': 'XMM_DR3_MASKED_Ks-stellar-cut_peters-LP-photo-z_LPz-AGN_chi2_LPclass_ModBest_LPmasses_LPsfr_LPssfr_flag_25-05-2025.fits',
    },
    'COSMOS': {
        'subdir': 'COSMOS_DR6',
        'filename': 'COSMOS_DR6_HSC-DR3_MASKED_Ks-stellar-cut_peters-LP-photo-z_LPz-AGN_chi2_LPclass_ModBest_LPmasses_LPsfr_LPssfr_flag_25-05-2025.fits',
    },
    'CDFS': {
        'subdir': 'E-CDFS',
        'filename': 'ECDFS_MASKED_Ks-stellar-cut_LP-photo-z_LPz-AGN_chi2_LPclass_ModBest_LPmasses_LPsfr_LPssfr_flag_Spec-z_25-05-2025.fits',
    },
}

#! Read in master redshift catalogue, from Natalia.
natalia_dir = Path.cwd().parents[3] / 'natalia' / 'photo-z_catalogues'
zphot_dir = natalia_dir / field_info[field_name]['subdir']
cat_file = zphot_dir / field_info[field_name]['filename']

t = Table.read(cat_file)
print(t)
print(t.colnames)


#! Read in my candidate galaxies
candidate_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
candidate_dict = {
    'XMM': 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits',
    'COSMOS': 'COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06.fits',
    'CDFS': 'CDFS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits',
}

cands = Table.read(candidate_dir / candidate_dict[field_name])
cands.sort('Muv')
print(cands)
print(cands.colnames)

ID_to_check = 1258338  # Change this to the ID of the candidate you want to check

input_file = Path.cwd() / f'low_z_with_prior_IDs_{field_name}.txt'

with open(input_file, 'r') as f:
    low_z_ids = [int(line.strip()) for line in f]
low_z_ids = np.array([int(z) for z in low_z_ids])

print(f"Low redshift IDs from file: {low_z_ids}")

# Restrict candidates to those with IDs in the low_z_ids list
mask = np.isin(cands['ID'], low_z_ids)
cands = cands[mask]
print(cands)

#! Function to find all galaxies in Nat's catalogue within 1 arcminute
def find_nearby_galaxies(candidate, t, radius=1 * u.arcminute):
    """
    Find galaxies in Natalia's catalogue within a specified radius of the candidate galaxy.
    
    Parameters:
    candidate (TableRow): A row from the candidates table.
    t (Table): The full redshift catalogue.
    radius (Quantity): The search radius in angular units.
    
    Returns:
    Table: A table of nearby galaxies.
    """
    coord = SkyCoord(ra=candidate['RA'], dec=candidate['DEC'], unit=(u.deg, u.deg))
    nearby = t[(SkyCoord(ra=t['RA'], dec=t['DEC'], unit=(u.deg, u.deg)).separation(coord) < radius)]
    return nearby

# Loop through each candidate galaxy and find nearby galaxies
for i, candidate in enumerate(cands):
    print(f"Processing candidate {i + 1}/{len(cands)}: {candidate['ID']}")
    nearby_galaxies = find_nearby_galaxies(candidate, t, radius=2 * u.arcminute)

    # Plot the RA DEC of the candidate and nearby galaxies
    # plt.figure(figsize=(8, 8))
    # plt.scatter(candidate['RA'], candidate['DEC'], color='red', label='Candidate', s=100, edgecolor='black')
    # plt.scatter(nearby_galaxies['RA'], nearby_galaxies['DEC'], 
    #             color='blue', label='Nearby Galaxies', s=10, alpha=0.5)
    # plt.xlabel('RA (degrees)')
    # plt.ylabel('DEC (degrees)')
    # plt.show()
    # plt.close()
    
    if len(nearby_galaxies) > 0:
        print(f"Found {len(nearby_galaxies)} nearby galaxies for candidate {candidate['ID']}")

        # This galaxy primary and secondary redshift
        z_lbg = candidate['Zphot']
        z_dusty = candidate['Zphot_sec']

        # Plot histogram of their redshifts
        z_nearby = nearby_galaxies['Z_BEST_peak']
        plt.figure(figsize=(10, 6))
        plt.hist(z_nearby, bins=np.arange(0,2,0.05), alpha=0.7, color='tab:red', histtype='step', lw=3, density=True, label='Within 1 arcminute')

        # Plot also the hist of all galaxies in this range
        z_all = t['Z_BEST_peak']
        plt.hist(z_all, bins=np.arange(0,2,0.01), alpha=0.5, color='gray', histtype='step', lw=3, density=True, label='All Galaxies')

        plt.xlabel(r'$z$')
        plt.ylabel('Number of Galaxies')
        plt.title(f'ID: {candidate["ID"]}' +  r', $M_{\rm{UV}}=$' + f'{candidate["Muv"]:.2f}' +  r', $z_{\rm{LBG}}=$' + f'{z_lbg:.2f}' +  r', $z_{\rm{dusty}}=$' + f'{z_dusty:.2f}')
        plt.axvline(z_dusty, color='green', linestyle='--', label=r'$z_{\rm{dusty}}$')
        plt.legend()
        plt.show()

    else:
        print(f"No nearby galaxies found for candidate {candidate['ID']}")


