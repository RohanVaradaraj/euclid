import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import os
import glob
import matplotlib.pyplot as plt
from pathlib import Path
import astropy.units as u
from tqdm import tqdm
from scipy.spatial import cKDTree

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 3
plt.rcParams['figure.dpi'] = 100

cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits'
t = Table.read(cat_dir / cat_name)

# if there isnt a completeness column in table, make it
if 'completeness' not in t.colnames:
    t['completeness'] = np.zeros(len(t))

def crossmatch_and_compute_completeness(input_cat_dir, output_cat_dir, Muv_bins, z_bins, match_radius_pix=2):
    """
    Computes the completeness (recovery rate) for each bin of Muv and z using RA and DEC.
    Optimized using cKDTree for faster spatial matching.
    """
    completeness_matrix = np.zeros((len(Muv_bins) - 1, len(z_bins) - 1))
    total_sources_per_bin = np.zeros_like(completeness_matrix)
    recovered_sources_per_bin = np.zeros_like(completeness_matrix)
    
    input_cats = glob.glob(os.path.join(input_cat_dir, '*.fits'))
    input_cat_names = [os.path.basename(cat) for cat in input_cats]

    for input_cat_name in tqdm(input_cat_names, desc="Processing input catalogs"):
        input_cat = Table.read(os.path.join(input_cat_dir, input_cat_name))
        output_cat_name = input_cat_name.split('_input_values.fits')[0] + '_cat.fits'

        # Check if output catalog exists, if not then skip
        if not os.path.exists(os.path.join(output_cat_dir, output_cat_name)):
            print(f"Output catalog {output_cat_name} not found. Skipping.")
            continue
        else:
            output_cat = Table.read(os.path.join(output_cat_dir, output_cat_name))

            # Positions from input and output catalogs
            input_positions = np.column_stack((input_cat['x'], input_cat['y']))
            output_positions = np.column_stack((output_cat['X_IMAGE'], output_cat['Y_IMAGE']))

            # Build KDTree for output positions
            output_tree = cKDTree(output_positions)

            # Query the tree for matches within the match radius
            distances, indices = output_tree.query(input_positions, distance_upper_bound=match_radius_pix)

            # Identify matched sources (valid distances are less than the match radius)
            matched = distances < match_radius_pix

            # Binning based on Muv and z
            Muv_values = input_cat['Muv']
            z_values = input_cat['z']

            for i in range(len(Muv_bins) - 1):
                for j in range(len(z_bins) - 1):
                    in_Muv_bin = (Muv_values >= Muv_bins[i]) & (Muv_values < Muv_bins[i + 1])
                    in_z_bin = (z_values >= z_bins[j]) & (z_values < z_bins[j + 1])
                    selected_input_sources = in_Muv_bin & in_z_bin

                    total_sources_per_bin[i, j] += np.sum(selected_input_sources)
                    recovered_sources_per_bin[i, j] += np.sum(matched[selected_input_sources])

    with np.errstate(divide='ignore', invalid='ignore'):
        completeness_matrix = recovered_sources_per_bin / total_sources_per_bin
        completeness_matrix[np.isnan(completeness_matrix)] = 0

    # Flip matrix in muv direction
    completeness_matrix = np.flip(completeness_matrix, axis=0)

    return completeness_matrix



def plot_completeness(completeness_matrix, Muv_bins, z_bins):
    """
    Plots the completeness matrix as a 2D heatmap.
    
    Parameters:
        completeness_matrix (2D array): The completeness values for each Muv, z bin.
        Muv_bins (array-like): The Muv bin edges.
        z_bins (array-like): The z bin edges.
    """
    plt.imshow(completeness_matrix, aspect='auto', origin='lower', cmap='viridis',
               extent=[z_bins[0], z_bins[-1], Muv_bins[-1], Muv_bins[0]])  # Invert Muv axis

    plt.tight_layout()    
    plt.colorbar()
    plt.xlabel(r'$z$')
    plt.ylabel(r'$M_{\rm UV}$')
    #plt.show()

def main():
    input_cat_dir =  Path.cwd() / 'catalogues' / 'input'
    output_cat_dir = Path.cwd() / 'catalogues' / 'output'
    
    # Define Muv and z bins
    Muv_bins = np.arange(-23, -20, 0.1)
    z_bins = np.arange(6.5, 7.5, 0.05)
    
    # Compute completeness
    #completeness_matrix = crossmatch_and_compute_completeness(input_cat_dir, output_cat_dir, Muv_bins, z_bins)

    # Read in completeness matrix
    completeness_matrix = np.load('completeness_matrix_2.npy')

    # Flip matrix in Muv direction
    #completeness_matrix = np.flip(completeness_matrix, axis=0)

    # for i, obj in enumerate(t):

    #     Muv = obj['Muv']
    #     z = obj['Zphot']

    #     print(f'z: {z}, Muv: {Muv}')

    #     Muv_bin = np.digitize([Muv], Muv_bins)[0] - 1
    #     z_bin = np.digitize([z], z_bins)[0] - 1

    #     # If z>7.49 and z<7.5, set z_bin to 18
    #     if z_bin >18:
    #         z_bin = 18

    #     print(f'Muv bin: {Muv_bin}, z bin: {z_bin}')

    #     if z > 7.495:
    #         completeness = 0
    #     else:
    #         completeness = completeness_matrix[z_bin, Muv_bin]
    #     print(f'Completeness: {completeness}')

    #     t['completeness'][i] = completeness

    # t.write(cat_dir / cat_name, overwrite=True)
        
    # Plot the completeness
    plot_completeness(completeness_matrix, Muv_bins, z_bins)

    # Save the completeness matrix as npy file
    #np.save('completeness_matrix_2.npy', completeness_matrix)


if __name__ == "__main__":
    main()
    plot_dir = Path.cwd().parents[1] / 'plots' / 'completeness'
    #plot_dir.mkdir(exist_ok=True)
    #plt.savefig(plot_dir / 'completeness_heatmap_2.pdf')
    plt.show()