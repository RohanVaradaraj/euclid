import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import os
import glob
import matplotlib.pyplot as plt
from pathlib import Path
import astropy.units as u

def crossmatch_and_compute_completeness(input_cat_dir, output_cat_dir, Muv_bins, z_bins, match_radius_pix=6.67):
    """
    Computes the completeness (recovery rate) for each bin of Muv and z using RA and DEC.
    """
    # Initialize completeness matrix, with z on x axis, Muv on y axis
    completeness_matrix = np.zeros((len(Muv_bins) - 1, len(z_bins) - 1))
    
    # Initialize counters for total and recovered sources per bin
    total_sources_per_bin = np.zeros_like(completeness_matrix)
    recovered_sources_per_bin = np.zeros_like(completeness_matrix)
    
    # List all input catalogs (which should match output catalogs by name)
    input_cats = glob.glob(os.path.join(input_cat_dir, '*.fits'))
    input_cat_names = [os.path.basename(cat) for cat in input_cats]

    for input_cat_name in input_cat_names:
        # Read catalogs
        input_cat = Table.read(os.path.join(input_cat_dir, input_cat_name))
        output_cat_name = input_cat_name.split('_input_values.fits')[0] + '_cat.fits'
        output_cat = Table.read(os.path.join(output_cat_dir, output_cat_name))

        # Extract coordinates
        input_x, input_y = input_cat['x'], input_cat['y']
        output_x, output_y = output_cat['X_IMAGE'], output_cat['Y_IMAGE']

        # Compute distances
        dist_x = input_x[:, None] - output_x
        dist_y = input_y[:, None] - output_y
        d2 = np.sqrt(dist_x**2 + dist_y**2)
        matched = np.any(d2 < match_radius_pix, axis=1)  # Check if each input source has at least one match

        # Bin assignment
        Muv_values = input_cat['Muv']
        z_values = input_cat['z']

        for i in range(len(Muv_bins) - 1):
            for j in range(len(z_bins) - 1):
                in_Muv_bin = (Muv_values >= Muv_bins[i]) & (Muv_values < Muv_bins[i+1])
                in_z_bin = (z_values >= z_bins[j]) & (z_values < z_bins[j+1])
                selected_input_sources = in_Muv_bin & in_z_bin

                # Update counters for this bin
                total_sources_per_bin[i, j] += np.sum(selected_input_sources)
                recovered_sources_per_bin[i, j] += np.sum(matched[selected_input_sources])

    # Compute completeness for each bin
    with np.errstate(divide='ignore', invalid='ignore'):  # Avoid division errors
        completeness_matrix = recovered_sources_per_bin / total_sources_per_bin
        completeness_matrix[np.isnan(completeness_matrix)] = 0  # Replace NaN values with 0

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
    #plt.gca().invert_yaxis()  # Make brighter (negative) Muv values appear at the top
    plt.colorbar(label='Completeness')
    plt.xlabel('Redshift (z)')
    plt.ylabel('Muv')
    plt.title('Source Injection Recovery Completeness')
    plt.show()

def main():
    input_cat_dir =  Path.cwd() / 'catalogues' / 'input'
    output_cat_dir = Path.cwd() / 'catalogues' / 'output'
    
    # Define Muv and z bins
    Muv_bins = np.arange(-25, -20, 0.1)
    z_bins = np.arange(6.5, 7.5, 0.05)
    
    # Compute completeness
    completeness_matrix = crossmatch_and_compute_completeness(input_cat_dir, output_cat_dir, Muv_bins, z_bins)
    
    # Plot the completeness
    plot_completeness(completeness_matrix, Muv_bins, z_bins)


if __name__ == "__main__":
    main()
