#!/bin/bash

# Define the filter names and tile suffixes
filters=("f115w" "f150w" "f277w" "f444w")
tiles=("0A" "0B" "1A" "1B" "2A" "2B" "3A" "3B" "4A" "4B" "5A" "5B" "6A" "6B" "7A" "7B")

# Paths for psfex and input/output directories
psfex_path="/mnt/zfsusers/varadaraj/psfex/bin/psfex"
input_dir="/mnt/vardy/vardygroupshare/rohan/euclid/data/psf/COSMOS/catalogues"
config_file="/mnt/vardy/vardygroupshare/HSC_SSP_DR3/config_files/default_cube.psfex"
output_dir="/mnt/vardy/vardygroupshare/rohan/euclid/data/psf/COSMOS/results"
psf_size="75,75"
psf_var_degrees=5
psf_var_nsnap=1
psf_sampling=1.0
fwhm_range="1,10"

# Loop through each filter and tile
for filter in "${filters[@]}"; do
    for tile in "${tiles[@]}"; do
        # Construct input catalog filename and output filenames
        input_file="${input_dir}/${filter}_${tile}.fits"
        snap_output="${output_dir}/snap_${filter}_${tile}.fits"
        samples_output="${output_dir}/samples_${filter}_${tile}.fits"
        residuals_output="${output_dir}/resi_${filter}_${tile}.fits"

        # Run the psfex command with the appropriate arguments
        $psfex_path $input_file -c $config_file -PSF_DIR $output_dir -PSF_SIZE $psf_size \
        -PSFVAR_DEGREES $psf_var_degrees -PSFVAR_NSNAP $psf_var_nsnap -CHECKIMAGE_TYPE SNAPSHOTS,SAMPLES,RESIDUALS \
        -CHECKIMAGE_NAME $snap_output,$samples_output,$residuals_output -SAMPLE_AUTOSELECT N -PSF_SAMPLING $psf_sampling \
        -SAMPLE_FWHMRANGE $fwhm_range

        echo "Processed $filter $tile"
    done
done

echo "All done!"
