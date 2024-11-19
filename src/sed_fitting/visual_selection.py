"""
visual_selection.py

Goes through the SED/cutout plot order and prompts the user whether to keep each object.

Created: Sunday 27th October 2024.
"""

from pathlib import Path
import glob
import sys
import json
import shutil

#! Main object type to plot
object_type = 'highz'

if len(sys.argv) > 1:
    filters_json = sys.argv[1]
    filters = json.loads(filters_json)
    bools_json = sys.argv[2]
    bools = json.loads(bools_json)
    all_filters_json = sys.argv[3]
    all_filters = json.loads(all_filters_json)

#! Output PDF name setup from detection filters
det_list = [f for f, t in filters.items() if t['type'] == 'detection']
if len(det_list) != 0:
    base_det = 'det_' + '_'.join(det_list)
else:
    stack_list = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    stack_filters = stack_list[0].split('+')
    base_det = 'det_' + '_'.join(stack_filters)
    det_list = stack_filters

# Directory setup
zphot_folder = base_det + f'_best_{object_type}'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / zphot_folder
print('Running visual selection on:', zphot_folder)

# Output directories
good_dir = zphot_dir.parents[0] / (zphot_folder + '_good')
bad_dir = zphot_dir.parents[0] / (zphot_folder + '_bad')
maybe_dir = zphot_dir.parents[0] / (zphot_folder + '_maybe')

# Create directories if they don't exist
good_dir.mkdir(parents=True, exist_ok=True)
bad_dir.mkdir(parents=True, exist_ok=True)
maybe_dir.mkdir(parents=True, exist_ok=True)

# Progress tracking file
progress_file = Path("progress.json")

# Load previous progress if available
if progress_file.exists():
    with open(progress_file, "r") as f:
        progress = json.load(f)
else:
    progress = {"last_index": 0}

# Collect and sort .spec files
spec_files = glob.glob(str(zphot_dir / '*.spec'))
spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

# Check if the user wants to start from a custom ID
custom_start = input("Do you want to start from a specific ID? (y/n): ").strip().lower()
if custom_start == 'y':
    custom_id = input("Enter the ID to start from: ").strip()
    
    # Find the index of the file with the matching ID
    custom_start_index = next(
        (i for i, file in enumerate(spec_files) if custom_id in file),
        None
    )
    
    if custom_start_index is not None:
        # Update progress with the custom starting index
        progress["last_index"] = custom_start_index
        print(f"Starting from ID {custom_id}.")
    else:
        print(f"ID {custom_id} not found. Starting from saved progress or the beginning.")


# Resume from the last processed file
for i, spec_file in enumerate(spec_files[progress["last_index"]:], start=progress["last_index"]):
    print(f'Object {i + 1} of {len(spec_files)}')
    ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]
    print('ID:', ID)

    # Prompt for input
    user_input = input("Press 'Q' for good, 'W' for maybe, Enter for bad, or type 'STOP' to save and exit: ").strip().lower()

    # Handle user input
    if user_input == 'q':
        shutil.copy2(spec_file, good_dir)
        print(f"File copied to {good_dir}")
    elif user_input == 'w':
        shutil.copy2(spec_file, maybe_dir)
        print(f"File copied to {maybe_dir}")
    elif user_input == '':
        shutil.copy2(spec_file, bad_dir)
        print(f"File copied to {bad_dir}")
    elif user_input == 'stop':
        progress["last_index"] = i
        with open(progress_file, "w") as f:
            json.dump(progress, f)
        print("Progress saved. Exiting.")
        sys.exit(0)
    else:
        print("Invalid input. Please try again.")

# Clear progress after all files are processed
if progress["last_index"] >= len(spec_files) - 1:
    progress_file.unlink(missing_ok=True)
    print("All files processed. Progress file cleared.")
