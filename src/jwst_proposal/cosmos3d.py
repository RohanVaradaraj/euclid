"""
Check if any of the sample exists in the COSMOS-3D survey.

Created: Friday 26th September 2025.
"""

import os
import importlib
import subprocess
import sys

# Set environment variables for JWST calibration + NIRCam configuration
os.environ['CRDS_CONTEXT'] = "jwst_1322.pmap"
os.environ['NIRCAM_CONF_VERSION'] = "V9"

def pip_install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

try:
    import grizli
except ImportError:
    # Install required packages
    pip_install("grizli[aws,jwst]")
    pip_install("msaexp")
    pip_install("git+https://github.com/karllark/dust_attenuation.git")

    # Try re-import after install
    importlib.invalidate_caches()
    import grizli