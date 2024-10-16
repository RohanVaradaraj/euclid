import sys

imageName = sys.argv[1]
whtName = sys.argv[2]
whtType = sys.argv[3]
zeropoint = float(sys.argv[4])
outputDir = sys.argv[5]
filterName = sys.argv[6]
overwriteSt = sys.argv[7]
starsOnlySt = sys.argv[8]
fieldName = sys.argv[9]
depthDir = sys.argv[10]
use_cat = sys.argv[11]

if starsOnlySt == 'True':
    starsOnly = True
elif starsOnlySt == 'False':
    starsOnly = False
else:
    print("BAD stuff happening... ")
    exit()

if overwriteSt == 'True':
    overwrite = True
elif overwriteSt == 'False':
    overwrite = False
else:
    print("BAD stuff happening... ")
    exit()
    
from new_psf_codes import psfex

print("Running STUPID PSF")
print(imageName)

# segseg = whtName.split('.')
# segName = segseg[0] + '_seg.fits'
psfex(imageName, filterName, fieldName, zeropoint, depthDir, wht_name = whtName, wht_type = whtType, output_dir = outputDir, overwrite = overwrite, stars_only = starsOnly, use_cat = use_cat)

# Print the command: psfex seems to fail on queue. So run this in python3 on command line!
command = f"psfex('{imageName}', '{filterName}', '{fieldName}', {zeropoint}, '{depthDir}', wht_name='{whtName}', wht_type='{whtType}', output_dir='{outputDir}', overwrite={overwriteSt}, stars_only={starsOnlySt})"
print(command)