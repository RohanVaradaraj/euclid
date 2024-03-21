import sys
import numpy as np

imageName = sys.argv[1]
whtName = sys.argv[2]
whtType = sys.argv[3]
zeropoint = float(sys.argv[4])
outputDir = sys.argv[5]
stripSt = sys.argv[6]
filterName = sys.argv[7]
overwriteSt = sys.argv[8]
mask = sys.argv[9]
gridSepAS = float(sys.argv[10])
apDi = sys.argv[11]

# make the apDi a float array again
diameters = apDi.split(',')
apDiametersAS = np.array(diameters)
apDiametersAS = apDiametersAS.astype(float)

if stripSt == 'True':
    strip = True
elif stripSt == 'False':
    strip = False
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
    
from new_depth_codes import *

print("Running STUPID")
print(imageName)

#segseg = whtName.split('.')
#segName = segseg[0] + '_seg.fits'

image_depth(imageName, zeropoint, wht_name = strwhtName, wht_type = whtType, output_dir = outputDir, strips = strip, filter_name = filterName, overwrite = overwrite, mask = mask, gridSepAS = gridSepAS, ap_diametersAS = apDiametersAS)    
