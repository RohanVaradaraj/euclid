"""
Initial selection of high-z Lyman break candidates.

Created: Friday 12th July 2024.

"""

import numpy as np
from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

cat_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS' / 'other'
cat_name = 'COSMOS_detYJH_masked_1.8as_Euclid_CWEB_2024_07_12.fits'

t = Table.read(cat_dir / cat_name)

print('INITIAL LENGTH: ', len(t))

#! Assert 5sigma detection in Ye
t = t[t['flux_Ye']/t['err_Ye'] > 5]
print('5sigma cut in Ye: ', len(t))

# Convert t['RA'] column from string to float
t['RA'] = t['RA'].astype(float)

#! Assert 3sigma detection in HSC-Y and VISTA-Y
#t = t[(t['flux_HSC-Y_DR3']/t['err_HSC-Y_DR3'] > 2) & (t['flux_Y']/t['err_Y'] > 2)]

# Or just in VISTA-Y
t = t[(t['flux_Y']/t['err_Y'] > 2)]
print('3sigma cut in HSC-Y and VISTA-Y: ', len(t))

#! Assert <2sigma in HSC-G
t = t[t['flux_HSC-G_DR3']/t['err_HSC-G_DR3'] < 2]

print('2sigma non-detection in HSC-G: ', len(t))

#! Assert <2sigma in HSC-R
t = t[t['flux_HSC-R_DR3']/t['err_HSC-R_DR3'] < 2]

print('2sigma non-detection in HSC-R: ', len(t))

#! Assert <2sigma in HSC-I
t = t[t['flux_HSC-I_DR3']/t['err_HSC-I_DR3'] < 2]

print('2sigma non-detection in HSC-I: ', len(t))

#! Assert <2sigma in VIS
#t = t[t['flux_VIS']/t['err_VIS'] < 2]

#print('2sigma non-detection in VIS: ', len(t))

#! Save table
save_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
save_name = 'COSMOS_5sig_Ye_2sig_VISTA_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'
t.write(save_dir / save_name, overwrite=True)



