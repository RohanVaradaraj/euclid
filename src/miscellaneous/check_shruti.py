"""
check weird I-band sources for Shruti!

Created: Friday 21st June 2024.
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.patches import Rectangle, Ellipse

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

linewidth=4

t = Table.read('CDFS_DES_BAGPIPES.csv', format='ascii.csv')

print(t.colnames)

voice_i = t['FLUX_AUTO_VOICE_i_µJy']
hsc_i = t['FLUX_AUTO_HSC_I_µJy']

#! Percentage difference
# Find where voice_i and hsc_i differ by more than 10%
diff = np.abs(voice_i - hsc_i) / voice_i
print(diff)

# Find the indices of the sources that differ by more than some threshold
thresh = 0.2
indices = np.where((diff > thresh) & (t['SN_RA'] < 54.))
print(indices)

#! Absolute difference
# Find where voice_i and hsc_i differ by more than 10%
diff = np.abs(voice_i - hsc_i)

# Find the indices of the sources that differ by more than some threshold
# thresh = 3
# indices = np.where( (diff > thresh) & (diff < 10))
# print(indices)

# # Restrict table
t = t[indices]
print(t)

fig, ax = plt.subplots()

ax.scatter(t['SN_RA'], t['SN_DEC'], c='r', s=10)

# Reverse x axis
plt.gca().invert_xaxis()

plt.title(f'i/I band fluxes differ by > {int(thresh * 100)}%')
#plt.title(f'i/I band fluxes differ by > {thresh} uJy')

plt.xlabel('RA')
plt.ylabel('DEC')

voice1 = SkyCoord('03h33m34.5s -27d34m10.8s', unit=(u.hourangle, u.deg))
voice2 = SkyCoord('03h29m02.7s -27d34m00.7s', unit=(u.hourangle, u.deg))
voice3 = SkyCoord('03h29m01.2s -28d34m15.0s', unit=(u.hourangle, u.deg))
voice4 = SkyCoord('03h33m35.2s -28d34m30.8s', unit=(u.hourangle, u.deg))
voice = SkyCoord('03h32m20s -27d48m30s', unit=(u.hourangle, u.deg))

hsc1 = SkyCoord('03h34m52.1s -27d18m57.5s', unit=(u.hourangle, u.deg))
hsc2 = SkyCoord('03h29m59.0s -27d18m43.7s', unit=(u.hourangle, u.deg))
hsc3 = SkyCoord('03h30m15.1s -28d19m21.0s', unit=(u.hourangle, u.deg))
hsc4 = SkyCoord('03h34m45.6s -28d16m31.0s', unit=(u.hourangle, u.deg))

v1_deg = np.array([float(voice1.ra.deg), float(voice1.dec.deg)])
vlen_deg = 0.5
vlen_hr = 60
v1_deg[0] = v1_deg[0] - vlen_deg/np.cos(np.deg2rad(v1_deg[1]))
v1_deg[1] = v1_deg[1] - vlen_deg

v2_deg = np.array([float(voice2.ra.deg), float(voice2.dec.deg)])
v2_deg[0] = v2_deg[0] - vlen_deg/np.cos(np.deg2rad(v2_deg[1]))
v2_deg[1] = v2_deg[1] - vlen_deg

v3_deg = np.array([float(voice3.ra.deg), float(voice3.dec.deg)])
v3_deg[0] = v3_deg[0] - vlen_deg/np.cos(np.deg2rad(v3_deg[1]))
v3_deg[1] = v3_deg[1] - vlen_deg

v4_deg = np.array([float(voice4.ra.deg), float(voice4.dec.deg)])
v4_deg[0] = v4_deg[0] - vlen_deg/np.cos(np.deg2rad(v4_deg[1]))
v4_deg[1] = v4_deg[1] - vlen_deg

v = np.array([float(voice.ra.deg), float(voice.dec.deg)])
v = [52.83, -28.072]

c1 = Ellipse((hsc1.ra.deg, hsc1.dec.deg), 1.5, 1.5, edgecolor='deepskyblue', facecolor='None', linewidth=linewidth)
c2 = Ellipse((hsc2.ra.deg, hsc2.dec.deg), 1.5, 1.5, edgecolor='deepskyblue', facecolor='None', linewidth=linewidth) 
c3 = Ellipse((hsc3.ra.deg, hsc3.dec.deg), 1.5, 1.5, edgecolor='deepskyblue', facecolor='None', linewidth=linewidth)
c4 = Ellipse((hsc4.ra.deg, hsc4.dec.deg), 1.5, 1.5, edgecolor='deepskyblue', facecolor='None', linewidth=linewidth, label='HSC')

vlen = 1.0
v[0] = v[0] - vlen / np.cos(np.deg2rad(v[1]))
v[1] = v[1] - vlen

rv = Rectangle((v[0], v[1]), 2*vlen/np.cos(np.deg2rad(v[1])), 2*vlen, edgecolor='red', facecolor='None', linewidth=linewidth, label='VOICE', linestyle='--')

ax.add_patch(c1)
ax.add_patch(c2)
ax.add_patch(c3)
ax.add_patch(c4)

ax.add_patch(rv)
plt.legend()


plt.show()
