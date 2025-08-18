"""
Make a latex table of the first few sources.

Created: Friday 27th June 2025.
"""

from astropy.table import Table
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

# === File paths ===
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
euclid_file = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits'

# === Load the catalog ===
t = Table.read(cat_dir / euclid_file)

# === Select columns of interest ===
columns_of_interest = ['ID', 'RA', 'DEC', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'Zphot_sec', 'Chi2_sec',
                       'Stellar_model', 'Chi2_star', 'Lyman_alpha_EW', 'Muv', 'dMuv_inf',
                       'dMuv_sup', 'zmax', 'Vmax']
t = t[columns_of_interest]

# === Calculate Zphot errors ===
zup = t['Zsup'] - t['Zphot']
zdown = t['Zphot'] - t['Zinf']

# === Helper: Truncate to N decimal places ===
def truncate(value, decimals):
    factor = 10 ** decimals
    return int(value * factor) / factor

# === Function to create IAU-compliant Euclid name ===
def make_euclid_name(ra_deg, dec_deg):
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    ra_h = c.ra.hms
    dec_d = c.dec.dms

    # Truncate RA to two decimal places (i.e., 0.01s)
    ra_sec_trunc = truncate(ra_h.s, 2)
    ra_str = f"{int(ra_h.h):02d}{int(ra_h.m):02d}{int(ra_sec_trunc):02d}.{int((ra_sec_trunc % 1) * 100):02d}"

    # Truncate Dec to one decimal place (i.e., 0.1")
    sign = '+' if dec_d.d >= 0 else '-'
    dec_d_abs = abs(dec_d.d)
    dec_m_abs = abs(dec_d.m)
    dec_s_trunc = truncate(abs(dec_d.s), 1)
    dec_str = f"{sign}{int(dec_d_abs):02d}{int(dec_m_abs):02d}{int(dec_s_trunc):02d}.{int((dec_s_trunc % 1) * 10):1d}"

    return f"EUCL\\,J{ra_str}${sign}${dec_str[1:]}"  # replace - with $-$ in LaTeX

print(make_euclid_name(150.11833152095758, 2.2522416552619746))
exit()



# === LaTeX Table Header ===
header = r"""\begin{table}[ht]
\centering
\caption{First five sources from the Euclid candidate catalog.}
\begin{tabular}{cccccccccccccccccc}
\hline
ID & RA & DEC & $z_\mathrm{phot}$ & $\chi^2$ & $z_\mathrm{phot,sec}$ & $\chi^2_\mathrm{sec}$ & Model & $\chi^2_\star$ & EW$_\mathrm{Ly\alpha}$ & $M_\mathrm{UV}$ \\
 & [deg] & [deg] &  &  &  &  &  &  & [\AA] &  &  &  \\
\hline
"""

# === LaTeX Table Rows with IAU IDs ===
rows = ""
for i in range(5):
    row = t[i]
    zphot_fmt = f"${row['Zphot']:.2f}^{{+{zup[i]:.2f}}}_{{-{zdown[i]:.2f}}}$"
    muv_fmt = f"${row['Muv']:.2f}^{{+{row['dMuv_sup']:.2f}}}_{{-{row['dMuv_inf']:.2f}}}$"
    
    euclid_name = make_euclid_name(row['RA'], row['DEC'])
    
    rows += (
        f"{euclid_name} & "
        f"{row['RA']:.5f} & {row['DEC']:.5f} & {zphot_fmt} & "
        f"{row['Chi2']:.2f} & {row['Zphot_sec']:.2f} & {row['Chi2_sec']:.2f} & "
        f"{row['Stellar_model']} & {row['Chi2_star']:.2f} & {row['Lyman_alpha_EW']:.1f} & "
        f"{muv_fmt} \\\\\n"
    )

# === Footer ===
footer = r"""\hline
\end{tabular}
\end{table}
"""

# === Output LaTeX ===
latex_table = header + rows + footer
print(latex_table)