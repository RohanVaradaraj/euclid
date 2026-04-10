"""
Format the sample fits files into the required format for uploading to CDS with the accepted paper.

Created: Wednesday 14th January 2026. Very sleepy today.
"""

from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def truncate(value, decimals):
    factor = 10 ** decimals
    return int(value * factor) / factor

def make_euclid_name(ra_deg, dec_deg, uvista=False):
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    ra_h = c.ra.hms
    dec_d = c.dec.dms

    # Truncate RA to two d.p.
    ra_sec_trunc = truncate(ra_h.s, 2)
    ra_str = f"{int(ra_h.h):02d}{int(ra_h.m):02d}{int(ra_sec_trunc):02d}.{int((ra_sec_trunc % 1) * 100):02d}"

    # Truncate DEC to one d.p.
    sign = '+' if dec_d.d >= 0 else '-'
    dec_d_abs = abs(dec_d.d)
    dec_m_abs = abs(dec_d.m)
    dec_s_trunc = truncate(abs(dec_d.s), 1)
    dec_str = f"{sign}{int(dec_d_abs):02d}{int(dec_m_abs):02d}{int(dec_s_trunc):02d}.{int((dec_s_trunc % 1) * 10):1d}"

    if uvista:
        return f"UVISTA J{ra_str}{sign}{dec_str[1:]}"
    else:
        return f"EUCL J{ra_str}{sign}{dec_str[1:]}"  

U_plus_E = True
U_only = False

cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'

if U_plus_E:
    ue_cat = Table.read(cat_dir / 'Euclid_UltraVISTA_z7_sample_kron_piecewise.fits')

    # Create new columns
    euclid_names = []
    for row in ue_cat:
        ra = row['RA']
        dec = row['DEC']
        euclid_name = make_euclid_name(ra, dec)
        euclid_names.append(euclid_name)
    ue_cat['SOURCE_ID'] = euclid_names

    # Change all flux_ column units to erg/s/cm2/Hz
    for col in ue_cat.colnames:
        if col.startswith('flux_'):
            ue_cat[col].unit = u.erg / u.s / u.cm**2 / u.Hz

    # Reorder columns
    new_order = ['SOURCE_ID', 'RA', 'DEC'] + [col for col in ue_cat.colnames if col not in ['SOURCE_ID', 'RA', 'DEC']]
    ue_cat = ue_cat[new_order]

    # Remove X_IMAGE, Y_IMAGE
    if 'X_IMAGE' in ue_cat.colnames:
        ue_cat.remove_column('X_IMAGE')
    if 'Y_IMAGE' in ue_cat.colnames:
        ue_cat.remove_column('Y_IMAGE')
    if 'ID' in ue_cat.colnames:
        ue_cat.remove_column('ID')  
    if 'zmax' in ue_cat.colnames:  
        ue_cat.remove_column('zmax')
    if 'Vmax' in ue_cat.colnames:
        ue_cat.remove_column('Vmax')
    if 'CLASS' in ue_cat.colnames:
        ue_cat.remove_column('CLASS')
    if 'separation' in ue_cat.colnames:
        ue_cat.remove_column('separation')
    if 'ra_irac' in ue_cat.colnames:
        ue_cat.remove_column('ra_irac')
    if 'dec_irac' in ue_cat.colnames:
        ue_cat.remove_column('dec_irac')


    # Remove fluxes which have CFHT in column name
    for col in ue_cat.colnames:
        if 'CFHT' in col:
            ue_cat.remove_column(col)

    # If column name contains HSC, remove _DR3 and replace - with _
    for col in ue_cat.colnames:
        if 'HSC' in col:
            new_col_name = col.replace('_DR3', '').replace('-', '_')
            ue_cat.rename_column(col, new_col_name)

    # Rename VIS, Ye, Je, He to IE, YE, JE, HE, if they exist in a column name
    for col in ue_cat.colnames:
        if 'VIS' in col:
            new_col_name = col.replace('VIS', 'IE')
            ue_cat.rename_column(col, new_col_name)
        if 'Ye' in col:
            new_col_name = col.replace('Ye', 'YE')
            ue_cat.rename_column(col, new_col_name)
        if 'Je' in col:
            new_col_name = col.replace('Je', 'JE')
            ue_cat.rename_column(col, new_col_name)
        if 'He' in col:
            new_col_name = col.replace('He', 'HE')
            ue_cat.rename_column(col, new_col_name)

    # If ch1cds and ch2cds are in column names, rename to IRAC1 and IRAC2
    for col in ue_cat.colnames:
        if 'ch1cds' in col:
            new_col_name = col.replace('ch1cds', 'IRAC1')
            ue_cat.rename_column(col, new_col_name)
        if 'ch2cds' in col:
            new_col_name = col.replace('ch2cds', 'IRAC2')
            ue_cat.rename_column(col, new_col_name)

    # Remove filters f115w, f150w, f277w, f444w if they exist
    for col in ue_cat.colnames:
        if 'f115w' in col or 'f150w' in col or 'f277w' in col or 'f444w' in col:
            ue_cat.remove_column(col)


    # Order: SOURCE_ID, RA, DEC, Zphot, Zsup, Zinf, Chi2, Zphot_sec, Chi2_sec, Stellar_model, Chi2_star, Lyman_alpha_EW, Muv, dMuv_sup, dMuv_inf, then everything else.
    final_order = ['SOURCE_ID', 'RA', 'DEC', 'Zphot', 'dZphot_sup', 'dZphot_inf', 'Chi2', 'Zphot_sec', 'Chi2_sec', 'Stellar_model', 'Chi2_star', 'Lyman_alpha_EW', 'Muv', 'dMuv_sup', 'dMuv_inf',
    'flux_HSC_G', 'flux_HSC_R', 'flux_HSC_I', 'flux_HSC_NB0816', 'flux_HSC_Z', 'flux_HSC_NB0921', 'flux_HSC_Y', 'flux_IE', 'flux_YE', 'flux_JE', 'flux_HE', 'flux_Y', 'flux_J', 'flux_H', 'flux_Ks', 'flux_IRAC1', 'flux_IRAC2',
    'err_HSC_G', 'err_HSC_R', 'err_HSC_I', 'err_HSC_NB0816', 'err_HSC_Z', 'err_HSC_NB0921', 'err_HSC_Y', 'err_IE', 'err_YE', 'err_JE', 'err_HE', 'err_Y', 'err_J', 'err_H', 'err_Ks', 'err_IRAC1', 'err_IRAC2']

    # Convert Zsup into Zerr_up and Zinf into Zerr_lo by subtracting Zphot 
    ue_cat['dZphot_sup'] = ue_cat['Zsup'] - ue_cat['Zphot']
    ue_cat['dZphot_inf'] = ue_cat['Zphot'] - ue_cat['Zinf']

    # Place these after Zphot and remove Zsup and Zinf
    # final_order.insert(4, 'dZphot_sup')
    # final_order.insert(5, 'dZphot_inf')

    # Remove Zsup and Zinf
    ue_cat.remove_column('Zsup')
    ue_cat.remove_column('Zinf')

    # Also set err column units to erg/s/cm2/Hz
    for col in final_order:
        if col.startswith('err_'):
            ue_cat[col].unit = u.erg / u.s / u.cm**2 / u.Hz

    ue_cat = ue_cat[final_order]

    # Sort rows by increasing DEC
    ue_cat.sort('DEC')

    # Set any flux/err IRAC1/2 which is < -80 to -99.0
    for row in ue_cat:
        if row['flux_IRAC1'] < -80:
            row['flux_IRAC1'] = -99.0
        if row['flux_IRAC2'] < -80:
            row['flux_IRAC2'] = -99.0
        if row['err_IRAC1'] < -80:
            row['err_IRAC1'] = -99.0
        if row['err_IRAC2'] < -80:
            row['err_IRAC2'] = -99.0

    # Round dZphot_sup and dZphot_inf to 5 decimal places
    ue_cat['dZphot_sup'] = np.round(ue_cat['dZphot_sup'], 5)
    ue_cat['dZphot_inf'] = np.round(ue_cat['dZphot_inf'], 5)

    print(ue_cat.colnames)
    print(ue_cat)

    # Remove Muv faintest source - was removed in internal review
    max_muv_index = np.argmax(ue_cat['Muv'])
    ue_cat.remove_row(max_muv_index)

    # Write out the final table
    name = 'UltraVISTA_plus_Euclid_z7_sample.fits'

    # Save the table here
    ue_cat.write(Path.cwd() / name, overwrite=True)

    # Save in .ascii format too
    ue_cat.write(Path.cwd() / name.replace('.fits', '.ascii'), format='ascii.fixed_width', overwrite=True)


if U_only:
    uo_cat = Table.read(cat_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_kron_piecewise.fits')

    uo_cat = uo_cat[uo_cat['Zphot'] < 7.5]

    # Sort by Muv
    uo_cat.sort('Muv')

    # Remove sources discounted after paper corrections.
    uo_cat.remove_row(-1)
    uo_cat.remove_row(-1)

    names = []
    for row in uo_cat:
        ra = row['RA']
        dec = row['DEC']
        euclid_name = make_euclid_name(ra, dec, uvista=True)
        names.append(euclid_name)
    
    uo_cat['SOURCE_ID'] = names

    # Change all flux_ column units to erg/s/cm2/Hz
    for col in uo_cat.colnames:
        if col.startswith('flux_'):
            uo_cat[col].unit = u.erg / u.s / u.cm**2 / u.Hz

    # Reorder columns
    new_order = ['SOURCE_ID', 'RA', 'DEC'] + [col for col in uo_cat.colnames if col not in ['SOURCE_ID', 'RA', 'DEC']]
    uo_cat = uo_cat[new_order]

    # Remove X_IMAGE, Y_IMAGE
    if 'X_IMAGE' in uo_cat.colnames:
        uo_cat.remove_column('X_IMAGE')
    if 'Y_IMAGE' in uo_cat.colnames:
        uo_cat.remove_column('Y_IMAGE')
    if 'ID' in uo_cat.colnames:
        uo_cat.remove_column('ID')  
    if 'zmax' in uo_cat.colnames:  
        uo_cat.remove_column('zmax')
    if 'Vmax' in uo_cat.colnames:
        uo_cat.remove_column('Vmax')
    if 'CLASS' in uo_cat.colnames:
        uo_cat.remove_column('CLASS')
    if 'separation' in uo_cat.colnames:
        uo_cat.remove_column('separation')
    if 'ra_irac' in uo_cat.colnames:
        uo_cat.remove_column('ra_irac')
    if 'dec_irac' in uo_cat.colnames:
        uo_cat.remove_column('dec_irac')


    # Remove fluxes which have CFHT in column name
    for col in uo_cat.colnames:
        if 'CFHT' in col:
            uo_cat.remove_column(col)

    # If column name contains HSC, remove _DR3 and replace - with _
    for col in uo_cat.colnames:
        if 'HSC' in col:
            new_col_name = col.replace('_DR3', '').replace('-', '_')
            uo_cat.rename_column(col, new_col_name)

    # Rename VIS, Ye, Je, He to IE, YE, JE, HE, if they exist in a column name
    for col in uo_cat.colnames:
        if 'VIS' in col or 'Ye' in col or 'Je' in col or 'He' in col:
            uo_cat.remove_column(col)

    # If ch1cds and ch2cds are in column names, rename to IRAC1 and IRAC2
    for col in uo_cat.colnames:
        if 'ch1cds' in col:
            new_col_name = col.replace('ch1cds', 'IRAC1')
            uo_cat.rename_column(col, new_col_name)
        if 'ch2cds' in col:
            new_col_name = col.replace('ch2cds', 'IRAC2')
            uo_cat.rename_column(col, new_col_name)

    # Remove filters f115w, f150w, f277w, f444w if they exist
    for col in uo_cat.colnames:
        if 'f115w' in col or 'f150w' in col or 'f277w' in col or 'f444w' in col:
            uo_cat.remove_column(col)
    
    # Order: SOURCE_ID, RA, DEC, Zphot, Zsup, Zinf, Chi2, Zphot_sec, Chi2_sec, Stellar_model, Chi2_star, Lyman_alpha_EW, Muv, dMuv_sup, dMuv_inf, then everything else.
    final_order = ['SOURCE_ID', 'RA', 'DEC', 'Zphot', 'dZphot_sup', 'dZphot_inf', 'Chi2', 'Zphot_sec', 'Chi2_sec', 'Stellar_model', 'Chi2_star', 'Lyman_alpha_EW', 'Muv', 'dMuv_sup', 'dMuv_inf',
    'flux_HSC_G', 'flux_HSC_R', 'flux_HSC_I', 'flux_HSC_NB0816', 'flux_HSC_Z', 'flux_HSC_NB0921', 'flux_HSC_Y', 'flux_Y', 'flux_J', 'flux_H', 'flux_Ks', 'flux_IRAC1', 'flux_IRAC2',
    'err_HSC_G', 'err_HSC_R', 'err_HSC_I', 'err_HSC_NB0816', 'err_HSC_Z', 'err_HSC_NB0921', 'err_HSC_Y', 'err_Y', 'err_J', 'err_H', 'err_Ks', 'err_IRAC1', 'err_IRAC2']

    # Convert Zsup into Zerr_up and Zinf into Zerr_lo by subtracting Zphot 
    uo_cat['dZphot_sup'] = uo_cat['Zsup'] - uo_cat['Zphot']
    uo_cat['dZphot_inf'] = uo_cat['Zphot'] - uo_cat['Zinf']

    # Place these after Zphot and remove Zsup and Zinf
    # final_order.insert(4, 'dZphot_sup')
    # final_order.insert(5, 'dZphot_inf')

    # Remove Zsup and Zinf
    uo_cat.remove_column('Zsup')
    uo_cat.remove_column('Zinf')

    # Also set err column units to erg/s/cm2/Hz
    for col in final_order:
        if col.startswith('err_'):
            uo_cat[col].unit = u.erg / u.s / u.cm**2 / u.Hz

    uo_cat = uo_cat[final_order]

    # Sort rows by increasing DEC
    uo_cat.sort('DEC')

    # Set any flux/err IRAC1/2 which is < -80 to -99.0
    for row in uo_cat:
        if row['flux_IRAC1'] < -80:
            row['flux_IRAC1'] = -99.0
        if row['flux_IRAC2'] < -80:
            row['flux_IRAC2'] = -99.0
        if row['err_IRAC1'] < -80:
            row['err_IRAC1'] = -99.0
        if row['err_IRAC2'] < -80:
            row['err_IRAC2'] = -99.0

    # Write out the final table
    name = 'UltraVISTA_only_z7_sample.fits'

    print(uo_cat)
    print(uo_cat.colnames)

    print(uo_cat[(uo_cat['err_IRAC1'] < 1e-33) & (uo_cat['err_IRAC1'] > 0)])
    # Get indices where err_IRAC1 is < 1e-33 but > 0
    indices = np.where((uo_cat['err_IRAC1'] < 1e-33) & (uo_cat['err_IRAC1'] > 0))[0]
    print(indices)


    for row in uo_cat:
        if (row['err_IRAC1'] < 1e-33) & (row['err_IRAC1'] > 0):

            # Set minimum as np.abs * 20% of flux
            row['err_IRAC1'] = np.abs(row['flux_IRAC1']) * 0.2

    # Round dZphot_sup and dZphot_inf to 5 decimal places
    uo_cat['dZphot_sup'] = np.round(uo_cat['dZphot_sup'], 5)
    uo_cat['dZphot_inf'] = np.round(uo_cat['dZphot_inf'], 5)

    # Print same indices again to check
    print(uo_cat[indices])

    # Save the table here
    uo_cat.write(Path.cwd() / name, overwrite=True)

    # Save in .ascii format too
    uo_cat.write(Path.cwd() / name.replace('.fits', '.ascii'), format='ascii.fixed_width', overwrite=True)
