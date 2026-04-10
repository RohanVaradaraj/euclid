"""
Based on an input exposure time, what masses/SFR can we see?

"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

#! Plotting configuration
plt.rcParams.update({
    'axes.linewidth': 2.5,
    'font.size': 15,
    'figure.dpi': 100
})

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    # 'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 2, 'ytick.major.width': 2,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 1.5, 'ytick.minor.width': 1.5,
})

#! #################### store results here ###########################
#z_6p52 = [nan, nan, nan, nan, nan, 343.90013224416276, 611.5505242794572, 1087.5077055246672, 1933.8925609931878, 247.57404956964152, 250.0213233415704, 244.96299563641782, 228.28017660072487, 241.8366212521438]
#z_5p38 = [nan, nan, nan, nan, nan, 343.90013224416276, 611.5505242794572, 1087.5077055246672, 186.89180011592623, 247.57404956964152, 190.2739644554019, 189.05880865024926, 178.40836785516177, 241.8366212521438]

# #! SSFR
# z_6p52 = [1.9338925609931878e-07, 1.9338925609931876e-07, 1.9338925609931884e-07, 1.784083678551627e-08, 9.779017801680628e-09, 5.402514320671059e-09, 3.0434426206665054e-09, 1.8650294449763956e-09, 1.3599472607448706e-09, 8.064816928306754e-10, 3.5109115137307125e-10, 3.510911513730712e-10, 3.510911513730712e-10, 3.510911513730712e-10]
# z_5p38 = [np.nan, 1.9338925609931876e-07, 1.9338925609931884e-07, 1.9338925609931884e-07, 1.933892560993188e-07, 1.3922111911773334e-08, 7.90636845368783e-09, 4.356126513616963e-09, 2.7844223823546677e-09, 1.3599472607448712e-09, 8.064816928306756e-10, 4.688570516748192e-10, 3.510911513730712e-10, 3.510911513730712e-10]
# z_5p12 = [np.nan, np.nan, 1.9338925609931884e-07, 1.9338925609931884e-07, 1.933892560993188e-07, 1.8689180011592623e-08, 1.3922111911773336e-08, 6.0169910710898985e-09, 3.3619938670922686e-09, 1.865029444976395e-09, 1.359947260744871e-09, 8.064816928306754e-10, 3.510911513730712e-10, 3.510911513730712e-10]
# z_4p76 = [1.9338925609931878e-07, 1.9338925609931876e-07, 1.7840836785516274e-08, 8.712253027233928e-09, 4.556788418556052e-09, 2.7844223823546677e-09, 1.7840836785516179e-09, 8.064816928306755e-10, 4.688570516748193e-10, 3.510911513730712e-10, 3.5109115137307125e-10, 3.510911513730712e-10, 3.510911513730712e-10, 2.691828092752804e-11]
# z_4p64 = [1.9338925609931878e-07, 1.9338925609931876e-07, 1.3922111911773343e-08, 7.906368453687833e-09, 4.356126513616962e-09, 2.7844223823546677e-09, 1.3599472607448708e-09, 8.064816928306755e-10, 4.688570516748193e-10, 3.510911513730712e-10, 3.5109115137307125e-10, 3.510911513730712e-10, 3.510911513730712e-10, 2.691828092752804e-11]

# M_grid = np.arange(8.0, 11.5, 0.25)

# results = [z_6p52, z_5p38, z_5p12, z_4p76, z_4p64]
# redshifts = [6.52, 5.38, 5.12, 4.76, 4.64]

# # foR EACH result, restrict to values less than 1e-8, then take max as representative value
# sfr_limits = []
# for res in results:
#     res = np.array(res)
#     res = res[res < 1e-8]
#     if len(res) > 0:
#         sfr_limits.append(np.max(res))
#     else:
#         sfr_limits.append(np.nan)
# print(sfr_limits)
# # Plot 


# limiting_fluxes = [5e-18, 2e-17, 3e-17, 5e-18, 5e-18]

# # Plot limiting flux vs redshift and colour by sSFR
# plt.figure(figsize=(10,5))

# # Convert from /yr to /Gyr
# sfr_limits = np.array(sfr_limits) * 1e9

# # Log sSFR
# log_ssfr = np.log10(sfr_limits)

# # Log limiting fluxes
# limiting_fluxes = np.array(limiting_fluxes)

# sc = plt.scatter(redshifts, limiting_fluxes, c=sfr_limits, s=500, cmap='inferno', marker='*', edgecolor='black')
# #sc = plt.scatter(redshifts, log_ssfr, c=np.log10(limiting_fluxes), s=100, cmap='viridis')
# plt.yscale('log')
# plt.xlabel('z')
# plt.ylabel(r'H$\alpha$/$\rm [OIII]$ flux limit ($\rm erg \, s^{-1} \, cm^{-2}$)')
# plt.colorbar(sc, label=r'sSFR ($\rm{Gyr}^{-1})$')
# plt.tight_layout()
# #plt.colorbar(sc, label='H-alpha flux limit (erg/s/cm2)')
# # Vertical line at z=5.25
# plt.axvline(5.25, color='gray', linestyle='--')
# plt.text(5.11, 5.14e-18, r'$\rm H\alpha$', rotation=0, verticalalignment='bottom', fontsize=15)
# plt.text(5.3, 5.14e-18, r'$\rm [OIII]$', rotation=0, verticalalignment='bottom', fontsize=15)
# plt.savefig('grism_flux_limits_ssfr.pdf', dpi=100)
# plt.show()

# exit()

z = 5.38
#line ="halpha"# "oiii" 
line='oiii'

line_dict = {'oiii':'oiii_flux_erg_s_cm2',
                'halpha':'halpha_flux_erg_s_cm2'
}

# Limiting flux
flux_fivesig = 1e-16
flux_fivesig = 2e-17

t = Table.read(f'bagpipes_delayed_tau_grid_z{z}.fits')

# restrict to line fluxes above limit 
t = t[t[line_dict[line]] > flux_fivesig]
print(t)


log_massformed = np.arange(8.0, 11.5, 0.25)

# At each mass, find the minimum SFR that is detected
sfr_limit = []
for lm in log_massformed:
    t_sub = t[np.abs(t['log_massformed'] - lm) < 0.1]
    if len(t_sub) > 0:
        min_sfr = np.min(t_sub['sfr_inst_msun_per_yr'])
        sfr_limit.append(min_sfr)
    else:
        sfr_limit.append(np.nan)

#print(sfr_limit)

# Now compute specific SFR
sfr_limit = np.array(sfr_limit)
ssfr_limit = sfr_limit / 10**log_massformed
print(list(ssfr_limit))

# Convert to per Gyr
ssfr_limit = ssfr_limit * 1e9

# Plot sfr limit in y axis mass, x axis log mass formed
plt.figure()
plt.plot(log_massformed, ssfr_limit, marker='o')
plt.yscale('log')
plt.xlabel('M*/Msun')
plt.ylabel('sSFR')
plt.show()