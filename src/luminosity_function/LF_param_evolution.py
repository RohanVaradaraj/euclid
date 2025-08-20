"""
Plot the evolution of LF parameters with redshift

Created: Thursday 6th February 2025
"""

import numpy as np
from astropy.table import Table, vstack, Column
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit


plt.rcParams.update({'font.size': 25})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 10, 'ytick.major.size': 10,
    'xtick.major.width': 3, 'ytick.major.width': 3,

    # Minor tick size and width
    'xtick.minor.size': 5, 'ytick.minor.size': 5,
    'xtick.minor.width': 2, 'ytick.minor.width': 2,
})

ms = 12
elinewidth = 2.5
zorder=2


plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'


# My LF parameters

# my_z = np.array([7])
# my_phi = np.array([9.1146e-05])
# my_phi_err_up = np.array([6.654507642353416e-05])
# my_phi_err_lo = np.array([3.8252481522955595e-05])

# my_M_star = np.array([-21.1304])
# my_M_star_err_up = np.array([0.27047277962511274])
# my_M_star_err_lo = np.array([0.24948880288177477])

# my_alpha = np.array([-2.1085])
# my_alpha_err_up = np.array([0.21369413904090906])
# my_alpha_err_lo = np.array([0.17156830499011377])

# my_beta = np.array([-4.5987])
# my_beta_err_up = np.array([0.32097933349011587])
# my_beta_err_lo = np.array([0.3703468864815038])

#? My LF parameters after kron and larger errors
my_z = np.array([7])
my_phi = np.array([9.0598e-05])
my_phi_err_up = np.array([6.704259904493972e-05])
my_phi_err_lo = np.array([3.806543369065727e-05])

my_M_star = np.array([-21.1431])
my_M_star_err_up = np.array([0.27714872059357987])
my_M_star_err_lo = np.array([0.25054261288004653])

my_alpha = np.array([-2.1029])
my_alpha_err_up = np.array([0.21467644551429])
my_alpha_err_lo = np.array([0.17300217202221013])

my_beta = np.array([-4.6266])
my_beta_err_up = np.array([0.3357791319020107])
my_beta_err_lo = np.array([0.385521158264428])


# Donnan+23
donnan23_z = np.array([8.0])
donnan23_phi = np.array([3.30e-4])
donnan23_phi_err = np.array([3.41e-4])
donnan23_M_star = np.array([-20.02])
donnan23_M_star_err = np.array([0.55])
donnan23_alpha = np.array([-2.04])
donnan23_alpha_err = np.array([0.29])
donnan23_beta = np.array([-4.26])
donnan23_beta_err = np.array([0.50])

# Donnan+24
donnan24_z = np.array([9.0, 10.0, 11.0, 12.5])
donnan24_phi = np.array([23.5e-5, 14.5e-5, 3.27e-5, 0.99e-5])
donnan24_phi_err = np.array([39.2e-5, 16.4e-5, 9.86e-5, 0.99e-5])
donnan24_M_star = np.array([-19.70, -19.98, -20.73, -20.82])
donnan24_M_star_err = np.array([0.96, 0.61, 1.61, 0.71])

donnan24_alpha = np.array([-2.00, -1.98, -2.19, -99])
donnan24_alpha_err = np.array([0.47, 0.40, 0.69, 0.0])

donnan24_beta = np.array([-3.81, -99, -4.29, -99])
donnan24_beta_err = np.array([0.49, 0., 1.30, 0.0])

# Bowler+20
bowler20_z = np.array([8.05, 9.05])
bowler20_phi = np.array([4.83e-4, 2.85e-4])
bowler20_phi_err = np.array([2.25e-4, 1.39e-4])

bowler20_M_star = np.array([-19.8, -19.67])
bowler20_M_star_err = np.array([0.26, 0.33])

bowler20_alpha = np.array([-1.96, -99])
bowler20_alpha_err = np.array([0.15, 0.])

bowler20_beta = np.array([-3.98, -3.75])
bowler20_beta_err = np.array([0.14, 0.22])

# Bowler+15
bowler15_z = np.array([6.0])
bowler15_phi = np.array([1.9e-4])
bowler15_phi_err_up = np.array([1.2e-4])
bowler15_phi_err_lo = np.array([0.8e-4])

bowler15_M_star = np.array([-21.20])
bowler15_M_star_err = np.array([0.22])

bowler15_alpha = np.array([-2.10])
bowler15_alpha_err_up = np.array([0.16])
bowler15_alpha_err_lo = np.array([0.14])

bowler15_beta = np.array([-5.1])
bowler15_beta_err_up = np.array([0.5])
bowler15_beta_err_lo = np.array([0.6])

# Adams+23
adams23_z = np.array([3.1, 4.0, 4.8])
adams23_phi = np.array([-2.63, -3.00, -3.88])
adams23_phi_err = np.array([-0.03, -0.04, -0.05])

# Compute actual errors from these \pm values
adams23_phi_err = np.array([adams23_phi - 10**(np.log10(adams23_phi) - adams23_phi_err), 10**(np.log10(adams23_phi) + adams23_phi_err) - adams23_phi])

adams23_M_star = np.array([-20.59, -21.11, -20.86])
adams23_M_star_err = np.array([0.03, 0.04, 0.05])

adams23_alpha = np.array([-1.52, -1.80, -1.57])
adams23_alpha_err = np.array([0.03, 0.03, 0.06])

adams23_beta = np.array([-4.95, -5.29, -5.73])
adams23_beta_err_up = np.array([0.08, 0.14, 0.18])
adams23_beta_err_lo = np.array([0.09, 0.15, 0.19])

# Finkelstein+15
finkelstein15_z = np.array([4.0, 5.0, 6.05, 7.05, 8.1])
finkelstein15_phi = np.array([14.1e-4, 8.95e-4, 1.86e-4, 1.57e-4, 0.72e-4])
finkelstein15_phi_err_up = np.array([2.05e-4, 1.92e-4, 0.94e-4, 1.49e-4, 2.52e-4])
finkelstein15_phi_err_lo = np.array([1.85e-4, 1.31e-4, 0.8e-4, 0.95e-4, 0.65e-4])

finkelstein15_M_star = np.array([-20.73, -20.81, -21.13, -21.03, -20.89])
finkelstein15_M_star_err_up = np.array([0.09, 0.13, 0.25, 0.37, 0.74])
finkelstein15_M_star_err_lo = np.array([0.09, 0.12, 0.31, 0.50, 1.08])     

finkelstein15_alpha = np.array([-1.56, -1.67, -2.02, -2.03, -2.36])
finkelstein15_alpha_err_up = np.array([0.06, 0.05, 0.10, 0.21, 0.54])
finkelstein15_alpha_err_lo = np.array([0.05, 0.06, 0.10, 0.20, 0.40])

# Bouwens+21
bouwens21_z = np.array([2.1, 2.9, 3.8, 4.9, 5.9, 6.8, 7.9, 8.9])
bouwens21_phi = np.array([4e-3, 2.1e-3, 1.69e-3, 0.79e-3, 0.51e-3, 0.19e-3, 0.09e-3, 0.021e-3])
bouwens21_phi_err_up = np.array([0.5, 0.3, 0.22, 0.16, 0.12, 0.08, 0.09, 0.014]) * 1e-3
bouwens21_phi_err_lo = np.array([0.4, 0.3, 0.2, 0.13, 0.1, 0.06, 0.05, 0.009]) * 1e-3

bouwens21_M_star = np.array([-20.28, -20.87, -20.93, -21.1, -20.93, -21.15, -20.93, -99])
bouwens21_M_star_err = np.array([0.09, 0.09, 0.08, 0.11, 0.09, 0.13, 0.28, 0.0])

bouwens21_alpha = np.array([-1.52, -1.61, -1.69, -1.74, -1.93, -2.06, -2.23, -2.33])
bouwens21_alpha_err = np.array([0.03, 0.03, 0.03, 0.06, 0.08, 0.11, 0.20, 0.19])

# Bowler+17
bowler17_z = np.array([7.1])
bowler17_phi = np.array([2.3e-4])
bowler17_phi_err_up = np.array([1.8e-4])
bowler17_phi_err_lo = np.array([0.9e-4])

bowler17_M_star = np.array([-20.6])
bowler17_M_star_err_up = np.array([0.33])
bowler17_M_star_err_lo = np.array([0.27])

bowler17_alpha = np.array([-2.19])
bowler17_alpha_err_up = np.array([0.12])
bowler17_alpha_err_lo = np.array([0.1])

bowler17_beta = np.array([-4.6])
bowler17_beta_err_up = np.array([0.4])
bowler17_beta_err_lo = np.array([0.5])

# Harikane+24
harikane24_z = np.array([7.15, 8.15, 10.05, 12])
harikane24_phi = np.array([-3.74, -4.1, -4.5, -4.82])
harikane24_phi_err_up = np.array([0.23, 0.22, 0.66, 0.52])
harikane24_phi_err_lo = np.array([0.22, 0.65, 0.68, 0.40])

# Compute actual errors from these \pm values
# harikane24_phi_err = np.array([harikane24_phi - 10**(np.log10(harikane24_phi) - harikane24_phi_err), 10**(np.log10(harikane24_phi) + harikane24_phi_err) - harikane24_phi])
# # First value is lower error, second is upper error

harikane24_M_star = np.array([-21.01, -20.88, -20.61, -99])
harikane24_M_star_err_up = np.array([0.3, 0.25, 0.71, 0.0])
harikane24_M_star_err_lo = np.array([0.26, 0.23, 0.90, 0.0])

harikane24_alpha = np.array([-2.08, -2.27, -99, -99])
harikane24_alpha_err_up = np.array([0.12, 0.16, 0.0, 0.0])
harikane24_alpha_err_lo = np.array([0.11, 0.25, 0.0, 0.0])

harikane24_beta = np.array([-4.81, -4.45, -99, -99])
harikane24_beta_err_up = np.array([0.50, 0.31, 0.0, 0.0])
harikane24_beta_err_lo = np.array([0.56, 2.04, 0.0, 0.0])

# McLeod+24
mcleod24_z = np.array([11.05])

mcleod24_phi = np.array([0.43e-5])
mcleod24_phi_err = np.array([1.15e-5])

mcleod24_M_star = np.array([-21.62])
mcleod24_M_star_err = np.array([1.28])

mcleod24_alpha = np.array([-2.59])
mcleod24_alpha_err = np.array([0.35])

mcleod24_beta = np.array([-4.99])
mcleod24_beta_err = np.array([2.63])

# Chemerynska+24
chemerynska24_z = np.array([10.5])

chemerynska24_phi = np.array([-4.22])
chemerynska24_phi_err = np.array([0.71])
# Convert to actual errors
#chemerynska24_phi_err = np.array([chemerynska24_phi - 10**(np.log10(chemerynska24_phi) - chemerynska24_phi_err), 10**(np.log10(chemerynska24_phi) + chemerynska24_phi_err) - chemerynska24_phi])

chemerynska24_M_star = np.array([-99])
chemerynska24_M_star_err = np.array([0.])

chemerynska24_alpha = np.array([-99])
chemerynska24_alpha_err = np.array([0.])

chemerynska24_beta = np.array([-2.66])
chemerynska24_beta_err = np.array([1.04])

# Whitler+25
whitler25_z = np.array([9.8])

whitler25_phi = np.array([4.849e-5])
whitler25_phi_err_up = np.array([5.324e-5])
whitler25_phi_err_lo = np.array([2.147e-5])

whitler25_M_star = np.array([-20.54])
whitler25_M_star_err_up = np.array([0.39])
whitler25_M_star_err_lo = np.array([0.30])

whitler25_alpha = np.array([-2.60])
whitler25_alpha_err_up = np.array([0.17])
whitler25_alpha_err_lo = np.array([0.19])

whitler25_beta = np.array([-3.49])
whitler25_beta_err_up = np.array([0.75])
whitler25_beta_err_lo = np.array([1.38])

# Adams+24
adams24_z = np.array([7.55, 9.1, 10.55])
adams24_phi = np.array([-3.88, -4.14, -4.93])
adams24_phi_err_up = np.array([0.54, 0.44, 0.62])
adams24_phi_err_lo = np.array([0.42, 0.87, 0.82])

adams24_M_star = np.array([-20.45, -20.60, -21.10])
adams24_M_star_err_up = np.array([0.70, 0.43, 0.78])
adams24_M_star_err_lo = np.array([0.43, 0.24, 0.64])

adams24_alpha = np.array([-2.20, -99, -99])
adams24_alpha_err_up = np.array([0.32, 0, 0])
adams24_alpha_err_lo = np.array([0.23, 0, 0])

adams24_beta = np.array([-4.56, -5.35, -4.45])
adams24_beta_err_up = np.array([0.66, 1.00, 0.97])
adams24_beta_err_lo = np.array([0.62, 1.08, 1.02])

# Perez-gonzalez 23
perez23_z = np.array([9, 10.75, 12.25])
perez23_M_star = np.array([-21, -20.74, -20.81])
perez23_M_star_err_up = np.array([0.34, 0.57, 0.77])
perez23_M_star_err_lo = np.array([0.45, 0.54, 0.67])
perez23_phi = np.array([1.8e-5, 3.9e-5, 1.6e-5])
perez23_phi_err_up = np.array([1.1e-5, 4.0e-5, 1.6e-5])
perez23_phi_err_lo = np.array([1.7e-5, 3.1e-5, 0.9e-5])
perez23_alpha = np.array([-2.3, -2.14, -2.19])
perez23_alpha_err_up = np.array([0.16, 0.24, 0.26])
perez23_alpha_err_lo = np.array([0.24, 0.38, 0.39])


# Convert all phi errors into log space
my_phi_err_up = my_phi_err_up / (my_phi * np.log(10))
my_phi_err_lo = my_phi_err_lo / (my_phi * np.log(10))
donnan23_phi_err = donnan23_phi_err / (donnan23_phi * np.log(10))
donnan24_phi_err = donnan24_phi_err / (donnan24_phi * np.log(10))
bowler20_phi_err = bowler20_phi_err / (bowler20_phi * np.log(10))
bowler15_phi_err_up = bowler15_phi_err_up / (bowler15_phi * np.log(10))
bowler15_phi_err_lo = bowler15_phi_err_lo / (bowler15_phi * np.log(10))
bowler17_phi_err_lo = bowler17_phi_err_lo / (bowler17_phi * np.log(10))
bowler17_phi_err_up = bowler17_phi_err_up / (bowler17_phi * np.log(10))
adams23_phi_err = adams23_phi_err / (adams23_phi * np.log(10))
finkelstein15_phi_err_up = finkelstein15_phi_err_up / (finkelstein15_phi * np.log(10))
finkelstein15_phi_err_lo = finkelstein15_phi_err_lo / (finkelstein15_phi * np.log(10))
bouwens21_phi_err_up = bouwens21_phi_err_up / (bouwens21_phi * np.log(10))
bouwens21_phi_err_lo = bouwens21_phi_err_lo / (bouwens21_phi * np.log(10))
#harikane24_phi_err = harikane24_phi_err / (harikane24_phi * np.log(10))
mcleod24_phi_err = mcleod24_phi_err / (mcleod24_phi * np.log(10))
#chemerynska24_phi_err = chemerynska24_phi_err / (chemerynska24_phi * np.log(10))
whitler25_phi_err_up = whitler25_phi_err_up / (whitler25_phi * np.log(10))
whitler25_phi_err_lo = whitler25_phi_err_lo / (whitler25_phi * np.log(10))
perez23_phi_err_up = perez23_phi_err_up / (perez23_phi * np.log(10))
perez23_phi_err_lo = perez23_phi_err_lo / (perez23_phi * np.log(10))


# take log of all phi values
my_phi = np.log10(my_phi)
donnan23_phi = np.log10(donnan23_phi)
donnan24_phi = np.log10(donnan24_phi)
bowler20_phi = np.log10(bowler20_phi)
bowler15_phi = np.log10(bowler15_phi)
# adams23_phi = np.log10(adams23_phi)
finkelstein15_phi = np.log10(finkelstein15_phi)
bouwens21_phi = np.log10(bouwens21_phi)
bowler17_phi = np.log10(bowler17_phi)
# harikane24_phi = np.log10(harikane24_phi)
mcleod24_phi = np.log10(mcleod24_phi)
# chemerynska24_phi = np.log10(chemerynska24_phi)
whitler25_phi = np.log10(whitler25_phi)
perez23_phi = np.log10(perez23_phi)


# Create figure and axes
fig, axs = plt.subplots(2, 2, figsize=(20, 12))

# Define a lambda function to make plotting easier
plot_data = lambda ax, z, val, err, label, style, zorder: ax.errorbar(
    z, val, yerr=err, fmt=style["marker"], label=label,
    markersize=style["ms"], markeredgecolor='black', 
    elinewidth=style["elinewidth"], color=style["color"],
    zorder=zorder
)

# Define styles for different datasets
styles = {
    "my":          {"color": "tab:red",      "marker": "o",  "ms": 18, "elinewidth": 4},
    "donnan23":    {"color": "blue",     "marker": "s",  "ms": 10,  "elinewidth": 3},
    "donnan24":    {"color": "green",    "marker": "^",  "ms": 10,  "elinewidth": 3},
    "bowler20":    {"color": "orange",   "marker": "d",  "ms": 10,  "elinewidth": 3},
    "bowler15":    {"color": "purple",   "marker": "v",  "ms": 10,  "elinewidth": 3},
    "adams23":     {"color": "gray",    "marker": "p",  "ms": 10,  "elinewidth": 3},
    "finkelstein15":{"color": "deepskyblue",    "marker": "*",  "ms": 10, "elinewidth": 3},
    "bouwens21":   {"color": "pink",     "marker": "X",  "ms": 10,  "elinewidth": 3},
    "bowler17":    {"color": "magenta",     "marker": "h",  "ms": 10,  "elinewidth": 3},
    "harikane24":  {"color": "black",    "marker": "H",  "ms": 10,  "elinewidth": 3},
    "mcleod24":    {"color": "coral",    "marker": ">",  "ms": 10,  "elinewidth": 3},
    "chemerynska24":    {"color": "steelblue",    "marker": "P",  "ms": 10,  "elinewidth": 3},
    "whitler25":    {"color": "plum",    "marker": "8",  "ms": 10,  "elinewidth": 3},
    "Adams24":     {"color": "gold",    "marker": "h",  "ms": 10,  "elinewidth": 3},
    "Perez-Gonzalez23": {"color": "dodgerblue", "marker":"<", "ms": 10, "elinewidth": 3}
}

#! Plot Phi*
plot_data(axs[0, 0], my_z, my_phi, [my_phi_err_lo, my_phi_err_up], "This work", styles["my"], zorder=10)
plot_data(axs[0, 0], donnan23_z, donnan23_phi, donnan23_phi_err, "Donnan+23", styles["donnan23"], zorder=zorder)
plot_data(axs[0, 0], donnan24_z, donnan24_phi, donnan24_phi_err, "Donnan+24", styles["donnan24"], zorder=zorder)
plot_data(axs[0, 0], bowler20_z, bowler20_phi, bowler20_phi_err, "Bowler+20", styles["bowler20"], zorder=zorder)
plot_data(axs[0, 0], bowler15_z, bowler15_phi, [bowler15_phi_err_lo, bowler15_phi_err_up], "Bowler+15", styles["bowler15"], zorder=zorder)
plot_data(axs[0, 0], adams23_z, adams23_phi, adams23_phi_err, "Adams+23", styles["adams23"], zorder=zorder)
#plot_data(axs[0, 0], finkelstein15_z, finkelstein15_phi, [finkelstein15_phi_err_lo, finkelstein15_phi_err_up], "Finkelstein+15", styles["finkelstein15"], zorder=zorder)
#plot_data(axs[0, 0], bouwens21_z, bouwens21_phi, [bouwens21_phi_err_lo, bouwens21_phi_err_up], "Bouwens+21", styles["bouwens21"], zorder=zorder)
plot_data(axs[0, 0], bowler17_z, bowler17_phi, [bowler17_phi_err_lo, bowler17_phi_err_up], "Bowler+17", styles["bowler17"], zorder=zorder)
plot_data(axs[0, 0], harikane24_z, harikane24_phi, [harikane24_phi_err_lo, harikane24_phi_err_up], "Harikane+24", styles["harikane24"], zorder=zorder)
plot_data(axs[0, 0], mcleod24_z, mcleod24_phi, mcleod24_phi_err, "McLeod+24", styles["mcleod24"], zorder=zorder)
plot_data(axs[0, 0], chemerynska24_z, chemerynska24_phi, chemerynska24_phi_err, "Chemerynska+24", styles["chemerynska24"], zorder=zorder)
plot_data(axs[0, 0], whitler25_z, whitler25_phi, [whitler25_phi_err_lo, whitler25_phi_err_up], "Whitler+25", styles["whitler25"], zorder=zorder)
plot_data(axs[0, 0], adams24_z, adams24_phi, [adams24_phi_err_lo, adams24_phi_err_up], "Adams+24", styles["Adams24"], zorder=zorder)
#plot_data(axs[0, 0], perez23_z, perez23_phi, [perez23_phi_err_lo, perez23_phi_err_up], "Perez-Gonzalez+23", styles["Perez-Gonzalez23"], zorder=zorder)

axs[0, 0].set_xlabel('z', fontsize=30)
axs[0, 0].set_ylabel(r'$\log_{10}(\phi^{*}\,/\,\rm{mag}^{-1}\,\rm{Mpc}^{-3})$')
axs[0, 0].set_xlim(5.5, 13)
axs[0, 0].set_ylim(-6., -2.7)



#! Plot M*
plot_data(axs[0, 1], my_z, my_M_star, [my_M_star_err_lo, my_M_star_err_up], "This work", styles["my"], zorder=10)
plot_data(axs[0, 1], donnan23_z, donnan23_M_star, donnan23_M_star_err, "Donnan+23", styles["donnan23"], zorder=zorder)
plot_data(axs[0, 1], donnan24_z, donnan24_M_star, donnan24_M_star_err, "Donnan+24", styles["donnan24"], zorder=zorder)
plot_data(axs[0, 1], bowler20_z, bowler20_M_star, bowler20_M_star_err, "Bowler+20", styles["bowler20"], zorder=zorder)
plot_data(axs[0, 1], bowler15_z, bowler15_M_star, bowler15_M_star_err, "Bowler+15", styles["bowler15"], zorder=zorder)
plot_data(axs[0, 1], adams23_z, adams23_M_star, adams23_M_star_err, "Adams+23", styles["adams23"], zorder=zorder)
#plot_data(axs[0, 1], finkelstein15_z, finkelstein15_M_star, [finkelstein15_M_star_err_lo, finkelstein15_M_star_err_up], "Finkelstein+15", styles["finkelstein15"], zorder=zorder)
#plot_data(axs[0, 1], bouwens21_z, bouwens21_M_star, bouwens21_M_star_err, "Bouwens+21", styles["bouwens21"], zorder=zorder)
plot_data(axs[0, 1], bowler17_z, bowler17_M_star, [bowler17_M_star_err_lo, bowler17_M_star_err_up], "Bowler+17", styles["bowler17"], zorder=zorder)
plot_data(axs[0, 1], harikane24_z, harikane24_M_star, [harikane24_M_star_err_lo, harikane24_M_star_err_up], "Harikane+24", styles["harikane24"], zorder=zorder)
plot_data(axs[0, 1], mcleod24_z, mcleod24_M_star, mcleod24_M_star_err, "McLeod+24", styles["mcleod24"], zorder=zorder)
plot_data(axs[0, 1], chemerynska24_z, chemerynska24_M_star, chemerynska24_M_star_err, "Chemerynska+24", styles["chemerynska24"], zorder=zorder)
plot_data(axs[0, 1], whitler25_z, whitler25_M_star, [whitler25_M_star_err_lo, whitler25_M_star_err_up], "Whitler+25", styles["whitler25"], zorder=zorder)
plot_data(axs[0, 1], adams24_z, adams24_M_star, [adams24_M_star_err_lo, adams24_M_star_err_up], "Adams+24", styles["Adams24"], zorder=zorder)
#plot_data(axs[0, 1], perez23_z, perez23_M_star, [perez23_M_star_err_lo, perez23_M_star_err_up], "Perez-Gonzalez+23", styles["Perez-Gonzalez23"], zorder=zorder)
axs[0, 1].set_xlabel('z', fontsize=30)
axs[0, 1].set_ylabel(r'$M^{*}$')
axs[0, 1].set_xlim(5.5, 13)
axs[0, 1].set_ylim(-23, -18.8)

#! Plot Alpha
plot_data(axs[1, 0], my_z, my_alpha, [my_alpha_err_lo, my_alpha_err_up], "This work", styles["my"], zorder=10)
plot_data(axs[1, 0], donnan23_z, donnan23_alpha, donnan23_alpha_err, "Donnan+23", styles["donnan23"], zorder=zorder)
plot_data(axs[1, 0], donnan24_z, donnan24_alpha, donnan24_alpha_err, "Donnan+24", styles["donnan24"], zorder=zorder)
plot_data(axs[1, 0], bowler20_z, bowler20_alpha, bowler20_alpha_err, "Bowler+20", styles["bowler20"], zorder=zorder)
plot_data(axs[1, 0], bowler15_z, bowler15_alpha, [bowler15_alpha_err_lo, bowler15_alpha_err_up], "Bowler+15", styles["bowler15"], zorder=zorder)
plot_data(axs[1, 0], adams23_z, adams23_alpha, adams23_alpha_err, "Adams+23", styles["adams23"], zorder=zorder)
#plot_data(axs[1, 0], finkelstein15_z, finkelstein15_alpha, [finkelstein15_alpha_err_lo, finkelstein15_alpha_err_up], "Finkelstein+15", styles["finkelstein15"], zorder=zorder)
#plot_data(axs[1, 0], bouwens21_z, bouwens21_alpha, bouwens21_alpha_err, "Bouwens+21", styles["bouwens21"], zorder=zorder)
plot_data(axs[1, 0], bowler17_z, bowler17_alpha, [bowler17_alpha_err_lo, bowler17_alpha_err_up], "Bowler+17", styles["bowler17"], zorder=zorder)
plot_data(axs[1, 0], harikane24_z, harikane24_alpha, [harikane24_alpha_err_lo, harikane24_alpha_err_up], "Harikane+24", styles["harikane24"], zorder=zorder)
plot_data(axs[1, 0], mcleod24_z, mcleod24_alpha, mcleod24_alpha_err, "McLeod+24", styles["mcleod24"], zorder=zorder)
plot_data(axs[1, 0], chemerynska24_z, chemerynska24_alpha, chemerynska24_alpha_err, "Chemerynska+24", styles["chemerynska24"], zorder=zorder)
plot_data(axs[1, 0], whitler25_z, whitler25_alpha, [whitler25_alpha_err_lo, whitler25_alpha_err_up], "Whitler+25", styles["whitler25"], zorder=zorder)
plot_data(axs[1, 0], adams24_z, adams24_alpha, [adams24_alpha_err_lo, adams24_alpha_err_up], "Adams+24", styles["Adams24"], zorder=zorder)
#plot_data(axs[1, 0], perez23_z, perez23_alpha, [perez23_alpha_err_lo, perez23_alpha_err_up], "Perez-Gonzalez+23", styles["Perez-Gonzalez23"], zorder=zorder)
axs[1, 0].set_xlabel('z', fontsize=30)
axs[1, 0].set_xlim(5.5, 13)
axs[1, 0].set_ylabel(r'$\alpha$')
axs[1, 0].set_ylim(-3, -1.4)

#axs[1, 1].legend(loc='lower left', fontsize=15)

#! Plot Beta
plot_data(axs[1, 1], my_z, my_beta, [my_beta_err_lo, my_beta_err_up], "This work", styles["my"], zorder=10)
plot_data(axs[1, 1], donnan23_z, donnan23_beta, donnan23_beta_err, "Donnan+23", styles["donnan23"], zorder=zorder)
plot_data(axs[1, 1], donnan24_z, donnan24_beta, donnan24_beta_err, "Donnan+24", styles["donnan24"], zorder=zorder)
plot_data(axs[1, 1], bowler20_z, bowler20_beta, bowler20_beta_err, "Bowler+20", styles["bowler20"], zorder=zorder)
plot_data(axs[1, 1], bowler15_z, bowler15_beta, [bowler15_beta_err_lo, bowler15_beta_err_up], "Bowler+15", styles["bowler15"], zorder=zorder)
plot_data(axs[1, 1], adams23_z, adams23_beta, [adams23_beta_err_lo, adams23_beta_err_up], "Adams+23", styles["adams23"], zorder=zorder)
plot_data(axs[1, 1], bowler17_z, bowler17_beta, [bowler17_beta_err_lo, bowler17_beta_err_up], "Bowler+17", styles["bowler17"], zorder=zorder)
plot_data(axs[1, 1], harikane24_z, harikane24_beta, [harikane24_beta_err_lo, harikane24_beta_err_up], "Harikane+24", styles["harikane24"], zorder=zorder)
#plot_data(axs[1, 1], mcleod24_z, mcleod24_beta, mcleod24_beta_err, "McLeod+24", styles["mcleod24"], zorder=zorder)
plot_data(axs[1, 1], chemerynska24_z, chemerynska24_beta, chemerynska24_beta_err, "Chemerynska+24", styles["chemerynska24"], zorder=zorder)
#plot_data(axs[1, 1], whitler25_z, whitler25_beta, [whitler25_beta_err_lo, whitler25_beta_err_up], "Whitler+25", styles["whitler25"], zorder=zorder)
plot_data(axs[1, 1], adams24_z[:-1], adams24_beta[:-1], [adams24_beta_err_lo[:-1], adams24_beta_err_up[:-1]], "Adams+24", styles["Adams24"], zorder=zorder)
axs[1, 1].set_xlabel('z', fontsize=30)
axs[1, 1].set_xlim(5.5, 13)
axs[1, 1].set_ylabel(r'$\beta$')
axs[1, 1].set_ylim(-6.8, -1.5)


########! FITTING #############
def linear(x, a, b):
    return a + b*(x-6)



def power_law(x, A, alpha):
    return A * x ** alpha




def fit_and_plot_linear(ax, redshift_data, param_data, param_err_lo, param_err_up, label, color='black', linewidth=4):
    """
    Fits a straight line to the given parameter data within z = 6 to 11, 
    excluding data points with zero error, and plots the fit with uncertainty.

    Parameters:
    - ax: The axis on which to plot
    - redshift_data: List of redshift arrays from different studies
    - param_data: List of parameter arrays (e.g., beta, M*, alpha, etc.)
    - param_err_lo: List of lower errors on parameter values
    - param_err_up: List of upper errors on parameter values
    - label: Label for the fit line
    - color: Color of the fit line (default: black)
    - linewidth: Line width of the fit (default: 4)
    """

    # Concatenate all datasets
    redshifts = np.concatenate(redshift_data)
    params = np.concatenate(param_data)
    errs_lo = np.concatenate(param_err_lo)
    errs_up = np.concatenate(param_err_up)

    # Apply mask: exclude points with zero error
    mask = (errs_lo != 0) & (errs_up != 0)
    redshifts = redshifts[mask]
    params = params[mask]
    errs = np.array([errs_lo[mask], errs_up[mask]])

    # Take the maximum of the two errors for fitting
    errs = np.max(errs, axis=0)

    # Fit a straight line
    popt, pcov = curve_fit(linear, redshifts, params, sigma=errs, absolute_sigma=True)

    # Generate fit line
    z_fit = np.linspace(4.6, 14, 100)
    param_fit = linear(z_fit, *popt)
    ax.plot(z_fit, param_fit, color=color, label=label, linewidth=linewidth)

    print('Best fit parameters:', popt) 

    # Get errors from covariance matrix
    perr = np.sqrt(np.diag(pcov))
    print('Errors on parameters:', perr)

    # Monte Carlo sampling for error band
    num_samples = 1000
    param_samples = np.random.multivariate_normal(popt, pcov, num_samples)
    param_fits = np.array([linear(z_fit, *params) for params in param_samples])

    # Compute confidence intervals (16th-84th percentile)
    lower_bound = np.percentile(param_fits, 16, axis=0)
    upper_bound = np.percentile(param_fits, 84, axis=0)

    # Plot shaded uncertainty region
    ax.fill_between(z_fit, lower_bound, upper_bound, color=color, alpha=0.3)

    # Write equation of the line in LaTeX format with errors
    eqn = f"$y = ({popt[0]:.2f} \pm {perr[0]:.2f}) + ({popt[1]:.2f} \pm {perr[1]:.2f})(z-6)$"
    print(eqn)

    return popt  # Returns fit parameters if needed
# * BETA
beta_fit = fit_and_plot_linear(
    axs[1, 1], 
    [donnan23_z, donnan24_z, bowler20_z, bowler15_z, bowler17_z, harikane24_z, chemerynska24_z, my_z, adams24_z[:-1]], 
    [donnan23_beta, donnan24_beta, bowler20_beta, bowler15_beta, bowler17_beta, harikane24_beta, chemerynska24_beta, my_beta, adams24_beta[:-1]], 
    [donnan23_beta_err, donnan24_beta_err, bowler20_beta_err, bowler15_beta_err_lo, bowler17_beta_err_lo, harikane24_beta_err_lo, chemerynska24_beta_err, my_beta_err_lo, adams24_beta_err_lo[:-1]], 
    [donnan23_beta_err, donnan24_beta_err, bowler20_beta_err, bowler15_beta_err_up, bowler17_beta_err_up, harikane24_beta_err_up, chemerynska24_beta_err, my_beta_err_up, adams24_beta_err_up[:-1]], 
    label=None
)


# For alpha, phi*, M*, restrict bouwens and finkelstein to z>6
min_z = 6

bouwens21_alpha = bouwens21_alpha[bouwens21_z > min_z]
bouwens21_alpha_err = bouwens21_alpha_err[bouwens21_z > min_z]
bouwens21_phi = bouwens21_phi[bouwens21_z > min_z]
bouwens21_phi_err_lo = bouwens21_phi_err_lo[bouwens21_z > min_z]
bouwens21_phi_err_up = bouwens21_phi_err_up[bouwens21_z > min_z]
bouwens21_M_star = bouwens21_M_star[bouwens21_z > min_z]
bouwens21_M_star_err = bouwens21_M_star_err[bouwens21_z > min_z]
bouwens21_z = bouwens21_z[bouwens21_z > min_z]

finkelstein15_alpha = finkelstein15_alpha[finkelstein15_z > min_z]
finkelstein15_alpha_err_up = finkelstein15_alpha_err_up[finkelstein15_z > min_z]
finkelstein15_alpha_err_lo = finkelstein15_alpha_err_lo[finkelstein15_z > min_z]
finkelstein15_phi = finkelstein15_phi[finkelstein15_z > min_z]
finkelstein15_phi_err_lo = finkelstein15_phi_err_lo[finkelstein15_z > min_z]
finkelstein15_phi_err_up = finkelstein15_phi_err_up[finkelstein15_z > min_z]
finkelstein15_M_star = finkelstein15_M_star[finkelstein15_z > min_z]
finkelstein15_M_star_err_up = finkelstein15_M_star_err_up[finkelstein15_z > min_z]
finkelstein15_M_star_err_lo = finkelstein15_M_star_err_lo[finkelstein15_z > min_z]
finkelstein15_z = finkelstein15_z[finkelstein15_z > min_z]


# * ALPHA
alpha_fit = fit_and_plot_linear(
    axs[1, 0],
    [donnan23_z, donnan24_z, bowler20_z, bowler15_z, bowler17_z, harikane24_z, mcleod24_z, chemerynska24_z, whitler25_z, my_z, adams24_z],
    [donnan23_alpha, donnan24_alpha, bowler20_alpha, bowler15_alpha, bowler17_alpha, harikane24_alpha, mcleod24_alpha, chemerynska24_alpha, whitler25_alpha, my_alpha, adams24_alpha],
    [donnan23_alpha_err, donnan24_alpha_err, bowler20_alpha_err, bowler15_alpha_err_lo, bowler17_alpha_err_lo, harikane24_alpha_err_lo, mcleod24_alpha_err, chemerynska24_alpha_err, whitler25_alpha_err_lo, my_alpha_err_lo, adams24_alpha_err_lo],
    [donnan23_alpha_err, donnan24_alpha_err, bowler20_alpha_err, bowler15_alpha_err_up, bowler17_alpha_err_up, harikane24_alpha_err_up, mcleod24_alpha_err, chemerynska24_alpha_err, whitler25_alpha_err_up, my_alpha_err_up, adams24_alpha_err_up],
    label=None
)

#print(donnan23_phi_err.shape, donnan24_phi_err.shape, bowler20_phi_err.shape, bowler15_phi_err_lo.shape, bowler17_phi_err_lo.shape, harikane24_phi_err.shape, mcleod24_phi_err.shape, chemerynska24_phi_err.shape, whitler25_phi_err_lo.shape, my_phi_err_lo.shape)

# Split the harikane24 errors into upper and lower. First element is lower, second is upper, of each entry
# harikane24_phi_err_up = harikane24_phi_err[1]
# harikane24_phi_err_lo = harikane24_phi_err[0]

# # Same with chemerynska24
# chemerynska24_phi_err_up = chemerynska24_phi_err[1]
# chemerynska24_phi_err_lo = chemerynska24_phi_err[0]

print('PHI')
def fit_and_plot_log_phi(ax, redshift_data, param_data, param_err_lo, param_err_up, label, color='black', linewidth=4):
    """
    Fits a straight line to the given parameter data within z = 6 to 11, 
    excluding data points with zero error, and plots the fit.

    Parameters:
    - ax: The axis on which to plot
    - redshift_data: List of redshift arrays from different studies
    - param_data: List of parameter arrays (e.g., beta, M*, alpha, etc.)
    - param_err_lo: List of lower errors on parameter values
    - param_err_up: List of upper errors on parameter values
    - label: Label for the fit line
    - color: Color of the fit line (default: black)
    - linewidth: Line width of the fit (default: 4)
    """

    # Concatenate all datasets
    redshifts = np.concatenate(redshift_data)
    params = np.concatenate(param_data)
    errs_lo = np.concatenate(param_err_lo)
    errs_up = np.concatenate(param_err_up)

    # Apply mask: exclude points with zero error
    mask = (errs_lo != 0) & (errs_up != 0)
    redshifts = redshifts[mask]
    params = params[mask]
    errs = np.array([errs_lo[mask], errs_up[mask]])

    # Take the maximum of the two errors for fitting
    errs = np.max(errs, axis=0)

    params = np.log10(params)

    # Fit a straight line
    popt, _ = curve_fit(linear, redshifts, params, sigma=errs, absolute_sigma=True)

    # Generate fit line
    z_fit = np.linspace(4.6, 14, 100)
    ax.plot(z_fit, linear(z_fit, *popt), color=color, label=label, linewidth=linewidth)
    print('Best fit parameters:', popt) 
    return popt  # Returns fit parameters if needed

phi_fit = fit_and_plot_linear(
    axs[0, 0],
    [donnan23_z, donnan24_z, bowler20_z, bowler15_z, bowler17_z, harikane24_z, mcleod24_z, chemerynska24_z, whitler25_z, my_z, adams24_z],
    [donnan23_phi, donnan24_phi, bowler20_phi, bowler15_phi, bowler17_phi, harikane24_phi, mcleod24_phi, chemerynska24_phi, whitler25_phi, my_phi, adams24_phi],
    [donnan23_phi_err, donnan24_phi_err, bowler20_phi_err, bowler15_phi_err_lo, bowler17_phi_err_lo, harikane24_phi_err_lo, mcleod24_phi_err, chemerynska24_phi_err, whitler25_phi_err_lo, my_phi_err_lo, adams24_phi_err_lo],
    [donnan23_phi_err, donnan24_phi_err, bowler20_phi_err, bowler15_phi_err_up, bowler17_phi_err_up, harikane24_phi_err_up, mcleod24_phi_err, chemerynska24_phi_err, whitler25_phi_err_up, my_phi_err_up, adams24_phi_err_up],
    label=None
)




#* M*
M_fit = fit_and_plot_linear(
    axs[0, 1],
    [donnan23_z, donnan24_z, bowler15_z, bowler17_z, harikane24_z, mcleod24_z, chemerynska24_z, my_z, adams24_z, whitler25_z, bowler20_z],
    [donnan23_M_star, donnan24_M_star, bowler15_M_star, bowler17_M_star, harikane24_M_star, mcleod24_M_star, chemerynska24_M_star, my_M_star, adams24_M_star, whitler25_M_star, bowler20_M_star],
    [donnan23_M_star_err, donnan24_M_star_err, bowler15_M_star_err, bowler17_M_star_err_lo, harikane24_M_star_err_lo, mcleod24_M_star_err, chemerynska24_M_star_err, my_M_star_err_lo, adams24_M_star_err_lo, whitler25_M_star_err_lo, bowler20_M_star_err],
    [donnan23_M_star_err, donnan24_M_star_err, bowler15_M_star_err, bowler17_M_star_err_up, harikane24_M_star_err_up, mcleod24_M_star_err, chemerynska24_M_star_err, my_M_star_err_up, adams24_M_star_err_up, whitler25_M_star_err_up, bowler20_M_star_err],
    label=None
)


##! FINAL PLOT DETAILS
# Add legend below the whole plot
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=25, bbox_to_anchor=(0.5, -0.01))

# Axes minor tickson
axs[0, 0].minorticks_on()
axs[0, 1].minorticks_on()
axs[1, 0].minorticks_on()
axs[1, 1].minorticks_on()

fig.subplots_adjust(bottom=0.25)

fig.tight_layout(rect=[0, 0.12, 1, 1])

#plt.tight_layout()
plt.savefig(plot_dir / 'LF_param_evolution.pdf')
