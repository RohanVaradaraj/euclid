"""
Generate BAGPIPES model galaxies at a fixed redshift over a grid of delayed-tau
star-formation-history parameters, then save SFR, stellar mass and H-alpha
line fluxes into an astropy Table.

Assumptions & notes:
- Uses Bagpipes (pip install bagpipes) to create model spectra and emission-line
  fluxes. The model_galaxy class exposes emission-line fluxes in
  `model.line_fluxes` (units: erg/s/cm^2) when a `nebular` component is
  included in model_components.
- The delayed-tau SFH used has the analytic form: SFR(t) = A * t * exp(-t/tau).
  `massformed` (in BAGPIPES) is the total mass formed up to `age` and must be
  provided as log10(M_solar). We compute the instantaneous SFR at the epoch
  of observation analytically from the massformed and SFH parameters.
- The script provides a simple estimate of the *living* stellar mass by
  applying a fixed return fraction (default: 0.30). If you want the exact
  living stellar mass from the SPS model you should extract derived
  quantities from a full BAGPIPES fit; this script is intended for fast
  forward-model grid generation.

Outputs:
- Writes an Astropy Table to disk (FITS and ASCII) with the following columns:
  redshift, tau_gyr, age_gyr, log_massformed, metallicity, sfr_inst_msun_per_yr,
  stellar_mass_living_msun, halpha_flux_erg_s_cm2 (NaN if nebular disabled or
  H-alpha not present).

Usage:
    python bagpipes_delayed_tau_grid.py

You can edit the grid parameters in the `CONFIG` block below.
"""

import numpy as np
import itertools
import argparse
import warnings
from tqdm import tqdm # progress bar
import bagpipes as pipes
import os

from astropy.table import Table, Column


# ----------------------------- CONFIG ------------------------------------
overwrite = True  # overwrite existing output files

# Fixed redshift for all models
redshifts = [6.52, 5.38, 5.1157, 4.76, 4.72, 4.64] #, 4.4899]
REDSHIFT = 4.64

# # ID of source with this redshift
# ID = 89279

# Nebular ionization parameter (log10 U)
LOGU = -2.0
# Default return fraction (mass returned to ISM) -> living mass = massformed * (1 - return_fraction)
DEFAULT_RETURN_FRACTION = 0.30
# Grids (you can change these lists)
TAU_GYR_LIST = np.array([0.1, 0.3, 0.5, 1.0, 3.0])        # tau in Gyr
AGE_GYR_LIST = np.array([0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0])   # age in Gyr (time since SF began)
LOG_MASSFORMED_LIST = np.arange(8.0, 11.5, 0.25)  # log10 total mass formed in M_sun
METALLICITY_LIST = np.array([0.02])    # Z (absolute; Z_sun ~ 0.02)
DUST_LIST = np.array([0.5]) #, 1.0, 2.0]) #, 3.0])
# Output files
OUT_BASENAME = "bagpipes_delayed_tau_grid_z{:.2f}".format(REDSHIFT)
# -------------------------------------------------------------------------

filt_lis = np.loadtxt('z6_filt_list.txt', dtype='str')

# Analytical functions for delayed-tau SFH
def delayed_tau_integral(age_gyr, tau_gyr):
    """Return integral I = \int_0^{age} t exp(-t/tau) dt
    with t in Gyr. Analytic result: tau^2 * (1 - (1 + age/tau) * exp(-age/tau)).
    """
    age = age_gyr
    tau = tau_gyr
    return tau**2.0 * (1.0 - (1.0 + age / tau) * np.exp(-age / tau))


def instantaneous_sfr_from_massformed(age_gyr, tau_gyr, massformed_total_msun):
    """Compute instantaneous SFR at t=age (M_sun / yr) for a delayed-tau SFH
    normalized to a total `massformed_total_msun` formed between t=0 and t=age.

    SFR(t) = A * t * exp(-t/tau)
    massformed = A * integral => A = massformed / integral
    SFR(age) = A * age * exp(-age/tau)

    Note: age and tau are in Gyr; convert SFR to M_sun/yr (divide by 1e9).
    """
    I = delayed_tau_integral(age_gyr, tau_gyr)
    if I <= 0:
        return 0.0
    A = massformed_total_msun / I
    sfr_gyr = A * age_gyr * np.exp(-age_gyr / tau_gyr)   # M_sun per Gyr
    sfr_yr = sfr_gyr / 1e9
    return sfr_yr


def find_halpha_key(line_fluxes_dict):
    """Return the key in bagpipes line_fluxes dict that corresponds to H-alpha if present.
    Bagpipes names lines similarly to Cloudy; H I 6562.81A or similar. We'll search
    for the substring '656' or 'H  1' together with '656'. Return None if not found.
    """
    for k in line_fluxes_dict.keys():
        if "656" in k.replace(" ", "") or "656" in k:
            # extra safeguard: require 'H' or 'H  1' in key (ignore numeric-only matches)
            if 'H' in k or 'H 1' in k or 'H  1' in k:
                return k
    # fallback: try keys containing 'H' and '656' separated by whitespace
    for k in line_fluxes_dict.keys():
        if 'H' in k and '656' in k:
            return k
    return None

def find_emission_line_keys(line_fluxes_dict):
    """
    Return the keys corresponding to key emission lines (H-alpha, [O III]) if present.
    
    Bagpipes line names typically follow Cloudy conventions, e.g.:
        'H  1 6562.81A'  for Hα
        'O  3 5006.84A'  for [O III]
    
    Returns
    -------
    lines_found : dict
        Dictionary with keys 'Halpha' and 'OIII5007' (if found), e.g.:
        {'Halpha': 'H  1 6562.81A', 'OIII5007': 'O  3 5006.84A'}
    """
    lines_found = {}

    # Normalize keys for easier searching
    for k in line_fluxes_dict.keys():
        k_clean = k.replace(" ", "").upper()

        # --- H-alpha (H I 6563 Å) ---
        if ("H" in k_clean or "HI" in k_clean) and "656" in k_clean:
            lines_found["Halpha"] = k

        # --- [O III] λ5007 ---
        elif ("OIII" in k_clean or "O3" in k_clean or "O" in k_clean) and (
            "5006" in k_clean or "5007" in k_clean
        ):
            lines_found["OIII5007"] = k

    return lines_found if lines_found else None

def make_model_and_extract(redshift, tau_gyr, age_gyr, log_massformed, metallicity, dust, logU=LOGU):
    """Build a BAGPIPES model_galaxy for a single delayed-tau parameter set and
    extract H-alpha line flux (if present). Returns (sfr_inst, stellar_mass_living, halpha_flux).

    Note: `massformed` in bagpipes is *log10* total mass formed in M_sun.
    """
    # Prepare model_components dict
    model_components = {}
    model_components['redshift'] = redshift

    # SFH component named e.g. 'delayed1'
    delayed = {}
    delayed['massformed'] = float(log_massformed)   # BAGPIPES expects log10 mass formed
    delayed['metallicity'] = float(metallicity)     # in absolute Z (0.02 ~= Z_sun)
    delayed['age'] = float(age_gyr)
    delayed['tau'] = float(tau_gyr)

    # Insert the sfh component; bagpipes expects the component dict under a key name
    model_components['delayed1'] = delayed

    # Nebular emission component to get lines
    nebular = {'logU': float(logU)}
    model_components['nebular'] = nebular

    # Minimal dust (required by some BAGPIPES versions) - use Calzetti with small Av
    dust = {'type': 'Calzetti', 'Av': dust}
    model_components['dust'] = dust

    # Create model_galaxy. We do not request photometry or spectroscopy sampling here,
    # only emission line fluxes (which are produced internally when nebular is present).
    model = pipes.model_galaxy(model_components, filt_list=filt_lis)

    # Extract H-alpha if available
    # halpha_flux = np.nan
    # try:
    #     lf = model.line_fluxes
    #     if lf is not None and len(lf) > 0:
    #         key = find_halpha_key(lf)
    #         if key is not None:
    #             halpha_flux = float(lf[key])  # erg / s / cm^2
    #         else:
    #             halpha_flux = np.nan
    #     else:
    #         halpha_flux = np.nan
    # except Exception:
    #     # Some bagpipes builds may not populate line_fluxes unless certain grids are present.
    #     halpha_flux = np.nan

    # Extract H-alpha and [O III] λ5007 if available
    halpha_flux = np.nan
    oiii_flux = np.nan

    try:
        lf = model.line_fluxes
        if lf is not None and len(lf) > 0:
            line_keys = find_emission_line_keys(lf)

            if line_keys is not None:
                if "Halpha" in line_keys:
                    halpha_flux = float(lf[line_keys["Halpha"]])  # erg/s/cm^2
                if "OIII5007" in line_keys:
                    oiii_flux = float(lf[line_keys["OIII5007"]])  # erg/s/cm^2
        # else: keep NaNs if no line_fluxes
    except Exception as e:
        # Some bagpipes builds may not populate line_fluxes unless certain grids are present.
        halpha_flux, oiii_flux = np.nan, np.nan

    # Compute instantaneous SFR analytically from SFH normalization
    massformed_total = 10.0 ** float(log_massformed)  # M_sun
    sfr_inst = instantaneous_sfr_from_massformed(age_gyr, tau_gyr, massformed_total)

    # Estimate living stellar mass using return fraction
    stellar_mass_living = massformed_total * (1.0 - DEFAULT_RETURN_FRACTION)

    return sfr_inst, stellar_mass_living, halpha_flux, oiii_flux


def main():
    import os
    fits_out = OUT_BASENAME + ".fits"
    ascii_out = OUT_BASENAME + ".ascii"

    if overwrite:
        os.remove(fits_out) if os.path.exists(fits_out) else None
        os.remove(ascii_out) if os.path.exists(ascii_out) else None

    files_exist = os.path.exists(fits_out) and os.path.exists(ascii_out)

    if files_exist:
        print(f"Output files already exist: {fits_out}, {ascii_out}. Skipping generation.")
    else:
        combos = list(itertools.product(TAU_GYR_LIST, AGE_GYR_LIST, LOG_MASSFORMED_LIST, METALLICITY_LIST, DUST_LIST))
        rows = []

        from tqdm import tqdm
        for tau_gyr, age_gyr, log_mf, Z, dust_Av in tqdm(combos, desc="Generating models", unit="model"):
            try:
                sfr_inst, stellar_mass_living, halpha_flux, oiii_flux = make_model_and_extract(REDSHIFT, tau_gyr, age_gyr, log_mf, Z, dust_Av)
            except Exception as e:
                warnings.warn(f"Failed for tau={tau_gyr}, age={age_gyr}, logM={log_mf}, Z={Z}: {e}")
                sfr_inst, stellar_mass_living, halpha_flux, oiii_flux = np.nan, np.nan, np.nan, np.nan
            rows.append((REDSHIFT, tau_gyr, age_gyr, log_mf, Z, sfr_inst, stellar_mass_living, halpha_flux, oiii_flux))

        t = Table(rows=rows, names=("redshift", "tau_gyr", "age_gyr", "log_massformed", "metallicity",
                                    "sfr_inst_msun_per_yr", "stellar_mass_living_msun", "halpha_flux_erg_s_cm2", "oiii_flux_erg_s_cm2"))

        t.write(fits_out, overwrite=True)
        t.write(ascii_out, format='ascii.fixed_width', overwrite=True)
        print(f"Wrote {len(t)} models to {fits_out} and {ascii_out}")

    return fits_out  # return filename so plotting can always use it




def plot_halpha_2d(filename, xcol, ycol, logx=False, logy=False):
    """Read an Astropy table and make a 2D colour plot of H-alpha flux."""
    import matplotlib.pyplot as plt
    from astropy.table import Table
    import numpy as np

    # Read table
    t = Table.read(filename)
    if "halpha_flux_erg_s_cm2" not in t.colnames:
        raise ValueError("Table must contain 'halpha_flux_erg_s_cm2' column")

    # Extract requested axes
    x = np.array(t[xcol])
    y = np.array(t[ycol])
    halpha = np.array(t["halpha_flux_erg_s_cm2"])
    halpha = np.where(halpha > 0, halpha, np.nan)  # mask invalids

    c = np.log10(halpha)  # log scale for visibility

    plt.figure(figsize=(7, 5))
    sc = plt.scatter(x, y, c=c, cmap="viridis", s=80, edgecolor="k", lw=0.3)
    plt.colorbar(sc, label=r"log$_{10}$(H$\alpha$ flux [erg s$^{-1}$ cm$^{-2}$])")
    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.title(f"Hα flux colour map: {ycol} vs {xcol}")
    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename.replace(".fits", f"_halpha_{xcol}_vs_{ycol}.pdf"), dpi=100)
    #plt.show()

def plot_halpha_vs_sfr(filename, logx=True, logy=True, logc=True):
    """
    Plot H-alpha flux vs instantaneous SFR, colour-coded by stellar mass.
    
    Parameters
    ----------
    filename : str
        Astropy table FITS file from BAGPIPES grid.
    logx : bool
        Use log scale for SFR axis.
    logy : bool
        Use log scale for H-alpha flux axis.
    logc : bool
        Use log scale for colour (stellar mass).
    """
    import matplotlib.pyplot as plt
    from astropy.table import Table
    import numpy as np

    t = Table.read(filename)
    for col in ["halpha_flux_erg_s_cm2", "sfr_inst_msun_per_yr", "stellar_mass_living_msun"]:
        if col not in t.colnames:
            raise ValueError(f"Table must contain '{col}' column")

    x = np.array(t["sfr_inst_msun_per_yr"])
    y = np.array(t["halpha_flux_erg_s_cm2"])
    c = np.array(t["stellar_mass_living_msun"])

    # Mask invalid or non-positive values for log scale
    x_plot = np.where(x>0, x, np.nan)
    y_plot = np.where(y>0, y, np.nan)
    c_plot = np.where(c>0, c, np.nan)

    if logc:
        c_plot = np.log10(c_plot)

    plt.figure(figsize=(8,6))
    sc = plt.scatter(x_plot, y_plot, c=c_plot, cmap="viridis", s=80, edgecolor='k', lw=0.3)
    cbar = plt.colorbar(sc)
    cbar.set_label("log10(Stellar mass [M_sun])" if logc else "Stellar mass [M_sun]")

    plt.xlabel("Instantaneous SFR [M_sun / yr]")
    plt.ylabel("H-alpha flux [erg / s / cm^2]")
    plt.title("H-alpha flux vs SFR, coloured by stellar mass")

    if logx:
        plt.xscale("log")
    if logy:
        plt.yscale("log")

    plt.grid(alpha=0.3)
    plt.tight_layout()
    outname = filename.replace(".fits", "_halpha_vs_sfr.pdf")
    plt.savefig(outname, dpi=150)
    #plt.show()
    print(f"Saved plot to {outname}")

    # Repeat with Oiii
    x = np.array(t["sfr_inst_msun_per_yr"])
    y = np.array(t["oiii_flux_erg_s_cm2"])
    c = np.array(t["stellar_mass_living_msun"])

    # Mask invalid or non-positive values for log scale
    x_plot = np.where(x>0, x, np.nan)
    y_plot = np.where(y>0, y, np.nan)
    c_plot = np.where(c>0, c, np.nan)

    if logc:
        c_plot = np.log10(c_plot)

    plt.figure(figsize=(8,6))
    sc = plt.scatter(x_plot, y_plot, c=c_plot, cmap="viridis", s=80, edgecolor='k', lw=0.3)
    cbar = plt.colorbar(sc)
    cbar.set_label("log10(Stellar mass [M_sun])" if logc else "Stellar mass [M_sun]")

    plt.xlabel("Instantaneous SFR [M_sun / yr]")
    plt.ylabel("[OIII] flux [erg / s / cm^2]")
    plt.title("[OIII] flux vs SFR, coloured by stellar mass")

    if logx:
        plt.xscale("log")
    if logy:
        plt.yscale("log")

    plt.grid(alpha=0.3)
    plt.tight_layout()
    outname = filename.replace(".fits", "_oiii_vs_sfr.pdf")
    plt.savefig(outname, dpi=150)
    #plt.show()
    print(f"Saved plot to {outname}")



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Generate BAGPIPES grid and optionally plot results.")
    parser.add_argument("--plot", action="store_true", help="If set, plot results after generation.")
    parser.add_argument("--type", choices=["tau_age", "halpha_sfr"], default="tau_age",
                        help="Type of plot: 'tau_age' (default) or 'halpha_sfr'")
    parser.add_argument("--xcol", default="tau_gyr", help="Column for x-axis when plotting (tau_age only).")
    parser.add_argument("--ycol", default="age_gyr", help="Column for y-axis when plotting (tau_age only).")
    parser.add_argument("--logx", action="store_true", help="Use log scale for x-axis.")
    parser.add_argument("--logy", action="store_true", help="Use log scale for y-axis.")
    parser.add_argument("--logc", action="store_true", help="Use log scale for colour (halpha_sfr only).")
    args = parser.parse_args()

    fits_out = main()

    if args.plot:
        if args.type == "tau_age":
            plot_halpha_2d(fits_out, args.xcol, args.ycol, logx=args.logx, logy=args.logy)
        elif args.type == "halpha_sfr":
            plot_halpha_vs_sfr(fits_out, logx=args.logx, logy=args.logy, logc=args.logc)

