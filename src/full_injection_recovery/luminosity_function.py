#!/usr/bin/env python3

"""
Contains the LuminosityFunction class for drawing objects from this M distribution.

Created: Wednesday 4th December 2024.
"""

import numpy as np
import matplotlib.pyplot as plt
from utils import load_config

class LuminosityFunction:
    def __init__(self, params):
        self.params = params
        self.alpha = params['alpha']
        self.beta = params['beta']
        self.M_star = params['M_star']
        self.phi_star = eval(params['phi_star'])
        self.Muv_range = params['Muv_range']
        self.n_samples = params['n_samples']
        self.dMuv = params['dMuv']
        self.type = params['type']


    def phi(self, M):
        """
        Compute luminosity function value at M (DPL or Schechter).
        """

        if self.type.lower() == 'dpl':
            const = np.log(10) * self.phi_star / 2.5
            alpha_cpt = 10 ** (0.4 * (self.alpha + 1) * (M - self.M_star))
            beta_cpt = 10 ** (0.4 * (self.beta + 1) * (M - self.M_star))
            denominator = alpha_cpt + beta_cpt
            return const / denominator

        if self.type.lower() == 'schechter':
            coeff = np.log(10) / 2.5
            faint = (10 ** (0.4 * (self.M_star - M))) ** (self.alpha + 1)
            bright_exponent = -10 ** (0.4 * (self.M_star - M))
            bright = np.exp(bright_exponent)
            return coeff * self.phi_star * faint * bright

        raise ValueError(f"Unsupported luminosity function type: {self.type}")



    def distribution(self, num_points=1000):
        """
        Compute the luminosity function over the magnitude range, for plotting.
        """
        M_values = np.arange(self.Muv_range[0], self.Muv_range[1], self.dMuv)
        phi_values = [self.phi(M) for M in M_values]
        return M_values, phi_values



    def pdf(self, M):
        """
        Convert the luminosity function to a probability density function.
        """
        return self.phi(M) / self._integrate_phi()



    def _integrate_phi(self):
        """
        Numerically integrate the luminosity function over the magnitude range.
        """
        return np.trapz(
            [self.phi(M) for M in np.arange(self.Muv_range[0], self.Muv_range[1], self.dMuv)],
            np.arange(self.Muv_range[0], self.Muv_range[1], self.dMuv))



    def sample_luminosities(self, uniform=True):
        """
        Sample Muv from the luminosity function using inverse transform sampling.

        Parameters:
            uniform (bool): Whether to sample from a uniform distribution.

        Returns:
            samples_Muv (array): The sampled Muv values.
        """
        M_values = np.arange(self.Muv_range[0], self.Muv_range[1], self.dMuv)
        phi_values = np.array([self.phi(M) for M in M_values])

        if uniform:
            samples_Muv = np.random.uniform(self.Muv_range[0], self.Muv_range[1], self.n_samples)

            # Snap samples_Muv to the Muv grid
            samples_Muv = np.round(samples_Muv / self.dMuv) * self.dMuv
        
        else:
            # Normalize the phi values to create a cumulative distribution
            cdf = np.cumsum(phi_values) / np.sum(phi_values)

            # Use vectorized interpolation to generate samples
            uniform_randoms = np.random.uniform(0, 1, self.n_samples)
            samples_Muv = np.interp(uniform_randoms, cdf, M_values)

            # Snap samples_Muv to the Muv grid
            samples_Muv = np.round(samples_Muv / self.dMuv) * self.dMuv

        return samples_Muv



if __name__ == "__main__":
    # Load configuration from YAML or other source
    config = load_config("config.yaml")
    
    luminosity_function_params = config['luminosity_function']
    luminosity_function = LuminosityFunction(luminosity_function_params)

    # Sample Muv values
    M_samples = luminosity_function.sample_luminosities(uniform=True)

    # Plot histogram of M_samples
    plt.hist(M_samples, bins=100, alpha=0.6, color='g', label='Muv Samples', density=False)

    # Plot normalised luminosity function for Muv to compare with the histogram
    M_values, phi_values = luminosity_function.distribution()
    plt.plot(M_values, np.array(phi_values)/np.max(phi_values), label='Luminosity Function')
    plt.yscale('log')
    plt.xlabel("Magnitude (Muv)")
    plt.ylabel("Probability Density")
    plt.legend()
    plt.show()

    # Print some sample Muv values for verification
    print(f"Sample Muv values: {M_samples[:10]}")
