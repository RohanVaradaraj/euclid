"""
Contains the LuminosityFunction class for drawing objects from this M distribution.

Created: Wednesday 4th December 2024.
"""

from utils import load_config
import numpy as np
import matplotlib.pyplot as plt

class LuminosityFunction:
    def __init__(self, params):
        self.params = params
        self.alpha = params['alpha']
        self.beta = params['beta']
        self.M_star = params['M_star']
        self.phi_star = eval(params['phi_star'])
        self.Muv_range = params['Muv_range']
        self.n_samples = params['n_samples']


    def phi(self, M):
        """
        Compute the value of the luminosity function at M.
        """

        const = np.log(10) * self.phi_star / 2.5
        alpha_cpt = 10 ** ( 0.4*(self.alpha+1) * (M-self.M_star))
        beta_cpt = 10 ** ( 0.4*(self.beta+1) * (M-self.M_star))

        denominator = alpha_cpt + beta_cpt
        
        return const / denominator



    def distribution(self, num_points=1000):
        """
        Compute the luminosity function over the magnitude range, for plotting.
        """
        M_values = np.linspace(self.Muv_range[0], self.Muv_range[1], num_points)
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
            [self.phi(M) for M in np.linspace(self.Muv_range[0], self.Muv_range[1], 1000)],
            np.linspace(self.Muv_range[0], self.Muv_range[1], 1000))



    def sample_luminosities(self, n_samples=self.n_samples):
        """
        Sample mags from the double power-law distribution.
        """
        # Generate the CDF
        M_values = np.linspace(self.Muv_range[0], self.Muv_range[1], 1000)
        phi_values = [self.phi(M) for M in M_values]
        cdf = np.cumsum(phi_values) / np.sum(phi_values)
        
        # Invert the CDF to draw samples
        uniform_randoms = np.random.uniform(0, 1, n_samples)
        samples = np.interp(uniform_randoms, cdf, M_values)
        
        return samples


if __name__ == "__main__":

    config = load_config("config.yaml")
    
    luminosity_function_params = config['luminosity_function']
    luminosity_function = LuminosityFunction(luminosity_function_params)

    print(luminosity_function.phi(-23))

    sample = luminosity_function.sample_luminosities(1000000)
    print(sample)
    plt.hist(sample, bins=1000, density=True)

    # Plot normalised lf dist to compare with hist
    M_values, phi_values = luminosity_function.distribution()
    plt.plot(M_values, phi_values/np.max(phi_values), label='Luminosity function')
    plt.yscale('log')
    plt.show()
