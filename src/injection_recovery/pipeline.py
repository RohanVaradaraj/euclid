"""
Combines all components into one high-level workflow.

Created: Wednesday 4th December 2024.
"""

from luminosity_function import LuminosityFunction
from source_injector import SourceInjector

#! Read in config file
config = load_config("config.yaml")

#! Draw sample from luminosity function
luminosity_function_params = config['luminosity_function']
luminosity_function = LuminosityFunction(luminosity_function_params)
sample = luminosity_function.sample_luminosities()

#! Inject sources into image