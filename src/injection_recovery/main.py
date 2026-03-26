#!/usr/bin/env python3

"""
Runs injection-recovery pipeline.

Created: Wednesday 4th December 2024.
"""

from pipeline import RunFullInjectionRecoveryPipeline

def main():
    RunFullInjectionRecoveryPipeline(base_image='HSC-Z_DR3.fits', overwrite=True)

if __name__ == "__main__":
    main()