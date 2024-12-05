"""
Runs injection-recovery pipeline.

Created: Wednesday 4th December 2024.
"""

from pipeline import RunFullInjectionRecoveryPipeline

def main():
    RunFullInjectionRecoveryPipeline(base_image='UVISTA_YJ_DR6.fits', overwrite=True)

if __name__ == "__main__":
    main()