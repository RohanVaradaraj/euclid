#!/usr/bin/env python3

"""Runs the full queue-driven injection-recovery pipeline."""

from pipeline import RunFullInjectionRecoveryPipeline


def main():
    RunFullInjectionRecoveryPipeline(overwrite=True)


if __name__ == "__main__":
    main()
