#!/usr/bin/env python3

"""Run the full injection-recovery pipeline sequentially for multiple fields."""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path
import tempfile

import yaml

from pipeline import RunFullInjectionRecoveryPipeline
from utils import load_config


PACKAGE_ROOT = Path(__file__).resolve().parent


def write_temp_field_config(base_config: dict, field_name: str) -> Path:
    config = deepcopy(base_config)
    config.setdefault("field", {})
    config["field"]["name"] = field_name

    tmp_file = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=f"_{field_name}.yaml",
        prefix="full_injrec_",
        dir=PACKAGE_ROOT,
        delete=False,
    )
    with tmp_file as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return Path(tmp_file.name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=PACKAGE_ROOT / "config.yaml",
        help="Base single-field config to clone and override per field.",
    )
    parser.add_argument(
        "--fields",
        nargs="+",
        required=True,
        help="Field names to run sequentially, e.g. XMM1 XMM2 XMM3 COSMOS.",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Do not pass overwrite=True into the single-field pipeline.",
    )
    args = parser.parse_args()

    base_config = load_config(args.config)
    overwrite = not args.no_overwrite

    for field_name in args.fields:
        print(f"\n================ Running field {field_name} ================\n")
        temp_config = write_temp_field_config(base_config, field_name)
        try:
            RunFullInjectionRecoveryPipeline(
                overwrite=overwrite,
                config_file=temp_config,
            )
        finally:
            temp_config.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
