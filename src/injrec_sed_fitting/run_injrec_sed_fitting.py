#!/usr/bin/env python3
"""
Run automated LePhare SED fitting on stacked injection-recovery matched catalogues.

The pipeline:
1. reads a stacked matched catalogue for one subfield,
2. records 5-sigma detection / 2-sigma non-detection pass flags,
3. runs LePhare for an LBG template set and a brown-dwarf template set on the
   selected rows only,
4. parses the outputs, computes an Muv from the best-fit LBG SED,
5. writes a merged FITS catalogue with the selection flags and fit summaries.
"""

from __future__ import annotations

import argparse
from astropy.cosmology import FlatLambdaCDM
from astropy.io import ascii
from astropy.table import Column, Table
import numpy as np
from pathlib import Path
import shlex
import shutil
import stat
import subprocess
import yaml

from sed_fitting_codes import (
    GenerateLePhareConfig,
    buildLePhareLibrary,
    parse_spec_file,
    runPhotometricRedshifts,
)


COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
PACKAGE_ROOT = Path(__file__).resolve().parent
LEPHARE_INPUT_DIR = Path.home().parents[1] / "hoy" / "temporaryFilesROHAN" / "lephare" / "inputs" / "euclid"
LEPHARE_CONFIG_DIR = Path.home() / "lephare" / "lephare_dev" / "config"
LEPHARE_TEST_DIR = Path.home() / "lephare" / "lephare_dev" / "test"
SHELL_SCRIPT_DIR = PACKAGE_ROOT / "shell_scripts"

STELLAR_TYPE_MAP = {
    1: "M4",
    2: "M5",
    3: "M6",
    4: "M7",
    5: "M8",
    6: "M9",
    7: "L0",
    8: "L1",
    9: "L2",
    10: "L3",
    11: "L4",
    12: "L5",
    13: "L6",
    14: "L7",
    15: "L8",
    16: "L9",
    17: "T0",
    18: "T1",
    19: "T2",
    20: "T3",
    21: "T4",
    22: "T5",
    23: "T6",
    24: "T7",
    25: "T8",
}


def load_yaml(path: Path):
    with open(path, "r") as handle:
        return yaml.safe_load(handle)


def zeropoint_flux_scale(zeropoint: float) -> float:
    return 10 ** (-0.4 * (float(zeropoint) + 48.6))


def zeropoint_map_from_config(config: dict) -> dict[str, float]:
    injection_cfg = config["source_injection"]
    filters = list(injection_cfg["filters"])
    zeropoints = list(injection_cfg["zeropoints"])
    if len(filters) != len(zeropoints):
        raise ValueError("source_injection.filters and source_injection.zeropoints must have the same length.")
    return {str(filter_name): float(zeropoint) for filter_name, zeropoint in zip(filters, zeropoints)}


def normalize_field_name(text: str) -> str:
    text = str(text)
    if text.startswith("XMM"):
        return "XMM"
    if text.startswith("CDFS"):
        return "CDFS"
    if text.startswith("COSMOS"):
        return "COSMOS"
    return text


def infer_selection_filters(field_name: str):
    base = normalize_field_name(field_name)
    if base == "CDFS":
        return "HSC-Z", ["HSC-G", "r"]
    return "HSC-Z_DR3", ["HSC-G_DR3", "HSC-R_DR3"]


def safe_suffix(text: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in str(text))


def flux_column(filter_name: str) -> str:
    return f"FLUX_APER_{safe_suffix(filter_name)}"


def fluxerr_column(filter_name: str) -> str:
    return f"FLUXERR_APER_{safe_suffix(filter_name)}"


def ensure_required_columns(table: Table, filters: list[str]):
    missing = []
    for filter_name in filters:
        for col in [flux_column(filter_name), fluxerr_column(filter_name)]:
            if col not in table.colnames:
                missing.append(col)
    if missing:
        raise KeyError(f"Missing required photometry columns: {missing}")


def compute_snr_flags(table: Table, detection_filter: str, non_detection_filters: list[str], detection_sigma: float, non_detection_sigma: float):
    flags = {}

    det_flux = np.asarray(table[flux_column(detection_filter)], dtype=float)
    det_err = np.asarray(table[fluxerr_column(detection_filter)], dtype=float)
    det_valid = np.isfinite(det_flux) & np.isfinite(det_err) & (det_err > 0)
    det_snr = np.full(len(table), np.nan, dtype=float)
    det_snr[det_valid] = det_flux[det_valid] / det_err[det_valid]
    flags[f"pass_{int(detection_sigma)}sig_det_{safe_suffix(detection_filter)}"] = det_valid & (det_snr >= detection_sigma)

    for filter_name in non_detection_filters:
        flux = np.asarray(table[flux_column(filter_name)], dtype=float)
        err = np.asarray(table[fluxerr_column(filter_name)], dtype=float)
        valid = np.isfinite(flux) & np.isfinite(err) & (err > 0)
        snr = np.full(len(table), np.nan, dtype=float)
        snr[valid] = flux[valid] / err[valid]
        flags[f"pass_{int(non_detection_sigma)}sig_nondet_{safe_suffix(filter_name)}"] = valid & (snr < non_detection_sigma)

    selection_mask = np.ones(len(table), dtype=bool)
    for values in flags.values():
        selection_mask &= values
    flags["pass_sed_selection"] = selection_mask

    return flags


def add_or_replace_column(table: Table, name: str, values):
    if name in table.colnames:
        table[name] = values
    else:
        table.add_column(Column(values, name=name))


def add_measured_flux_cgs_columns(table: Table, filters: list[str], zeropoints: dict[str, float]):
    for filter_name in filters:
        if filter_name not in zeropoints:
            raise KeyError(f"Missing zeropoint for filter {filter_name}")

        suffix = safe_suffix(filter_name)
        flux_scale = zeropoint_flux_scale(zeropoints[filter_name])
        measured_flux = np.asarray(table[flux_column(filter_name)], dtype=float) * flux_scale
        measured_fluxerr = np.asarray(table[fluxerr_column(filter_name)], dtype=float) * flux_scale

        add_or_replace_column(table, f"measured_flux_cgs_{suffix}", measured_flux)
        add_or_replace_column(table, f"measured_fluxerr_cgs_{suffix}", measured_fluxerr)

    return table


def prepare_catalogue(table: Table, detection_filter: str, non_detection_filters: list[str], detection_sigma: float, non_detection_sigma: float):
    working = table.copy(copy_data=True)
    if "injrec_sed_id" not in working.colnames:
        working.add_column(Column(np.arange(1, len(working) + 1, dtype=int), name="injrec_sed_id"), index=0)

    flags = compute_snr_flags(
        working,
        detection_filter=detection_filter,
        non_detection_filters=non_detection_filters,
        detection_sigma=detection_sigma,
        non_detection_sigma=non_detection_sigma,
    )
    for name, values in flags.items():
        add_or_replace_column(working, name, values.astype(bool))

    return working


def context_value(filters: list[str]) -> int:
    total = 0
    for i in range(len(filters)):
        total += 2**i
    return total


def write_lephare_input(table: Table, filters: list[str], zeropoints: dict[str, float], output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out = Table()
    out["ID"] = np.asarray(table["injrec_sed_id"], dtype=int)
    for filter_name in filters:
        if filter_name not in zeropoints:
            raise KeyError(f"Missing zeropoint for filter {filter_name}")
        flux_scale = zeropoint_flux_scale(zeropoints[filter_name])
        out[f"flux_{filter_name}"] = np.asarray(table[flux_column(filter_name)], dtype=float) * flux_scale
        out[f"err_{filter_name}"] = np.asarray(table[fluxerr_column(filter_name)], dtype=float) * flux_scale
    out["Context"] = np.full(len(out), context_value(filters), dtype=int)
    ascii.write(out, output_path, format="commented_header", overwrite=True)
    print(
        f"Wrote LePhare input with {len(out)} rows and {len(filters)} filters to {output_path}"
    )
    print(f"LePhare filter order: {filters}")


def lbg_filters(all_filters: list[str]) -> list[str]:
    return list(all_filters)


def bd_filters(all_filters: list[str], field_name: str) -> list[str]:
    filters = list(all_filters)
    for remove in ["HSC-G_DR3", "HSC-R_DR3", "f277w", "f444w", "ch1cds", "ch2cds"]:
        if remove in filters:
            filters.remove(remove)
    if normalize_field_name(field_name) == "CDFS":
        for remove in ["HSC-G", "HSC-R", "u", "g", "r"]:
            if remove in filters:
                filters.remove(remove)
    return filters


def lephare_out_path(input_name: str) -> Path:
    return LEPHARE_TEST_DIR / input_name.replace(".in", ".out").replace("_DR3", "")


def spec_path(zphot_dir: Path, obj_id: int) -> Path:
    return zphot_dir / f"Id{int(obj_id):09d}.spec"


def compute_muv_from_spec(spec_file: Path) -> float:
    parsed = parse_spec_file(spec_file)
    model = parsed["model"]
    sed_sections = parsed["sed"]
    if len(model) == 0 or len(sed_sections) == 0:
        return np.nan

    zphot = float(model["Zphot"][0])
    primary_sed = sed_sections[0]
    wave_obs = np.asarray(primary_sed["lambda"], dtype=float)
    mag_obs = np.asarray(primary_sed["flux"], dtype=float)
    target_wave = 1500.0 * (1.0 + zphot)
    if target_wave < np.nanmin(wave_obs) or target_wave > np.nanmax(wave_obs):
        return np.nan

    m_obs = float(np.interp(target_wave, wave_obs, mag_obs))
    dl_pc = COSMO.luminosity_distance(zphot).value * 1e6
    return m_obs - 5.0 * np.log10(dl_pc / 10.0) + 2.5 * np.log10(1.0 + zphot)


def stellar_type_from_model(model_number) -> str:
    try:
        model_number = int(model_number)
    except Exception:
        return "UNKNOWN"
    return STELLAR_TYPE_MAP.get(model_number, "UNKNOWN")


def parse_lbg_results(out_path: Path, zphot_dir: Path):
    results = {}
    if not out_path.is_file():
        return results

    table = ascii.read(out_path, format="no_header")
    for row in table:
        obj_id = int(row["col1"])
        spec_file = spec_path(zphot_dir, obj_id)
        muv = compute_muv_from_spec(spec_file) if spec_file.is_file() else np.nan
        results[obj_id] = {
            "sedfit_lbg_z": float(row["col2"]),
            "sedfit_lbg_chi2": float(row["col6"]),
            "sedfit_lbg_Muv": float(muv) if np.isfinite(muv) else np.nan,
        }
    return results


def parse_bd_results(out_path: Path, zphot_dir: Path):
    results = {}
    if not out_path.is_file():
        return results

    table = ascii.read(out_path, format="no_header")
    chi2_col = "col21" if "col21" in table.colnames else "col6"
    for row in table:
        obj_id = int(row["col1"])
        model_number = np.nan
        spec_file = spec_path(zphot_dir, obj_id)
        if spec_file.is_file():
            parsed = parse_spec_file(spec_file)
            model_table = parsed.get("model")
            if model_table is not None and len(model_table) > 0:
                model_number = model_table["Model"][-1]
        results[obj_id] = {
            "sedfit_bd_model": int(model_number) if np.isfinite(model_number) else -1,
            "sedfit_bd_type": stellar_type_from_model(model_number),
            "sedfit_bd_chi2": float(row[chi2_col]),
        }
    return results


def merge_results(table: Table, selected_mask, lbg_results, bd_results):
    n = len(table)
    add_or_replace_column(table, "sedfit_ran", np.asarray(selected_mask, dtype=bool))
    add_or_replace_column(table, "sedfit_lbg_z", np.full(n, np.nan, dtype=float))
    add_or_replace_column(table, "sedfit_lbg_chi2", np.full(n, np.nan, dtype=float))
    add_or_replace_column(table, "sedfit_lbg_Muv", np.full(n, np.nan, dtype=float))
    add_or_replace_column(table, "sedfit_bd_model", np.full(n, -1, dtype=int))
    add_or_replace_column(table, "sedfit_bd_type", np.full(n, "NONE", dtype="U16"))
    add_or_replace_column(table, "sedfit_bd_chi2", np.full(n, np.nan, dtype=float))

    id_to_index = {int(obj_id): idx for idx, obj_id in enumerate(np.asarray(table["injrec_sed_id"], dtype=int))}

    for obj_id, values in lbg_results.items():
        idx = id_to_index.get(int(obj_id))
        if idx is None:
            continue
        table["sedfit_lbg_z"][idx] = values["sedfit_lbg_z"]
        table["sedfit_lbg_chi2"][idx] = values["sedfit_lbg_chi2"]
        table["sedfit_lbg_Muv"][idx] = values["sedfit_lbg_Muv"]

    for obj_id, values in bd_results.items():
        idx = id_to_index.get(int(obj_id))
        if idx is None:
            continue
        table["sedfit_bd_model"][idx] = values["sedfit_bd_model"]
        table["sedfit_bd_type"][idx] = values["sedfit_bd_type"]
        table["sedfit_bd_chi2"][idx] = values["sedfit_bd_chi2"]

    return table


def default_full_injrec_config() -> Path:
    return PACKAGE_ROOT.parent / "full_injection_recovery" / "config.yaml"


def run_mode(
    mode_name: str,
    selected_table: Table,
    filters: list[str],
    zeropoints: dict[str, float],
    det_filters: list[str],
    field_name: str,
    output_root: Path,
):
    mode_root = output_root / mode_name
    mode_root.mkdir(parents=True, exist_ok=True)
    input_name = f"{output_root.name}_{mode_name}.in"
    input_path = LEPHARE_INPUT_DIR / input_name
    print(f"\n[{mode_name}] Preparing LePhare input")
    write_lephare_input(selected_table, filters, zeropoints, input_path)

    para_name = f"{output_root.name}_{mode_name}.para"
    para_path = LEPHARE_CONFIG_DIR / para_name
    print(f"[{mode_name}] Writing LePhare parameter file to {para_path}")
    GenerateLePhareConfig(
        field_name=field_name,
        all_filters=filters,
        det_filters=det_filters,
        run_type="",
        run_brown_dwarfs=(mode_name == "bd"),
        run_dusty=False,
        run_lya=False,
        file_name=para_path,
        spec_out=True,
        custom_name=input_name,
    )

    print(f"[{mode_name}] Building LePhare libraries/filters/mags")
    buildLePhareLibrary(para_name, build_libs=True, build_filters=True, build_mags=True)

    zphot_dir = mode_root / "zphot"
    zphot_dir.mkdir(parents=True, exist_ok=True)
    print(f"[{mode_name}] Running zphota; .spec files will go to {zphot_dir}")
    runPhotometricRedshifts(para_name, zphot_dir=zphot_dir)

    out_path = lephare_out_path(input_name)
    copied_out = mode_root / out_path.name
    if out_path.is_file():
        copied_out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(out_path, copied_out)
        print(f"[{mode_name}] Copied LePhare .out file to {copied_out}")
    else:
        print(f"[{mode_name}] Warning: expected LePhare .out file not found at {out_path}")

    return copied_out, zphot_dir


def default_output_root(stacked_catalogue: Path) -> Path:
    subfield_name = stacked_catalogue.parent.name
    return PACKAGE_ROOT / "catalogues" / subfield_name / stacked_catalogue.stem


def build_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stacked-catalogue", type=Path, required=True)
    parser.add_argument(
        "--full-injrec-config",
        type=Path,
        default=default_full_injrec_config(),
    )
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--detection-filter", default=None)
    parser.add_argument("--non-detection-filters", nargs="*", default=None)
    parser.add_argument("--detection-sigma", type=float, default=5.0)
    parser.add_argument("--non-detection-sigma", type=float, default=2.0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--submit-to-queue", action="store_true")
    parser.add_argument("--queue", default="cmb")
    parser.add_argument("--queue-memory-gb", type=int, default=8)
    parser.add_argument("--job-name", default=None)
    return parser


def _safe_job_name(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ["_", "-"] else "_" for ch in str(text))


def _build_local_command(args) -> list[str]:
    command = [
        "python3",
        str((PACKAGE_ROOT / "run_injrec_sed_fitting.py").resolve()),
        "--stacked-catalogue",
        str(args.stacked_catalogue.resolve()),
        "--full-injrec-config",
        str(args.full_injrec_config.resolve()),
    ]

    if args.output_root is not None:
        command.extend(["--output-root", str(args.output_root.resolve())])
    if args.detection_filter is not None:
        command.extend(["--detection-filter", str(args.detection_filter)])
    if args.non_detection_filters:
        command.append("--non-detection-filters")
        command.extend(str(item) for item in args.non_detection_filters)
    if float(args.detection_sigma) != 5.0:
        command.extend(["--detection-sigma", str(args.detection_sigma)])
    if float(args.non_detection_sigma) != 2.0:
        command.extend(["--non-detection-sigma", str(args.non_detection_sigma)])
    if args.overwrite:
        command.append("--overwrite")

    return command


def submit_to_queue(args):
    SHELL_SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    output_root = args.output_root or default_output_root(args.stacked_catalogue)
    inferred_name = args.job_name or f"injrec_sed_{args.stacked_catalogue.stem}"
    job_name = _safe_job_name(inferred_name)
    shell_path = SHELL_SCRIPT_DIR / f"{job_name}.sh"
    command = _build_local_command(args)
    command_str = " ".join(shlex.quote(part) for part in command)

    with open(shell_path, "w") as handle:
        handle.write("#!/bin/bash\n")
        handle.write(f"cd {shlex.quote(str(PACKAGE_ROOT))}\n")
        handle.write(command_str + "\n")

    shell_path.chmod(shell_path.stat().st_mode | stat.S_IXUSR)

    print(f"Queue shell script written to {shell_path}")
    print(f"Queue job output will be managed by addqueue")
    print(f"Queued run will write outputs under {output_root}")

    subprocess.run(
        [
            "addqueue",
            "-c",
            job_name,
            "-m",
            str(args.queue_memory_gb),
            "-q",
            str(args.queue),
            f"./{shell_path.relative_to(PACKAGE_ROOT)}",
        ],
        cwd=PACKAGE_ROOT,
        check=True,
    )
    print(
        f"Submitted queue job {job_name} to queue={args.queue} with memory={args.queue_memory_gb} GB"
    )


def main():
    parser = build_argument_parser()
    args = parser.parse_args()

    if args.submit_to_queue:
        submit_to_queue(args)
        return

    if not args.full_injrec_config.is_file():
        raise FileNotFoundError(
            f"Could not find full_injection_recovery config at {args.full_injrec_config}"
        )
    if not args.stacked_catalogue.is_file():
        raise FileNotFoundError(
            f"Could not find stacked matched catalogue at {args.stacked_catalogue}"
        )

    injrec_config = load_yaml(args.full_injrec_config)
    all_filters = list(injrec_config["source_injection"]["filters"])
    zeropoints = zeropoint_map_from_config(injrec_config)
    subfield_name = args.stacked_catalogue.parent.name
    field_name = normalize_field_name(subfield_name)
    print(f"Reading stacked matched catalogue: {args.stacked_catalogue}")
    print(f"Using full_injection_recovery config: {args.full_injrec_config}")
    print(f"Resolved subfield {subfield_name} -> field {field_name}")

    detection_filter, non_detection_filters = infer_selection_filters(subfield_name)
    if args.detection_filter is not None:
        detection_filter = args.detection_filter
    if args.non_detection_filters is not None and len(args.non_detection_filters) > 0:
        non_detection_filters = args.non_detection_filters
    print(f"5 sigma detection filter: {detection_filter}")
    print(f"2 sigma non-detection filters: {non_detection_filters}")

    required_filters = list(dict.fromkeys(all_filters + [detection_filter] + non_detection_filters))

    stacked = Table.read(args.stacked_catalogue)
    print(f"Loaded stacked catalogue with {len(stacked)} rows")
    ensure_required_columns(stacked, required_filters)
    prepared = prepare_catalogue(
        stacked,
        detection_filter=detection_filter,
        non_detection_filters=non_detection_filters,
        detection_sigma=args.detection_sigma,
        non_detection_sigma=args.non_detection_sigma,
    )
    add_measured_flux_cgs_columns(prepared, all_filters, zeropoints)
    det_col = f"pass_{int(args.detection_sigma)}sig_det_{safe_suffix(detection_filter)}"
    print(f"Sources passing {args.detection_sigma:g} sigma detection in {detection_filter}: {np.count_nonzero(prepared[det_col])}")
    for filter_name in non_detection_filters:
        nondet_col = f"pass_{int(args.non_detection_sigma)}sig_nondet_{safe_suffix(filter_name)}"
        print(
            f"Sources passing {args.non_detection_sigma:g} sigma non-detection in {filter_name}: "
            f"{np.count_nonzero(prepared[nondet_col])}"
        )
    print(f"Added measured cgs flux/error columns for filters: {all_filters}")

    output_root = args.output_root or default_output_root(args.stacked_catalogue)
    if output_root.exists() and args.overwrite:
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    print(f"Writing injrec SED-fitting outputs under {output_root}")

    selected_mask = np.asarray(prepared["pass_sed_selection"], dtype=bool)
    selected = prepared[selected_mask]
    print(f"Sources passing all SED-selection cuts: {len(selected)} / {len(prepared)}")

    selected_path = output_root / f"{args.stacked_catalogue.stem}_selection_flags.fits"
    prepared.write(selected_path, overwrite=True)
    print(f"Wrote selection-flag catalogue to {selected_path}")

    if len(selected) == 0:
        merged_path = output_root / f"{args.stacked_catalogue.stem}_with_sed_fits.fits"
        merge_results(prepared, selected_mask, {}, {})
        prepared.write(merged_path, overwrite=True)
        print("No sources passed the selection cuts; wrote flag-only merged catalogue.")
        print(f"Wrote merged catalogue to {merged_path}")
        return

    lbg_out, lbg_zphot_dir = run_mode(
        mode_name="lbg",
        selected_table=selected,
        filters=lbg_filters(all_filters),
        zeropoints=zeropoints,
        det_filters=[detection_filter],
        field_name=field_name,
        output_root=output_root,
    )
    bd_out, bd_zphot_dir = run_mode(
        mode_name="bd",
        selected_table=selected,
        filters=bd_filters(all_filters, field_name),
        zeropoints=zeropoints,
        det_filters=[detection_filter],
        field_name=field_name,
        output_root=output_root,
    )

    lbg_results = parse_lbg_results(lbg_out, lbg_zphot_dir)
    bd_results = parse_bd_results(bd_out, bd_zphot_dir)
    print(f"Parsed {len(lbg_results)} LBG fit rows from {lbg_out}")
    print(f"Parsed {len(bd_results)} BD fit rows from {bd_out}")
    merged = merge_results(prepared, selected_mask, lbg_results, bd_results)

    merged_path = output_root / f"{args.stacked_catalogue.stem}_with_sed_fits.fits"
    merged.write(merged_path, overwrite=True)
    print(f"Wrote merged SED-fitting catalogue to {merged_path}")


if __name__ == "__main__":
    main()
