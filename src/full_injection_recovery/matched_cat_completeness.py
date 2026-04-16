#!/usr/bin/env python3

"""
Stack matched catalogues and plot completeness heatmaps for:
1. LBGs as a function of Muv and z
2. M dwarfs as a function of HSC-Z apparent magnitude and subtype
3. L dwarfs as a function of HSC-Z apparent magnitude and subtype
4. T dwarfs as a function of HSC-Z apparent magnitude and subtype

By default, completeness is defined using the detection-band recovered flag:
`recovered_HSC_Z_DR3`.

Run like:
python3 matched_cat_completeness.py  \ 
--matched-dir catalogues/matched  \ 
--write-stacked catalogues/XMM3/all_matched_catalogues_XMM3.fits  \ 
--output-plot ../../plots/completeness/matched_catalogue_completeness_heatmaps_XMM3.pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack
from tqdm import tqdm

from utils import load_config


plt.rcParams.update({"font.size": 14})
plt.rcParams["axes.linewidth"] = 2.5
plt.rcParams["figure.dpi"] = 120


def counts_to_abmag(counts, zeropoint):
    counts = np.asarray(counts, dtype=float)
    mag = np.full(counts.shape, np.nan, dtype=float)
    good = counts > 0
    mag[good] = -2.5 * np.log10(counts[good]) + zeropoint
    return mag


def get_filter_zeropoint(config: dict, filter_name: str) -> float:
    injection_cfg = config["source_injection"]
    filters = list(injection_cfg["filters"])
    zeropoints = list(injection_cfg["zeropoints"])
    if filter_name not in filters:
        raise KeyError(f"Filter '{filter_name}' not found in source_injection.filters")
    return float(zeropoints[filters.index(filter_name)])


def as_boolean_array(values):
    arr = np.asarray(values)
    if arr.dtype.kind == "b":
        return arr.astype(bool)

    text = np.char.strip(arr.astype(str))
    return np.isin(text, ["True", "true", "T", "t", "1", "yes", "Y"])


def load_and_stack_matched_catalogues(matched_dir: Path) -> Table:
    files = sorted(matched_dir.glob("*.fits"))
    if not files:
        raise FileNotFoundError(f"No matched catalogues found in {matched_dir}")

    tables = [
        Table.read(fp)
        for fp in tqdm(files, desc="Reading matched catalogues", unit="file")
    ]
    return vstack(tables, metadata_conflicts="silent")


def extract_subtype_number(type_values):
    subtype = np.full(len(type_values), np.nan, dtype=float)
    for i, value in enumerate(type_values):
        text = str(value)
        if len(text) >= 2 and text[0] in {"M", "L", "T"}:
            try:
                subtype[i] = float(text[1:])
            except ValueError:
                pass
    return subtype


def is_stellar_subtype_label(text, family):
    text = str(text)
    if not text.startswith(family) or len(text) < 2:
        return False
    try:
        float(text[1:])
    except ValueError:
        return False
    return True


def sort_subtype_labels(type_values, family):
    labels = [str(value) for value in type_values if is_stellar_subtype_label(value, family)]
    unique_labels = sorted(set(labels), key=lambda label: float(label[1:]))
    return unique_labels


def compute_fraction_grid(x, y, recovered, x_bins, y_bins):
    total, _, _ = np.histogram2d(x, y, bins=[x_bins, y_bins])
    found, _, _ = np.histogram2d(x[recovered], y[recovered], bins=[x_bins, y_bins])

    with np.errstate(divide="ignore", invalid="ignore"):
        frac = found / total
    frac[~np.isfinite(frac)] = np.nan
    return frac, total


def build_mag_bins(values, default_bins):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return default_bins

    if np.any((values >= default_bins[0]) & (values <= default_bins[-1])):
        return default_bins

    step = default_bins[1] - default_bins[0]
    min_edge = step * np.floor(np.nanmin(values) / step)
    max_edge = step * np.ceil(np.nanmax(values) / step)
    if max_edge <= min_edge:
        max_edge = min_edge + step
    return np.arange(min_edge, max_edge + step, step)


def plot_heatmap(ax, grid, x_bins, y_bins, xlabel, ylabel, title, y_is_category=False):
    x_edges = np.asarray(x_bins, dtype=float)
    y_edges = np.asarray(y_bins, dtype=float)

    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        grid.T,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        shading="auto",
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if not y_is_category:
        ax.invert_yaxis()

    return mesh


def make_lbg_panel(ax, table, recovered_col, muv_bins, z_bins):
    lbg = np.asarray(table["type"].astype(str)) == "LBG"
    if not np.any(lbg):
        ax.set_visible(False)
        return None

    x = np.asarray(table["z"][lbg], dtype=float)
    y = np.asarray(table["Muv"][lbg], dtype=float)
    recovered = as_boolean_array(table[recovered_col][lbg])

    grid, _ = compute_fraction_grid(x, y, recovered, z_bins, muv_bins)
    return plot_heatmap(
        ax,
        grid,
        z_bins,
        muv_bins,
        xlabel=r"$z$",
        ylabel=r"$M_{\rm UV}$",
        title="LBG Completeness",
    )


def make_dwarf_panel(ax, table, recovered_col, family, mag_bins, hsc_z_zeropoint):
    types = np.asarray(table["type"].astype(str))
    family_mask = np.array([is_stellar_subtype_label(value, family) for value in types], dtype=bool)
    if not np.any(family_mask):
        ax.set_visible(False)
        return None

    family_types = types[family_mask]
    subtype_labels = sort_subtype_labels(family_types, family)
    if not subtype_labels:
        ax.set_visible(False)
        return None

    subtype_map = {label: idx for idx, label in enumerate(subtype_labels)}
    subtype_index = np.array([subtype_map.get(label, -1) for label in family_types], dtype=int)
    valid_subtype = subtype_index >= 0

    z_flux_col = "model_flux_HSC-Z_DR3"
    z_mag = counts_to_abmag(table[z_flux_col][family_mask], hsc_z_zeropoint)

    recovered = as_boolean_array(table[recovered_col][family_mask])
    valid = valid_subtype & np.isfinite(z_mag)

    if not np.any(valid):
        ax.set_visible(False)
        return None

    x = subtype_index[valid].astype(float)
    y = z_mag[valid]
    recovered = recovered[valid]

    mag_bins = build_mag_bins(y, mag_bins)
    subtype_bins = np.arange(len(subtype_labels) + 1, dtype=float) - 0.5
    grid, total = compute_fraction_grid(x, y, recovered, subtype_bins, mag_bins)
    if np.all(~np.isfinite(grid)) and np.any(total > 0):
        grid = np.where(total.T > 0, 0.0, np.nan)
    elif np.any(~np.isfinite(grid)):
        grid = np.where(np.isfinite(grid), grid, 0.0)

    mesh = plot_heatmap(
        ax,
        grid,
        subtype_bins,
        mag_bins,
        xlabel=f"{family} subtype",
        ylabel=r"$m_{\rm AB}(\mathrm{HSC\!-\!Z})$",
        title=f"{family} Dwarf Completeness",
    )

    subtype_ticks = np.arange(len(subtype_labels), dtype=float) + 0.5
    ax.set_xticks(subtype_ticks)
    ax.set_xticklabels(subtype_labels, rotation=45, ha="right")

    return mesh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path.cwd() / "config.yaml",
        help="Config file used to read filter zeropoints.",
    )
    parser.add_argument(
        "--matched-dir",
        type=Path,
        default=Path.cwd() / "catalogues" / "matched",
        help="Directory containing matched catalogue FITS files.",
    )
    parser.add_argument(
        "--recovered-col",
        default="recovered_HSC_Z_DR3",
        help="Recovered flag column to use for completeness.",
    )
    parser.add_argument(
        "--write-stacked",
        type=Path,
        default=None,
        help="Optional output FITS path for the vertically stacked catalogue.",
    )
    parser.add_argument(
        "--output-plot",
        type=Path,
        default=Path.cwd() / "matched_catalogue_completeness_heatmaps.pdf",
        help="Output path for the completeness heatmap figure.",
    )
    parser.add_argument(
        "--muv-min",
        type=float,
        default=-24.0,
    )
    parser.add_argument(
        "--muv-max",
        type=float,
        default=-19.0,
    )
    parser.add_argument(
        "--muv-step",
        type=float,
        default=0.2,
    )
    parser.add_argument(
        "--z-min",
        type=float,
        default=5.0,
    )
    parser.add_argument(
        "--z-max",
        type=float,
        default=7.0,
    )
    parser.add_argument(
        "--z-step",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--mag-min",
        type=float,
        default=22.0,
    )
    parser.add_argument(
        "--mag-max",
        type=float,
        default=27.5,
    )
    parser.add_argument(
        "--mag-step",
        type=float,
        default=0.25,
    )
    args = parser.parse_args()
    config = load_config(str(args.config))
    hsc_z_zeropoint = get_filter_zeropoint(config, "HSC-Z_DR3")

    if args.write_stacked is not None and args.write_stacked.exists():
        table = Table.read(args.write_stacked)
    else:
        table = load_and_stack_matched_catalogues(args.matched_dir)
        if args.write_stacked is not None:
            args.write_stacked.parent.mkdir(parents=True, exist_ok=True)
            table.write(args.write_stacked, overwrite=True)

    if args.recovered_col not in table.colnames:
        raise KeyError(f"Column '{args.recovered_col}' not found in stacked table.")
    if "type" not in table.colnames:
        raise KeyError("Column 'type' not found in stacked table.")
    if "model_flux_HSC-Z_DR3" not in table.colnames:
        raise KeyError("Column 'model_flux_HSC-Z_DR3' not found in stacked table.")

    muv_bins = np.arange(args.muv_min, args.muv_max + args.muv_step, args.muv_step)
    z_bins = np.arange(args.z_min, args.z_max + args.z_step, args.z_step)
    mag_bins = np.arange(args.mag_min, args.mag_max + args.mag_step, args.mag_step)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=True)

    meshes = []
    mesh = make_lbg_panel(axes[0, 0], table, args.recovered_col, muv_bins, z_bins)
    if mesh is not None:
        meshes.append(mesh)

    for ax, family in zip([axes[0, 1], axes[1, 0], axes[1, 1]], ["M", "L", "T"]):
        mesh = make_dwarf_panel(ax, table, args.recovered_col, family, mag_bins, hsc_z_zeropoint)
        if mesh is not None:
            meshes.append(mesh)

    if meshes:
        cbar = fig.colorbar(meshes[0], ax=axes.ravel().tolist(), shrink=0.95)
        cbar.set_label("Completeness")

    args.output_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_plot, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
