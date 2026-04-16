#!/usr/bin/env python3
"""
Copied/adapted LePhare helper functions for injection-recovery SED fitting.

This file is intentionally independent of src/sed_fitting so the injrec
pipeline can evolve without modifying the science pipeline code.
"""

from __future__ import annotations

from astropy.table import Table
import os
from pathlib import Path
import re


def remove_items(master_list, items_to_remove):
    remaining = list(master_list)
    for item in items_to_remove:
        if item in remaining:
            remaining.remove(item)
    return remaining


def skip_initial_lines(lines):
    count = 0
    for idx, line in enumerate(lines):
        if not line.startswith("#"):
            count += 1
        if count == 3:
            return lines[idx + 1 :]
    return lines


def detect_section_changes(lines):
    line_lengths = [len(re.split(r"\s+", line.strip())) for line in lines]
    section_indices = []

    for i in range(1, len(line_lengths)):
        if line_lengths[i] != line_lengths[i - 1]:
            section_indices.append(i)

    for i in range(1, len(lines)):
        try:
            prev_val = float(re.split(r"\s+", lines[i - 1].strip())[0])
            curr_val = float(re.split(r"\s+", lines[i].strip())[0])
            if curr_val - prev_val > 500:
                section_indices.append(i)
                break
        except ValueError:
            continue

    return sorted(set(section_indices))


def split_sed_section(sed_lines):
    sed_sections = []
    current_sed = []

    for line in sed_lines:
        try:
            wavelength = float(re.split(r"\s+", line.strip())[0])
            if current_sed and wavelength < float(re.split(r"\s+", current_sed[-1].strip())[0]):
                sed_sections.append(current_sed)
                current_sed = []
            current_sed.append(line)
        except ValueError:
            continue

    if current_sed:
        sed_sections.append(current_sed)

    return sed_sections


def parse_spec_file(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    data_lines = skip_initial_lines(lines)
    data_lines = data_lines[1:]

    section_indices = detect_section_changes(data_lines)
    section_indices.append(len(data_lines))

    sections = []
    start_idx = 0
    for end_idx in section_indices:
        sections.append(data_lines[start_idx:end_idx])
        start_idx = end_idx

    tables = {}
    section_names = ["model", "phot", "zpdf", "sed"]

    for idx, section in enumerate(sections):
        section_name = section_names[idx] if idx < len(section_names) else f"unknown_{idx}"

        if section_name == "model":
            column_names = [
                "ID",
                "Nline",
                "Model",
                "Library",
                "Nband",
                "Zphot",
                "Zinf",
                "Zsup",
                "Chi2",
                "PDF",
                "Extlaw",
                "EB-V",
                "Lir",
                "Age",
                "Mass",
                "SFR",
                "SSFR",
            ]
        elif section_name == "phot":
            column_names = [f"col{i+1}" for i in range(len(section[0].split()))]
        elif section_name == "zpdf":
            column_names = ["z", "P(z)"]
        elif section_name == "sed":
            sed_sections = split_sed_section(section)
            sed_tables = []
            for sed_section in sed_sections:
                data_rows = [re.split(r"\s+", line.strip()) for line in sed_section if line.strip()]
                data_rows = [[float(val) for val in row] for row in data_rows]
                sed_tables.append(Table(rows=data_rows, names=["lambda", "flux"]))
            tables[section_name] = sed_tables
            continue
        else:
            column_names = [f"col{i+1}" for i in range(len(section[0].split()))]

        data_rows = [re.split(r"\s+", line.strip()) for line in section if line.strip()]
        not_allowed = ["GAL-", "QSO", "STAR"]
        data_rows = [
            [val if any(na in val for na in not_allowed) else float(val) for val in row]
            for row in data_rows
        ]
        tables[section_name] = Table(rows=data_rows, names=column_names)

    return tables


def buildLePhareLibrary(
    parameter_file: str,
    parameter_dir: Path = Path.home() / "lephare" / "lephare_dev" / "config",
    build_libs: bool = False,
    build_filters: bool = False,
    build_mags: bool = False,
) -> None:
    code_dir = Path(__file__).resolve().parent
    os.chdir(str(parameter_dir))

    if build_libs:
        os.system(f"$LEPHAREDIR/source/sedtolib -t S -c $LEPHAREDIR/config/{parameter_file}")
        os.system(f"$LEPHAREDIR/source/sedtolib -t Q -c $LEPHAREDIR/config/{parameter_file}")
        os.system(f"$LEPHAREDIR/source/sedtolib -t G -c $LEPHAREDIR/config/{parameter_file}")

    if build_filters:
        os.system(f"$LEPHAREDIR/source/filter -c $LEPHAREDIR/config/{parameter_file}")

    if build_mags:
        os.system(f"$LEPHAREDIR/source/mag_star -c $LEPHAREDIR/config/{parameter_file}")
        os.system(f"$LEPHAREDIR/source/mag_gal -t Q -c $LEPHAREDIR/config/{parameter_file}")
        os.system(f"$LEPHAREDIR/source/mag_gal -t G -c $LEPHAREDIR/config/{parameter_file}")

    os.chdir(code_dir)


def runPhotometricRedshifts(
    parameter_file: str,
    zphot_dir: Path | None,
    parameter_dir: Path = Path.home().parent.parent / "lephare" / "lephare_dev" / "config",
) -> None:
    code_dir = Path(__file__).resolve().parent

    if zphot_dir is not None:
        os.chdir(zphot_dir)
        os.system(f"$LEPHAREDIR/source/zphota -c $LEPHAREDIR/config/{parameter_file}")
        os.chdir(code_dir)
    else:
        os.system(f"$LEPHAREDIR/source/zphota -c $LEPHAREDIR/config/{parameter_file}")


def filter_files(CDFS=False):
    filt_files = {
        "u": "cfht/megacam/up.pb",
        "g": "cfht/megacam/gp.pb",
        "r": "cfht/megacam/rp.pb",
        "i": "cfht/megacam/ip.pb",
        "z": "cfht/megacam/zp.pb",
        "HSC-G_DR3": "myfilters/HSC/g_HSC.txt",
        "HSC-R_DR3": "myfilters/HSC/r_HSC.txt",
        "HSC-I_DR3": "myfilters/HSC/i_HSC.txt",
        "HSC-NB0816_DR3": "myfilters/HSC/nb816_HSC.txt",
        "HSC-Z_DR3": "myfilters/HSC/z_HSC.txt",
        "HSC-NB0921_DR3": "myfilters/HSC/nb921_HSC.txt",
        "HSC-Y_DR3": "myfilters/HSC/y_HSC.txt",
        "Y": "myfilters/VISTA/VISTA_Y.txt",
        "J": "myfilters/VISTA/VISTA_J.txt",
        "H": "myfilters/VISTA/VISTA_H.txt",
        "Ks": "myfilters/VISTA/VISTA_Ks.txt",
        "VIS": "myfilters/Euclid/Euclid_VIS.txt",
        "Ye": "myfilters/Euclid/Euclid_Y.txt",
        "Je": "myfilters/Euclid/Euclid_J.txt",
        "He": "myfilters/Euclid/Euclid_H.txt",
        "ch1cds": "myfilters/SPITZER/irac_ch1.txt",
        "ch2cds": "myfilters/SPITZER/irac_ch2.txt",
    }

    if CDFS:
        new_filt_files = {}
        for k, v in filt_files.items():
            if k.startswith("HSC-") and k.endswith("_DR3"):
                k = k.replace("_DR3", "")
            if k == "VIS":
                k = "VIS_Q1"
            elif k == "Ye":
                k = "YE_Q1"
            elif k == "Je":
                k = "JE_Q1"
            elif k == "He":
                k = "HE_Q1"
            new_filt_files[k] = v
        filt_files = new_filt_files

    return filt_files


def GenerateLePhareConfig(
    field_name,
    all_filters,
    det_filters,
    run_type,
    run_brown_dwarfs,
    run_dusty,
    run_lya,
    file_name=Path.home() / "lephare" / "lephare_dev" / "config" / "euclid.para",
    filter_file_name="FILTERS.filt",
    z_step=(0.05, 10.0, 0.05),
    spec_out=True,
    custom_name="",
) -> None:
    expanded_filename = os.path.expanduser(str(file_name))

    with open(expanded_filename, "w") as f:
        f.write("##############################################################################\n")
        f.write("#                CREATION OF LIBRARIES FROM SEDs List                        #\n")
        f.write("##############################################################################\n")
        f.write("#\n")
        f.write("STAR_SED\t$LEPHAREDIR/sed/STAR/DWARFSTARS/DWARFSTARS_MOD.list\n")
        f.write("STAR_FSCALE\t1.\n")
        f.write("STAR_LIB\tLIB_STAR\n")
        f.write("QSO_SED\t\t$LEPHAREDIR/sed/QSO/QSO_MOD.list\n")
        f.write("QSO_FSCALE\t1.\n")
        f.write("QSO_LIB\tLIB_QSO\n")
        if not run_lya:
            f.write("GAL_SED\t$LEPHAREDIR/sed/GAL/BC03/ASCII/tmp_BC03_MOD.lis\n")
        else:
            f.write("GAL_SED\t$LEPHAREDIR/sed/GAL/BC03/ASCII/BC03_MOD.lis\n")
        f.write("GAL_FSCALE\t1.\n")
        f.write("GAL_LIB\tLIB_GAL\n")
        f.write("SEL_AGE\t$LEPHAREDIR/sed/GAL/DEFAULT_BC03_CHAB/BC03_AGE.list\n")
        f.write("AGE_RANGE\t1.0e7,13.8e9\n")
        f.write("\n")

        filter_names = list(all_filters)
        if run_brown_dwarfs:
            filter_names = remove_items(
                filter_names,
                ["HSC-G_DR3", "HSC-R_DR3", "f277w", "f444w", "ch1cds", "ch2cds"],
            )
            if field_name == "CDFS":
                filter_names = remove_items(filter_names, ["HSC-G", "HSC-R", "u", "g", "r"])

        if not run_dusty:
            if field_name != "XMM":
                filter_names = remove_items(filter_names, ["f444w", "ch1cds", "ch2cds"])
            else:
                filter_names = remove_items(filter_names, ["f444w", "ch1servs", "ch2servs"])

        filter_dict = filter_files(CDFS=field_name == "CDFS")
        if field_name == "XMM":
            if "ch1cds" in filter_dict:
                filter_dict["ch1servs"] = filter_dict.pop("ch1cds")
            if "ch2cds" in filter_dict:
                filter_dict["ch2servs"] = filter_dict.pop("ch2cds")

        filter_list = "FILTER_LIST " + ",".join(filter_dict[name] for name in filter_names)
        f.write(filter_list + "\n")
        f.write("TRANS_TYPE\t0\n")
        f.write("FILTER_CALIB\t0\n")
        f.write(f"FILTER_FILE\t{filter_file_name}\n")
        f.write("\n")

        f.write("STAR_LIB_IN\tLIB_STAR\n")
        f.write("STAR_LIB_OUT\tSTAR_EUC\n")
        f.write("QSO_LIB_IN\tLIB_QSO\n")
        f.write("QSO_LIB_OUT\tQSO_EUC\n")
        f.write("GAL_LIB_IN\tLIB_GAL\n")
        f.write("GAL_LIB_OUT\tGAL_EUC\n")
        f.write("MAGTYPE\tAB\n")
        f.write(f"Z_STEP\t{z_step[0]},{z_step[1]},{z_step[2]}\n")
        f.write("COSMOLOGY\t70,0.3,0.7\n")
        f.write("MOD_EXTINC\t0,12\n")
        f.write("EXTINC_LAW\tcalzetti.dat\n")

        base_Av = [
            0.0,
            0.025,
            0.05,
            0.075,
            0.1,
            0.125,
            0.15,
            0.175,
            0.2,
            0.225,
            0.25,
            0.3,
            0.35,
            0.4,
            0.45,
            0.5,
            0.55,
            0.6,
            0.65,
            0.7,
            0.75,
            0.8,
            0.85,
            0.9,
            0.95,
            1.0,
        ]
        if run_dusty:
            base_Av.extend([1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.4, 3.6, 3.8, 4.0])
        f.write("EB_V\t" + ",".join(str(x) for x in base_Av) + "\n")
        f.write("EM_LINES\tNO\n")

        det_string = custom_name or "injrec.in"
        if not det_string.endswith(".in"):
            det_string += ".in"
        cat_in = "/mnt/hoy/temporaryFilesROHAN/lephare/inputs/euclid/" + det_string.replace("_DR3", "")
        cat_out = det_string.replace(".in", ".out").replace("_DR3", "")

        f.write(f"CAT_IN\t{cat_in}\n")
        f.write("INP_TYPE\tF\n")
        f.write("CAT_MAG\tAB\n")
        f.write("CAT_FMT\tMEME\n")
        f.write("CAT_LINES\t1,2000000\n")
        f.write("CAT_TYPE\tSHORT\n")
        f.write(f"CAT_OUT\t$LEPHAREDIR/test/{cat_out}\n")
        f.write("PARA_OUT\t$LEPHAREDIR/config/zphot_output.para\n")
        f.write("BD_SCALE\t0\n")
        f.write("GLB_CONTEXT\t-1\n")
        f.write("ERR_FACTOR\t1.0\n")
        f.write("ZPHOTLIB\tGAL_EUC,STAR_EUC\n")
        f.write("ADD_EMLINES\tNO\n")
        f.write("MAG_ABS\t-10.,-30.\n")
        f.write("MAG_REF\t4\n")
        f.write("Z_RANGE\t0.,99.99\n")
        f.write("EBV_RANGE\t0,9\n")
        f.write("ZFIX\tNO\n")
        f.write("Z_INTERP\tYES\n")
        f.write("DZ_WIN\t0.25\n")
        f.write("MIN_THRES\t0.0\n")
        f.write("MABS_METHOD\t3\n")
        f.write("MABS_CONTEXT\t-1\n")
        f.write("MABS_REF\t4\n")
        f.write("MABS_FILT\t1,2,3,4\n")
        f.write("MABS_ZBIN\t0,0.5,1,1.5,2,3,3.5,4\n")
        f.write(f"SPEC_OUT\t{'YES' if spec_out else 'NO'}\n")
        f.write("CHI2_OUT\tNO\n")
        f.write("PDZ_OUT\tNONE\n")
        f.write("PDZ_MABS_FILT\t2,10,14\n")
        f.write("FAST_MODE\tNO\n")
        f.write("AUTO_ADAPT\tNO\n")
