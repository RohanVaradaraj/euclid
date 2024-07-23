#!/bin/bash
# Euclid data access requires login your ESA credentials.
# Usage: sh archive_access <username> <password>.
wget --keep-session-cookies --save-cookies cookies.txt --post-data "username=$1&password=$2" "https://easotf.esac.esa.int/sas-dd/login"
wget --load-cookies cookies.txt -O "EUC_MER_DETECTOR-LAYERING-NIR-Y_TILE101541377-6ED40C_20240201T184941.979007Z_00.00.fits" "https://easotf.esac.esa.int/sas-dd/data?file_name=EUC_MER_DETECTOR-LAYERING-NIR-Y_TILE101541377-6ED40C_20240201T184941.979007Z_00.00.fits&release=sedm&RETRIEVAL_TYPE=FILE"
wget --load-cookies cookies.txt "https://easotf.esac.esa.int/sas-dd/logout"
