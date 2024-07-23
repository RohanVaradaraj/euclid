#!/bin/bash
# Euclid data access requires login your ESA credentials.
# Usage: sh archive_access <username> <password>.
curl -k -c cookies.txt -X POST -d username=$1 -d password=$2 -L 'https://easotf.esac.esa.int/sas-dd/login'
curl -k -b cookies.txt -o 'EUC_MER_MOSAIC-NIR-Y-FLAG_TILE101545697-672320_20240202T054903.910504Z_00.00.fits' 'https://easotf.esac.esa.int/sas-dd/data?file_name=EUC_MER_MOSAIC-NIR-Y-FLAG_TILE101545697-672320_20240202T054903.910504Z_00.00.fits&release=sedm&RETRIEVAL_TYPE=FILE'
curl -k -b cookies.txt -X POST -d -L 'https://easotf.esac.esa.int/sas-dd/logout'
