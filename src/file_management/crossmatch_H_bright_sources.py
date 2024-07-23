"""
crossmatch_H_bright_sources.py

Get the Euclid MER IDs for our H-bright z=3 sources in COSMOS using command line queries.

Created: Tuesday 7th May 2024.

Run the following before the code:
curl -k -c cookies.txt -X POST -d username=USER -d password=PASSWORD -L "https://easotf.esac.esa.int/tap-server/login"
"""

import numpy as np
from astropy.table import Table
from pathlib import Path
import os 

# Load the H-bright sources
cat_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues'
cat_name = 'H_bright_z3_to_5.fits'

t = Table.read(cat_dir / cat_name)

for i in range(len(t)):

    # Get the RA and Dec
    ra = t['RA'][i]
    dec = t['DEC'][i]
    id = t['UID'][i]

    #### Make the query string ####

    #! synchronous query
    base_string = "https://easotf.esac.esa.int/tap-server/tap/sync?LANG=ADQL&FORMAT=csv&REQUEST=doQuery&" 

    #! asynchronous query
    #base_string = "PHASE=run&LANG=ADQL&REQUEST=doQuery&"                                                  

    # SELECT part
    select_string = "query=SELECT+mer_catalogue.object_id,mer_catalogue.right_ascension,mer_catalogue.declination+FROM+catalogue.mer_catalogue"

    # WHERE part
    where_string = f"+WHERE+CONTAINS\(POINT\('ICRS',mer_catalogue.right_ascension,mer_catalogue.declination\),CIRCLE\('ICRS',{ra},{dec},2/3600\)\)=1"

    # Combine the strings
    query_string = base_string + select_string + where_string

    #! synchronous query
    print(f'curl -k -b cookies.txt -X GET {query_string} > out_{id}.csv')
    os.system(f'curl -k -b cookies.txt -X GET {query_string} > out_{id}.csv') 

    #! asynchronous query
    #print(f'curl -i -k -b cookies.txt -X POST --data {query_string} https://easotf.esac.esa.int/tap-server/tap/async')
    #os.system(f'curl -i -k -b cookies.txt -X POST --data "{query_string}" "https://easotf.esac.esa.int/tap-server/tap/async"')                    

    # Write RA, DEC to a .csv file, with each row being RA, DEC with no space in between REPLACED BY PLUS, e.g. 250.4+0.3 for RA=250.4, DEC=0.3
    # with open(f'ra_dec.csv', 'a') as f:
    #     f.write(f'{ra}+{dec}\n')



