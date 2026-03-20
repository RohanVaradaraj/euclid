import numpy as np
import pyvo
import requests
from pyvo.dal.adhoc import DatalinkResults
from lsst.daf.butler import Butler

token = "gt-4qgBILM-wppoXGn4-IFh5Q.2AfwKN0PeTW26ONidp8EOg"
session = requests.Session()
session.headers["Authorization"] = f"Bearer {token}"

rsp_sia_url = "https://data.lsst.cloud/api/sia/dp1/query"
sia_service = pyvo.dal.SIA2Service(
    rsp_sia_url, session=session, check_baseurl=False
)

circle = (53.076, -28.110, 1.5)
max_rec = 5

results = sia_service.search(
    pos=circle,
    calib_level=3,
    maxrec=max_rec,
    data_type="image",
    res_format="lsst.deep_coadd",
)

table = results.to_table()
table = table[table["dataproduct_subtype"] == "lsst.deep_coadd"]

# Adjust these for your environment
butler = Butler("/path/to/dp1/repo", collections=["YOUR_COLLECTION"])

for row in table:
    tract = row["lsst_tract"]
    patch = row["lsst_patch"]
    band = row["lsst_band"]

    data_id = {"tract": tract, "patch": patch, "band": band}

    # This is the calibrated coadd exposure's photometric calibration
    photo_calib = butler.get("deepCoadd_calexp.photoCalib", dataId=data_id)

    flux_mag0 = photo_calib.getFluxMag0()
    zp_ab = 2.5 * np.log10(flux_mag0)

    print(f"tract={tract} patch={patch} band={band}")
    print(f"  fluxMag0 = {flux_mag0:.6g}")
    print(f"  zero point (AB) = {zp_ab:.6f} mag")