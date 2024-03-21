"""
cutout_codes.py

Functions and classes used to make cutouts

Created: Thursday 21st March 2024.
"""

from shapely.geometry import Point, Polygon
import numpy as np
from pathlib import Path
from astropy.io import ascii




def isCoordInSurveyFootprints(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """
    Check if a coordinate lies in the footprints of:
        1) Euclid On-the-fly COSMOS release
        2) COSMOS-Web
        3 )PRIMER

    Parameters
    ----------
    ra : np.array
        Right ascension of the coordinate in degrees.
    dec : np.array
        Declination of the coordinate in degrees.


    Returns
    -------
    np.ndarray
        Array with rows corresponding to the input coordinates and columns corresponding to the footprints.
        The value is '0' if it is not in the footprint.
        If it is in the footprint then: 
            for Euclid it takes a numeric string value of the tile label I have assigned (1-22).
            for CWEB it takes a numeric string value of the tile label Nathan has assigned (4A, 5A, 5B, 6A, etc.)
            for PRIMER it takes the value '1'.

    """

    # Check if ra and dec are arrays
    ra = np.array(ra)
    dec = np.array(dec)

    assert ra.shape == dec.shape, "ra and dec arrays must have the same shape"

    points = [Point(ra, dec) for ra, dec in zip(ra, dec)]

    # Load the footprints
    euclid_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_footprint.npy')
    primer_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'primer_footprint.npy')
    cweb_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_footprint.npy')

    # Load the labels for CWEB and Euclid
    euclid_labels = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_labels.npy')
    cweb_labels = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_labels.npy')

    # Create shapely polygons
    euclid_polygons = [Polygon(footprint) for footprint in euclid_footprint]
    cweb_polygons = [Polygon(footprint) for footprint in cweb_footprint]
    primer_polygons = Polygon(primer_footprint)

    points_in_euclid = []
    points_in_cweb = []
    points_in_primer = []

    for point in points:

        in_euclid = 0
        in_cweb = 0
        in_primer = 0

        #! Check if the point lies in any of the Euclid polygons
        for i, polygon in enumerate(euclid_polygons):
            if polygon.contains(point):
                in_euclid = euclid_labels[i]
                break

        #! Check if the point lies in any of the COSMOS-Web polygons
        for i, polygon in enumerate(cweb_polygons):
            if polygon.contains(point):
                in_cweb = cweb_labels[i]
                break

        #! Check if point lies in PRIMER footprint
        if primer_polygons.contains(point):
            in_primer = 1

        # Append results of search to lists
        points_in_euclid.append(in_euclid)
        points_in_cweb.append(in_cweb)
        points_in_primer.append(in_primer)
    
    # Convert lists to numpy arrays
    euclid_mask = np.array(points_in_euclid)
    cweb_mask = np.array(points_in_cweb)
    primer_mask = np.array(points_in_primer)

    return np.column_stack((euclid_mask, cweb_mask, primer_mask))




t = ascii.read(Path.cwd().parent.parent / 'data' / 'mosaic' / 'REBELS.csv', format='csv')

ra = t['RA']
dec = t['Dec']


class Cutout:
    """
    Class to make cutouts of images in our various fields
    """

    
