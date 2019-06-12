#!/usr/bin/env python

"""common_code.py: Loads the tables and provides common masks and others."""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits


__author__ = "Jeffrey Simpson"
__copyright__ = "Copyright 2019, Jeffrey Simpson"
__credits__ = ["Jeffrey Simpson"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey Simpson"
__email__ = "jeffrey.simpson@unsw.edu.au"
__status__ = "Development"


def good_photom_idx(table, FITS=False):
    """Use the photometric quality criteria from Evans+2018."""
    if FITS:
        bp_rp_excess = table[1].data['phot_bp_rp_excess_factor']
        bp_rp = table[1].data['bp_rp']
    else:
        bp_rp_excess = table.phot_bp_rp_excess_factor
        bp_rp = table.bp_rp
    return ((bp_rp_excess <
             1.3 + 0.06*bp_rp**2) &
            (bp_rp_excess >
             1.0 + 0.015*bp_rp**2))


def good_astrom_idx(table, FITS=False):
    """Require the star to have good astrometry."""
    if FITS:
        ruwe = table[1].data['ruwe']
    else:
        ruwe = table.ruwe
    return ruwe < 1.4


def common_idx(table):
    """Returns a couple of common masks."""
    c = SkyCoord(ra=table['ra']*u.degree,
                 dec=table['dec']*u.degree, frame='icrs')
    cluster_centre = SkyCoord(ra=262.806*u.degree,
                              dec=-39.822*u.degree, frame='icrs')
    cluster_pos_idx = c.separation(cluster_centre) < 0.2*u.deg
    cluster_pm_idx = np.sqrt((table['pmra']--2.85)**2 +
                             (table['pmdec']-2.55)**2) < 1.2
    return cluster_pos_idx, cluster_pm_idx
