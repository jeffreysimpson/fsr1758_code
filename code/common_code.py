#!/usr/bin/env python

"""common_code.py: Loads the tables and provides common masks and others."""

import numpy as np
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


def good_photom_func(table):
    """Use the photometric quality criteria from Evans+2018."""
    bp_rp_excess = table['phot_bp_rp_excess_factor']
    bp_rp = table['bp_rp']
    return ((bp_rp_excess <
             1.3 + 0.06*bp_rp**2) &
            (bp_rp_excess >
             1.0 + 0.015*bp_rp**2))


def good_astrom_func(table):
    """Require the star to have good astrometry."""
    return table['ruwe'] < 1.4


def open_file(file_name):
    gaia_1p5deg = fits.open(file_name)
    fsr1758 = gaia_1p5deg[1].data
    """Returns a couple of common masks."""
    c = SkyCoord(ra=fsr1758['ra']*u.degree,
                 dec=fsr1758['dec']*u.degree, frame='icrs')
    cluster_centre = SkyCoord(ra=262.806*u.degree,
                              dec=-39.822*u.degree, frame='icrs')
    cluster_pos_idx = c.separation(cluster_centre) < 0.2*u.deg
    cluster_pm_idx = np.sqrt((fsr1758['pmra']--2.85)**2 +
                             (fsr1758['pmdec']-2.55)**2) < 1.2
    good_photom_idx = good_photom_func(fsr1758)
    good_astrom_idx = good_astrom_func(fsr1758)
    return (fsr1758, cluster_pos_idx, cluster_pm_idx,
            good_photom_idx, good_astrom_idx)
