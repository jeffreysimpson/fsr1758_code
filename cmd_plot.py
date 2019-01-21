#!/usr/bin/env python

"""cmd_plot.py: Identifies and plots members of FSR1758."""

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import vaex

__author__ = "Jeffrey Simpson"
__copyright__ = "Copyright 2019, Jeffrey Simpson"
__credits__ = ["Jeffrey Simpson"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey Simpson"
__email__ = "jeffrey.simpson@unsw.edu.au"
__status__ = "Development"

sns.set_context("paper", font_scale=0.8)

gaia_1p5deg = fits.open("../data/gaia_fsr1758.fits")
fsr1758 = gaia_1p5deg[1].data
c = SkyCoord(ra=fsr1758['ra']*u.degree,
             dec=fsr1758['dec']*u.degree, frame='icrs')
cluster_centre = SkyCoord(ra=262.806*u.degree,
                          dec=-39.822*u.degree, frame='icrs')
cluster_pos_idx = c.separation(cluster_centre) < 0.2*u.deg
cluster_pm_idx = np.sqrt((fsr1758['pmra']--2.85)**2 +
                         (fsr1758['pmdec']-2.55)**2) < 1.2
large_pm_idx = (c.separation(cluster_centre) < 35*u.arcmin)
radial_velocity_members_idx = fsr1758['radial_velocity'] > 200
has_rv_idx = ~np.isnan(fsr1758['radial_velocity'])

likely_cluster_idx = cluster_pos_idx & cluster_pm_idx

fsr1758_vaex = vaex.from_arrays(ra=fsr1758['ra'],
                                dec=fsr1758['dec'],
                                pmra=fsr1758['pmra'],
                                pmdec=fsr1758['pmdec'],
                                radial_velocity=fsr1758['radial_velocity'],
                                cluster_distance=c.separation(cluster_centre),
                                bp_rp=fsr1758['bp_rp'],
                                phot_g_mean_mag=fsr1758['phot_g_mean_mag'])

def plot(axes, x, y, RV=False, PM=False):
    background_idx = cluster_pm_idx & ~likely_cluster_idx & ~has_rv_idx
    pm_field_idx = (likely_cluster_idx &
                    ~radial_velocity_members_idx & ~has_rv_idx)
    pm_rv_idx = cluster_pm_idx & radial_velocity_members_idx
    pm_not_rv_idx = cluster_pm_idx & has_rv_idx & ~radial_velocity_members_idx
    axes.scatter(x[background_idx],
                 y[background_idx],
                 alpha=0.5, s=3, c='#984ea3', lw=0)
    axes.scatter((x[pm_field_idx]),
                 (y[pm_field_idx]),
                 alpha=0.8, s=3, c='#4daf4a', lw=0)
    axes.scatter((x[pm_rv_idx]),
                 (y[pm_rv_idx]),
                 marker='*', alpha=1., s=60, c='#e41a1c', lw=0)
    axes.scatter((x[pm_not_rv_idx]),
                 (y[pm_not_rv_idx]),
                 marker='*', edgecolor='#e41a1c',
                 facecolor='None', alpha=1., s=30, lw=0.4)

ra = 262.806
dec = -39.822
box = 0.8
vaex_selection = f'(ra > {ra}-{box}) & (ra < {ra}+{box}) & (dec > {dec}-{box}) & (dec < {dec}+{box})'

fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(3.32*2, 5.5))
#                          constrained_layout=True)#,sharex='col',sharey='col')
plt.sca(axes[0, 0])
fsr1758_vaex.plot(fsr1758_vaex.ra, fsr1758_vaex.dec,
                  selection=[vaex_selection],
                  colorbar=False,
                  colormap="Greys")
plt.sca(axes[0, 1])
fsr1758_vaex.plot(fsr1758_vaex.pmra, fsr1758_vaex.pmdec,
                  selection=[vaex_selection],
                  colorbar=False,
                  colormap="Greys")
plt.sca(axes[1, 0])
fsr1758_vaex.plot(fsr1758_vaex.cluster_distance, fsr1758_vaex.radial_velocity,
                  selection=[vaex_selection],
                  colorbar=False,
                  colormap="Greys")
plt.sca(axes[1, 1])
fsr1758_vaex.plot(fsr1758_vaex.bp_rp, fsr1758_vaex.phot_g_mean_mag,
                  selection=[vaex_selection],
                  colorbar=False,
                  colormap="Greys")
plot(axes[0, 0], fsr1758['ra'], fsr1758['dec'])
plot(axes[0, 1], fsr1758['pmra'], fsr1758['pmdec'], PM=True, RV=True)
plot(axes[1, 0],
     c.separation(cluster_centre), fsr1758['radial_velocity'], RV=True)
plot(axes[1, 1], fsr1758['bp_rp'], fsr1758['phot_g_mean_mag'], RV=True)


axes[0, 0].set_xlim(ra-box, ra+box)
axes[0, 0].set_ylim(dec-box, dec+box)

axes[0, 1].set_xlim(-2.85-7, -2.85+7)
axes[0, 1].set_ylim(2.55-8, 2.55+5.5)

axes[1, 1].set_xlim(0.5, 3.5)
axes[1, 1].set_ylim(21, 10)
fig.align_labels()
plt.savefig("../fsr1758_paper/figures/cmd.pdf", bbox_inches='tight')

# plt.show()
