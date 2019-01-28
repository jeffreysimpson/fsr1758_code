#!/usr/bin/env python

"""cmd_plot.py: Identifies and plots members of FSR1758."""

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# import vaex

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
rrlyraeab_idx = fsr1758['Mode'] == "RRab"
rrlyraeac_idx = fsr1758['Mode'] == "RRc"


panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

xy_values = [['ra', 'dec'],
             ['pmra', 'pmdec'],
             ['bp_rp', 'phot_g_mean_mag'],
             ['g_i', 'gmag'],
             ['cluster_distance', 'radial_velocity'],
             ['parallax', ]]

axes_labels = [['RA (deg)', 'Dec (deg)'],
               [r'$\mu_\mathrm{RA}$ (mas yr$^{-1}$)',
                r'$\mu_\mathrm{Dec}$ (mas yr$^{-1}$)'],
               [r'$G_\mathrm{BP}-\mathrm{G_{RP}}$', r'$G$'],
               [r'$g-i$', r'$g$'],
               ['Angular distance (deg)', r'$v_r$ (km s$^{-1}$'],
               ['parallax (mas)', 'Number of stars']]

idx_list = [cluster_pm_idx & ~cluster_pos_idx & ~has_rv_idx,
            cluster_pm_idx & cluster_pos_idx & ~has_rv_idx,
            cluster_pm_idx & cluster_pos_idx & radial_velocity_members_idx,
            cluster_pm_idx & cluster_pos_idx & has_rv_idx & ~radial_velocity_members_idx,
            cluster_pm_idx & ~cluster_pos_idx & radial_velocity_members_idx,
            cluster_pm_idx & ~cluster_pos_idx & has_rv_idx & ~radial_velocity_members_idx,
            cluster_pm_idx & cluster_pos_idx & rrlyraeab_idx,
            cluster_pm_idx & ~cluster_pos_idx & rrlyraeab_idx,
            cluster_pm_idx & cluster_pos_idx & rrlyraeac_idx,
            cluster_pm_idx & ~cluster_pos_idx & rrlyraeac_idx]

plot_kwargs = [dict(alpha=0.5, s=2, c='#4daf4a', lw=0),
               dict(alpha=0.8, s=4, c='#984ea3', lw=0),
               dict(marker='*', alpha=1., s=60, c='#e41a1c', lw=0),
               dict(marker='*', edgecolor='#e41a1c',
                    facecolor='None', alpha=1., s=30, lw=0.4),
               dict(marker='*', alpha=1., s=60, c='k', lw=0),
               dict(marker='*', edgecolor='k',
                    facecolor='None', alpha=1., s=30, lw=0.4),
               dict(marker='s', alpha=1., s=10, c='C0', lw=0),
               dict(marker='s', edgecolor='C0',
                    facecolor='None', alpha=1., s=10, lw=0.4),
               dict(marker='<', alpha=1., s=10, c='k', lw=0),
               dict(marker='<', edgecolor='k',
                    facecolor='None', alpha=1., s=10, lw=0.4)]

ra = 262.806
dec = -39.822
box = 1.3
plot_limits = [[[ra-box, ra+box], [dec-(box-0.25), dec+(box-0.25)]],
               [[-2.85-1.3, -2.85+1.3], [2.55-1.3, 2.55+1.3]],
               [[0.7, 3.1], [20.5, 10]],
               [[0.5, 3.0], [21.5, 13]],
               [[0, 1.1], [-160, 250]],
               [[-2, 2], []]]

fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(3.32*2, 4.5),
                         constrained_layout=True)

for axes_count, ax in enumerate(axes.flatten()):
    ax.set_xlabel(axes_labels[axes_count][0])
    ax.set_ylabel(axes_labels[axes_count][1])
    ax.set_xlim(plot_limits[axes_count][0])
    ax.annotate(panel_labels[axes_count], (0.85, 0.92),
                xycoords='axes fraction')
    if axes_count != 5:
        ax.set_ylim(plot_limits[axes_count][1])
    for idx_count, idx in enumerate(idx_list):
        if axes_count != 5:
            ax.scatter(fsr1758[xy_values[axes_count][0]][idx],
                       fsr1758[xy_values[axes_count][1]][idx],
                       **plot_kwargs[idx_count])
        else:
            if idx_count < 2:
                ax.hist(fsr1758[xy_values[axes_count][0]][idx],
                        bins=np.arange(-2, 2, 0.05), histtype='step',
                        color=plot_kwargs[idx_count]['c'])
            if idx_count > 1:
                if idx_count in [2, 3]:
                    offset = 10
                if idx_count in [4, 5]:
                    offset = 15
                ax.scatter(fsr1758[xy_values[axes_count][0]][idx],
                           np.random.randn(np.sum(idx))+offset,
                           **plot_kwargs[idx_count], zorder=100)

# fig.align_labels()
plt.savefig("../fsr1758_paper/figures/cmd.pdf", bbox_inches='tight')
