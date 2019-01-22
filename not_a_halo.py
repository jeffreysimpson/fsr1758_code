#!/usr/bin/env python

"""not_a_halo.py: Purported halo around FSR1758."""

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

__author__ = "Jeffrey Simpson"
__copyright__ = "Copyright 2019, Jeffrey Simpson"
__credits__ = ["Jeffrey Simpson"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey Simpson"
__email__ = "jeffrey.simpson@unsw.edu.au"
__status__ = "Development"

sns.set_context("paper", font_scale=0.8)

gaia_1p5deg = fits.open("../data/decam_gaia_2deg.fits")

fsr1758 = gaia_1p5deg[1].data
c = SkyCoord(ra=fsr1758['ra']*u.degree,
             dec=fsr1758['dec']*u.degree, frame='icrs')
cluster_centre = SkyCoord(ra=262.806*u.degree,
                          dec=-39.822*u.degree, frame='icrs')
cluster_pos_idx = c.separation(cluster_centre) < 0.2*u.deg
cluster_pm_idx = np.sqrt((fsr1758['pmra']--2.85)**2 +
                         (fsr1758['pmdec']-2.55)**2) < 1.2
parallax_cut = fsr1758['parallax'] < 0.3
large_pm_idx = (c.separation(cluster_centre) < 35*u.arcmin)
radial_velocity_members_idx = fsr1758['radial_velocity'] > 200
has_rv_idx = ~np.isnan(fsr1758['radial_velocity'])
locus_sample_idx = fsr1758['locus_sample']

likely_cluster_idx = cluster_pos_idx & cluster_pm_idx


panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

xy_values = [['ra', 'dec'],
             ['pmra', 'pmdec'],
             ['bp_rp', 'phot_g_mean_mag'],
             ['g_i', 'gmag'],
             ['g_i', ],
             ['parallax', ]]

axes_labels = [['RA (deg)', 'Dec (deg)'],
               [r'$\mu_\mathrm{RA}$ (mas yr$^{-1}$)',
                r'$\mu_\mathrm{Dec}$ (mas yr$^{-1}$)'],
               [r'$G_\mathrm{BP}-\mathrm{G_{RP}}$', r'$G$'],
               [r'$g-i$', r'$g$'],
               [r'$g-i$', 'Number of stars'],
               ['parallax (mas)', 'Number of stars']]

idx_list = [(parallax_cut & cluster_pm_idx & ~likely_cluster_idx &
             ~has_rv_idx & ~locus_sample_idx),
            (parallax_cut & likely_cluster_idx & ~radial_velocity_members_idx &
             ~has_rv_idx & ~locus_sample_idx),
            (parallax_cut & cluster_pm_idx & ~likely_cluster_idx &
             ~has_rv_idx & locus_sample_idx),
            (parallax_cut & likely_cluster_idx & ~radial_velocity_members_idx &
             ~has_rv_idx & locus_sample_idx)]

plot_kwargs = [dict(alpha=0.5/3, s=2/2, c='#4daf4a', lw=0),
               dict(alpha=0.8/3, s=4/2, c='#984ea3', lw=0),
               dict(alpha=0.7, s=4, c='#4daf4a', lw=0),
               dict(alpha=0.7, s=4, c='#984ea3', lw=0)]

ra = 262.806
dec = -39.822
box = 2.5
plot_limits = [[[ra-box, ra+box], [dec-(box-0.1), dec+(box-0.1)]],
               [[-2.85-1.3, -2.85+1.3], [2.55-1.3, 2.55+1.3]],
               [[0.7, 3.0], [20.5, 13]],
               [[0.5, 3.0], [21.5, 16]],
               [[0.5, 3.0], [-160, 250]],
               [[-0.5, 0.5], []]]

fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(3.32*2, 7.),
                         constrained_layout=True)

for axes_count, ax in enumerate(axes.flatten()):
    ax.set_xlabel(axes_labels[axes_count][0])
    ax.set_ylabel(axes_labels[axes_count][1])
    ax.set_xlim(plot_limits[axes_count][0])
    ax.annotate(panel_labels[axes_count], (0.90, 0.90),
                xycoords='axes fraction')
    if axes_count < 4:
        ax.set_ylim(plot_limits[axes_count][1])
    for idx_count, idx in enumerate(idx_list):
        if axes_count < 4:
            ax.scatter(fsr1758[xy_values[axes_count][0]][idx],
                       fsr1758[xy_values[axes_count][1]][idx],
                       **plot_kwargs[idx_count])
        else:
            if (idx_count > 1) & (axes_count == 4):
                ax.hist(fsr1758[xy_values[axes_count][0]][idx],
                        bins=np.arange(0.5, 3.0, 0.07), histtype='step',
                        color=plot_kwargs[idx_count]['c'])
            if (idx_count > 1) & (axes_count == 5):
                ax.hist(fsr1758[xy_values[axes_count][0]][idx],
                        bins=np.arange(-0.5, 0.5, 0.05), histtype='step',
                        color=plot_kwargs[idx_count]['c'])

# fig.align_labels()
plt.savefig("../fsr1758_paper/figures/not_a_halo.pdf", bbox_inches='tight')
