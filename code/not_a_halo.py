#!/usr/bin/env python

"""not_a_halo.py: Purported halo around FSR1758."""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import common_code

__author__ = "Jeffrey Simpson"
__copyright__ = "Copyright 2019, Jeffrey Simpson"
__credits__ = ["Jeffrey Simpson"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey Simpson"
__email__ = "jeffrey.simpson@unsw.edu.au"
__status__ = "Development"

sns.set_context("paper", font_scale=0.8)

(fsr1758,
 cluster_pos_idx,
 cluster_pm_idx,
 good_photom_idx,
 good_astrom_idx) = common_code.open_file("../data/decam_gaia_2deg.fits")
parallax_cut = fsr1758['parallax'] < 0.3
locus_sample_idx = fsr1758['locus_sample']
likely_cluster_idx = cluster_pos_idx & cluster_pm_idx

field_not_locus_idx = ~cluster_pos_idx & ~locus_sample_idx
cluster_not_locus_idx = cluster_pos_idx & ~locus_sample_idx
field_locus_idx = ~cluster_pos_idx & locus_sample_idx
cluster_locus_idx = cluster_pos_idx & locus_sample_idx

# Require for all that they have good astrometry
idx_list = [idx & parallax_cut & cluster_pm_idx & good_astrom_idx for idx in [
    field_not_locus_idx, cluster_not_locus_idx,
    field_locus_idx, cluster_locus_idx]]

panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

xy_values = [['ra', 'dec'],
             ['pmra', 'pmdec'],
             ['bp_rp', 'phot_g_mean_mag'],
             ['g_i', 'gmag'],
             ['g_i', ],
             ['parallax', 'phot_g_mean_mag']]

axes_labels = [['RA (deg)', 'Dec (deg)'],
               [r'$\mu_\mathrm{RA}$ (mas yr$^{-1}$)',
                r'$\mu_\mathrm{Dec}$ (mas yr$^{-1}$)'],
               [r'$G_\mathrm{BP}-\mathrm{G_{RP}}$', r'$G$'],
               [r'$g-i$', r'$g$'],
               [r'$g-i$', 'Number of stars'],
               [r'$\varpi$ (mas)', r'$G$']]

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
               [[-0.5, 0.35], [19., 13]]]

fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(3.32*2, 4.5),
                         constrained_layout=True)

for axes_count, ax in enumerate(axes.flatten()):
    ax.set_xlabel(axes_labels[axes_count][0])
    ax.set_ylabel(axes_labels[axes_count][1])
    ax.set_xlim(plot_limits[axes_count][0])
    ax.annotate(panel_labels[axes_count], (0.85, 0.92),
                xycoords='axes fraction')
    if axes_count != 4:
        ax.set_ylim(plot_limits[axes_count][1])
    for idx_count, idx in enumerate(idx_list):
        if axes_count == 2:
            idx = idx & good_photom_idx
        if (axes_count == 5) and (idx_count < 2):
            continue
        if axes_count != 4:
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
plt.savefig("../paper/figures/not_a_halo.pdf", bbox_inches='tight')
