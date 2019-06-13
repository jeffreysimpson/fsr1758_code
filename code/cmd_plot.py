#!/usr/bin/env python

"""cmd_plot.py: Identifies and plots members of FSR1758."""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import common_code

__author__ = "Jeffrey Simpson"
__copyright__ = "Copyright 2019, Jeffrey Simpson"
__credits__ = ["Jeffrey Simpson"]
__license__ = "MIT"
__version__ = "0.0.2"
__maintainer__ = "Jeffrey Simpson"
__email__ = "jeffrey.simpson@unsw.edu.au"
__status__ = "Development"


(fsr1758,
 cluster_pos_idx,
 cluster_pm_idx,
 good_photom_idx,
 good_astrom_idx) = common_code.open_file("../data/gaia_fsr1758.fits")
fast_enough_rv_idx = fsr1758['radial_velocity'] > 200
good_rv_idx = fsr1758['rv_nb_transits'] >= 5
has_rv_idx = ~np.isnan(fsr1758['radial_velocity'])
likely_cluster_idx = cluster_pos_idx & cluster_pm_idx

# The field stars and cluster stars without radial velocities
field_stars_idx = ~cluster_pos_idx & ~has_rv_idx
cluster_stars_idx = cluster_pos_idx & ~has_rv_idx

# The cluster stars with RVs
cluster_star_rv_idx = cluster_pos_idx & has_rv_idx & good_rv_idx
cluster_stars_right_rv_idx = cluster_star_rv_idx & fast_enough_rv_idx
cluster_stars_wrong_rv_idx = cluster_star_rv_idx & ~fast_enough_rv_idx

# The field stars with RVs
field_star_rv_idx = ~cluster_pos_idx & has_rv_idx & good_rv_idx
field_stars_right_rv_idx = field_star_rv_idx & fast_enough_rv_idx
field_stars_wrong_rv_idx = field_star_rv_idx & ~fast_enough_rv_idx

# Require for all that they have good astrometry
idx_list = [idx & cluster_pm_idx & good_astrom_idx for idx in [
    field_stars_idx, cluster_stars_idx,
    cluster_stars_right_rv_idx, cluster_stars_wrong_rv_idx,
    field_stars_right_rv_idx, field_stars_wrong_rv_idx]]


print(np.sum(likely_cluster_idx & good_astrom_idx))

print("\\hline")
print(r"source\_id & RA & Dec & $\varpi$ & $v_r$ & \bprp & $G$ & $\mu_\mathrm{RA}$ & $\mu_\mathrm{Dec}$ & ang. dist \\")
print(r" & (deg) & (deg) & (mas) & (\kms) &  &  & (\masyr) & (\masyr) & (deg) \\")
print(r"\hline")
for idx in idx_list[2:]:
    for star in fsr1758[idx]:
        print_str = f"{star['source_id']} & "
        print_str += rf"${star['ra']:0.3f}$ & "
        print_str += rf"${star['dec']:0.3f}$ & "
        print_str += rf"${star['parallax']:0.3f}$ & "
        print_str += rf"${star['radial_velocity']:0.2f}$ & "
    #     print_str += f"{star['rv_guess']:0.1f} & "
        print_str += rf"${star['bp_rp']:0.2f}$ & "
        print_str += rf"${star['phot_g_mean_mag']:0.2f}$ & "
        print_str += rf"${star['pmra']:0.2f}$ & "
        print_str += rf"${star['pmdec']:0.2f}$ & "
        print_str += rf"${star['cluster_distance']:0.2f}$"
        print(f"{print_str} \\\\")
print(r"\hline")

panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

xy_values = [['ra', 'dec'],
             ['pmra', 'pmdec'],
             ['bp_rp', 'phot_g_mean_mag'],
             ['g_i', 'gmag'],
             ['cluster_distance', 'radial_velocity'],
             ['parallax', 'phot_g_mean_mag']]

axes_labels = [['RA (deg)', 'Dec (deg)'],
               [r'$\mu_\mathrm{RA}$ (mas yr$^{-1}$)',
                r'$\mu_\mathrm{Dec}$ (mas yr$^{-1}$)'],
               [r'$G_\mathrm{BP}-\mathrm{G_{RP}}$', r'$G$'],
               [r'$g-i$', r'$g$'],
               ['Angular distance (deg)', r'$v_r$ (km s$^{-1}$'],
               [r'$\varpi$ (mas)', r'$G$']]


plot_kwargs = [dict(alpha=0.5, s=2, c='#4daf4a', lw=0),
               dict(alpha=0.8, s=4, c='#984ea3', lw=0),
               dict(marker='*', alpha=1., s=60, c='#e41a1c', lw=0),
               dict(marker='*', edgecolor='#e41a1c',
                    facecolor='None', alpha=1., s=30, lw=0.4),
               dict(marker='*', alpha=1., s=60, c='k', lw=0),
               dict(marker='*', edgecolor='k',
                    facecolor='None', alpha=1., s=30, lw=0.4)
               ]

ra = 262.806
dec = -39.822
box = 1.3
plot_limits = [[[ra-box, ra+box], [dec-(box-0.25), dec+(box-0.25)]],
               [[-2.85-1.3, -2.85+1.3], [2.55-1.3, 2.55+1.3]],
               [[0.7, 3.1], [20.5, 10]],
               [[0.5, 3.0], [21.5, 13]],
               [[0.0, 1.1], [-160, 250]],
               [[-2, 2], [20.5, 10]]]

sns.set_context("paper", font_scale=0.8)
fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(3.32*2, 4.5),
                         constrained_layout=True)

for axes_count, ax in enumerate(axes.flatten()):
    ax.set_xlabel(axes_labels[axes_count][0])
    ax.set_ylabel(axes_labels[axes_count][1])
    ax.set_xlim(plot_limits[axes_count][0])
    ax.annotate(panel_labels[axes_count], (0.85, 0.92),
                xycoords='axes fraction')
    ax.set_ylim(plot_limits[axes_count][1])
    for idx_count, idx in enumerate(idx_list):
        if axes_count == 2:
            # Require that the CMD only has stars with good photometry.
            idx = idx & good_photom_idx
        ax.scatter(fsr1758[xy_values[axes_count][0]][idx],
                   fsr1758[xy_values[axes_count][1]][idx],
                   **plot_kwargs[idx_count])
plt.savefig("../paper/figures/cmd.pdf", bbox_inches='tight')
