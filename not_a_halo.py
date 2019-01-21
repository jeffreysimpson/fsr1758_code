
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table

sns.set_context("paper", font_scale=0.8)

gaia_decam = fits.open("../data/decam_gaia_2deg.fits")
gaia_decam = gaia_decam[1].data
gaia_decam_c = SkyCoord(ra=gaia_decam['ra']*u.degree,
                        dec=gaia_decam['dec']*u.degree, frame='icrs')

gaia_decam = Table(np.vstack((gaia_decam['ra'],
                              gaia_decam['dec'],
                              gaia_decam['gmag']-gaia_decam['imag'],
                              gaia_decam['gmag'],
                              gaia_decam['pmra'],
                              gaia_decam['pmdec'],
                              gaia_decam['bp_rp'],
                              gaia_decam['phot_g_mean_mag'],
                              gaia_decam['parallax'],
                              gaia_decam['locus_sample'])).T,
                   names=('ra', 'dec', 'g_i', 'gmag',
                          'pmra', 'pmdec', 'bp_rp', 'phot_g_mean_mag',
                          'parallax', 'locus'))

cluster_centre = SkyCoord(ra=262.806*u.degree,
                          dec=-39.822*u.degree, frame='icrs')
cluster_pos_idx = gaia_decam_c.separation(cluster_centre) < 0.2*u.deg
cluster_pm_idx = np.sqrt((gaia_decam['pmra']--2.85)**2 +
                         (gaia_decam['pmdec']-2.55)**2) < 1.2

parallax_idx = gaia_decam['parallax'] < 0.3
halo_idx = ~cluster_pos_idx & cluster_pm_idx
locus_idx = gaia_decam['locus'] == 1
likely_cluster_idx = cluster_pos_idx & cluster_pm_idx

cluster_locus_idx = parallax_idx & cluster_pm_idx & locus_idx & cluster_pos_idx
halo_locus_idx = parallax_idx & cluster_pm_idx & locus_idx & ~cluster_pos_idx

list_plots = [['bp_rp', 'phot_g_mean_mag'],
              ['g_i', 'gmag'],
              ['ra', 'dec'],
              ['pmra', 'pmdec'],
              ['g_i'],
              ['parallax']]

axes_labels = [[r'$G_\mathrm{BP}-\mathrm{G_{RP}}$', r'$G$'],
               [r'$g-i$', r'$g$'],
               ['RA (deg)', 'Dec (deg)'],
               [r'$\mu_\mathrm{RA}$ (mas yr$^{-1}$)',
                r'$\mu_\mathrm{Dec}$ (mas yr$^{-1}$)'],
               [r'$g-i$', 'Number of stars'],
               [r'parallax (mas)', 'Number of stars']]

idx_list = [parallax_idx & cluster_pm_idx & ~locus_idx,
            parallax_idx & likely_cluster_idx,
            cluster_locus_idx,
            halo_locus_idx]

plot_kwargs = [dict(s=0.5, color='k', alpha=0.2),
               dict(s=0.5, color='C3', alpha=0.0),
               dict(s=3, color='C0'),
               dict(s=3, color='C1')]

bins_list = [np.arange(0.5, 2.7, 0.07),
             np.arange(-1., 1, 0.05)]

fig, axes = plt.subplots(ncols=2, nrows=3, figsize=(3.32*2, 7.),
                         constrained_layout=True)

for axes_count, ax in enumerate(axes.flatten()):
    ax.set_xlabel(axes_labels[axes_count][0])
    ax.set_ylabel(axes_labels[axes_count][1])
    for idx_count, idx in enumerate(idx_list):
        if axes_count < 4:
            ax.scatter(gaia_decam[list_plots[axes_count][0]][idx],
                       gaia_decam[list_plots[axes_count][1]][idx],
                       **plot_kwargs[idx_count])
        else:
            if idx_count > 1:
                ax.hist(gaia_decam[idx].as_array()[list_plots[axes_count][0]],
                        histtype='step', bins=bins_list[axes_count-4])
#                        **plot_kwargs[idx_count])

axes[1, 0].set_xlim([265.5, 260.1])
axes[1, 0].set_ylim([-41.9, -37.8])

axes[0, 0].set_xlim([0.5, 3.5])
axes[0, 0].set_ylim([19, 13])

axes[0, 1].set_xlim([0.5, 3.0])
axes[0, 1].set_ylim([22, 15])

axes[2, 1].set_xlim([-1, 0.5])


# fig.align_labels()
plt.savefig("../fsr1758_paper/figures/not_a_halo.pdf", bbox_inches='tight')

# plt.show()
