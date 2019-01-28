#!/usr/bin/env python

"""cluster_comparison.py: Plot cluster orbit comparison for FSR1758 paper."""

import astropy.coordinates as coord
import astropy.units as u
import gala.dynamics as gd
import gala.potential as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

sns.set_context("paper", font_scale=0.8)

__author__ = "Jeffrey Simpson"
__copyright__ = "Copyright 2019, Jeffrey Simpson"
__credits__ = ["Jeffrey Simpson"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey Simpson"
__email__ = "jeffrey.simpson@unsw.edu.au"
__status__ = "Development"


def compute_orbit(icrs, dt=-0.1*u.Myr, n_steps=50000):
    c_icrs = icrs.transform_to(gc_frame).cartesian
    object_phase_space = gd.PhaseSpacePosition(
        pos=c_icrs.xyz,
        vel=c_icrs.differentials['s'].d_xyz)
    return gp.Hamiltonian(pot).integrate_orbit(object_phase_space, dt=dt, n_steps=n_steps)

pot = gp.MilkyWayPotential()
v_sun = coord.CartesianDifferential([11, 248, 7.25]*u.km/u.s)
gc_frame = coord.Galactocentric(galcen_distance=8.2*u.kpc,
                                z_sun=25.*u.pc,
                                galcen_v_sun=v_sun)

clusters = pd.read_csv("../data/vasilev_catalogue.csv")

clusters = clusters[~np.in1d(clusters.name,['E1',"Pal3","Pal4","Crater","Eridanus","NGC6121"])] #excluded due to large errors in PM

# clusters = clusters[np.in1d(clusters.name,['Crater'])]
clusters_icrs = coord.ICRS(ra=clusters.ra.values * u.deg,
                           dec=clusters.dec.values * u.deg,
                           distance=clusters.dist.values * u.kpc,
                           pm_ra_cosdec=clusters.pmra.values * u.mas/u.yr,
                           pm_dec=clusters.pmdec.values * u.mas/u.yr,
                           radial_velocity=clusters.rv.values * u.km/u.s)

fsr1758_icrs = coord.ICRS(ra=262.806 * u.deg, dec=-39.822 * u.deg,
                          distance=11.5 * u.kpc,
                          pm_ra_cosdec=-2.85 * u.mas/u.yr,
                          pm_dec=2.55 * u.mas/u.yr,
                          radial_velocity=227 * u.km/u.s)

cluster_orbit = compute_orbit(clusters_icrs, dt=-0.2*u.Myr, n_steps=50000)
fsr1758_orbit = compute_orbit(fsr1758_icrs, dt=-0.2*u.Myr, n_steps=50000)

cluster_peri = cluster_orbit.pericenter(approximate=True)
cluster_apo = cluster_orbit.apocenter(approximate=True)
cluster_e = cluster_orbit.eccentricity(approximate=True)
cluster_r = cluster_orbit.spherical.distance[0]

from matplotlib.ticker import ScalarFormatter

annotate_pos = {"lower_left": dict(horizontalalignment='right', verticalalignment='top',
                         xytext=(-1.1, -1.1), textcoords='offset pixels'),
                    "lower_right": dict(verticalalignment='top',
                         xytext=(1.1, 1.1), textcoords='offset pixels'),
                    "upper_left": dict(horizontalalignment='right', verticalalignment='bottom',
                         xytext=(1.1, 1.1), textcoords='offset pixels'),
                    "upper_right": dict(xytext=(1.1, 1.1), textcoords='offset pixels')}

cluster_name_dict = {"Terzan5": ["Terzan 5", annotate_pos['upper_right']],
                     "NGC5139": ["ω Cen", annotate_pos['lower_right']],
                     "NGC6715": ["M54", annotate_pos['lower_left']],
                     "NGC6656": ["M22", annotate_pos['upper_right']],
                     "NGC7089": ["M2", annotate_pos['upper_right']],
                     "Djorg1": ["Djorgovski 1", annotate_pos['lower_left']],
                     "Terzan10": ["Terzan 10", annotate_pos['lower_left']],
                     "IC4499": ["IC4499", annotate_pos['lower_left']],
                     "NGC2419": ["NGC2419", annotate_pos['lower_left']],
                     "FSR1758": ["FSR1758", annotate_pos['upper_right']],
                     "NGC6544": ["NGC6544", annotate_pos['upper_right']],
                     "NGC3201": ["NGC3201", annotate_pos['upper_left']],
                     "NGC362": ["NGC362", annotate_pos['upper_left']],
                     "NGC1851": ["NGC1851", annotate_pos['upper_left']],
                     "NGC1904": ["NGC1904", annotate_pos['upper_left']],
                     "NGC2298": ["NGC2298", annotate_pos['upper_left']],
                     "NGC6341": ["NGC6341", annotate_pos['upper_left']],
                     "NGC4833": ["NGC4833", annotate_pos['upper_left']],
                     "NGC5286": ["NGC5286", annotate_pos['upper_left']],
                     "NGC6205": ["NGC6205", annotate_pos['upper_left']],
                     "NGC4833": ["NGC4833", annotate_pos['upper_left']],
                     "NGC6779": ["NGC6779", annotate_pos['upper_left']]}

def plot_orbits(axes, peri, apo, e, Lz, r, s=1, marker='.', INTEREST = False, name=None):
    axes[0].scatter(peri, apo, s=s, marker=marker)
    axes[2].scatter(e, apo, s=s, marker=marker)
    axes[3].scatter(Lz, apo, s=s, marker=marker)
    axes[1].scatter(r, apo, s=s, marker=marker)
    if INTEREST:
        axes[0].annotate(cluster_name_dict[name][0], (peri.value, apo.value), fontsize=4,
                     bbox=dict(boxstyle='square', fc='w', ec='none',alpha=0.3, pad=-0.2),
                     **cluster_name_dict[name][1])


fig, axes = plt.subplots(1, 4, figsize=(3.32*2, 2.0), sharey=True)

panel_labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
for axes_count, ax in enumerate(axes.flatten()):
    ax.annotate(panel_labels[axes_count], (0.10, 0.90),
            xycoords='axes fraction')

plot_orbits(axes,
            cluster_peri, cluster_apo, cluster_e,
            cluster_orbit.angular_momentum()[2][0,:].to(u.kpc*u.km/u.s), cluster_r)


for ax in axes[0:2]:
    ax.plot([0,200],[0,200],color='k',lw=0.5)
    ax.axvspan(xmin=0,xmax=3,alpha=0.1,color='C1',lw=0.)

# for cluster_interest in []:
for cluster_interest in ["Terzan5","NGC5139","NGC6715","NGC6656",
                         "NGC7089","Djorg1","Terzan10","IC4499","NGC2419","NGC3201"]:

# for cluster_interest in ['NGC362','NGC1851','NGC1904','NGC2298','NGC4833','NGC5139',
#                          'NGC5286','NGC6205','NGC6341','NGC6779']:
    cluster_idx = np.in1d(clusters.name, [cluster_interest])
    plot_orbits(axes,
                cluster_peri[cluster_idx], cluster_apo[cluster_idx], cluster_e[cluster_idx],
                cluster_orbit.angular_momentum()[2][0,:].to(u.kpc*u.km/u.s)[cluster_idx], cluster_r[cluster_idx],
                s=3, marker='s', INTEREST=True, name=cluster_interest)

plot_orbits(axes,
            fsr1758_orbit.pericenter(approximate=True),
            fsr1758_orbit.apocenter(approximate=True),
            fsr1758_orbit.eccentricity(approximate=True),
            fsr1758_orbit.angular_momentum()[2][0].to(u.kpc*u.km/u.s),
            fsr1758_orbit.spherical.distance[0],
            s=25, marker='*', INTEREST=True, name="FSR1758")

axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[1].set_xscale('log')
axes[0].yaxis.set_major_formatter(ScalarFormatter())
axes[0].xaxis.set_major_formatter(ScalarFormatter())
axes[1].xaxis.set_major_formatter(ScalarFormatter())

axes[0].set_ylabel("apocentre (kpc)")
axes[0].set_xlabel("pericentre (kpc)")
axes[2].set_xlabel("eccentricity")
axes[3].set_xlabel(r"$L_z$ (kpc km/s)")
axes[1].set_xlabel(r"$R_{GC}$ (kpc)")

axes[3].annotate("retrograde", (1000, 1), fontsize=4)
axes[3].annotate("prograde", (-3000, 1), fontsize=4)

axes[0].set_xlim(0.1, 30)
axes[0].set_ylim(0.8, 110)
axes[1].set_xlim(0.3, 130)
axes[2].set_xlim(0, 1)
# axes[3].set_xlim(-4500, 3500)

fig.align_labels()
plt.savefig("../fsr1758_paper/figures/orbit_comparison_new.pdf", bbox_inches='tight')
plt.show()
