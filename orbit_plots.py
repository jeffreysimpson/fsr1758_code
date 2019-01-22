#!/usr/bin/env python

"""orbit_plots.py: Make the orbit plots for the FSR1758 paper."""

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import gala.coordinates as gc
import gala.dynamics as gd
import gala.potential as gp
from gala.units import galactic
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("paper", font_scale=0.8)


def create_error_samples(nominal_position, errors, n_samples=100):
    """Create a random sample of cluster positions using the errors."""
    dist = np.random.normal(
        nominal_position.distance.value,
        errors.distance.value,
        n_samples) * nominal_position.distance.unit

    pm_ra_cosdec = np.random.normal(
        nominal_position.pm_ra_cosdec.value,
        errors.pm_ra_cosdec.value,
        n_samples) * nominal_position.pm_ra_cosdec.unit

    pm_dec = np.random.normal(
        nominal_position.pm_dec.value,
        errors.pm_dec.value,
        n_samples) * nominal_position.pm_dec.unit

    rv = np.random.normal(
        nominal_position.radial_velocity.value,
        errors.radial_velocity.value,
        n_samples) * nominal_position.radial_velocity.unit

    ra = np.full(n_samples, nominal_position.ra.degree) * u.degree
    dec = np.full(n_samples, nominal_position.dec.degree) * u.degree

    return coord.ICRS(ra=ra, dec=dec, distance=dist,
                      pm_ra_cosdec=pm_ra_cosdec,
                      pm_dec=pm_dec, radial_velocity=rv)


def compute_orbit(icrs, dt=-0.1*u.Myr, n_steps=50000):
    """Integrate the orbit."""
    c_icrs = icrs.transform_to(gc_frame).cartesian
    object_phase_space = gd.PhaseSpacePosition(
        pos=c_icrs.xyz,
        vel=c_icrs.differentials['s'].d_xyz)
    return gp.Hamiltonian(pot).integrate_orbit(object_phase_space,
                                               dt=dt,
                                               n_steps=n_steps)


pot = gp.MilkyWayPotential()
v_sun = coord.CartesianDifferential([11, 248, 7.25]*u.km/u.s)
gc_frame = coord.Galactocentric(galcen_distance=8.2*u.kpc,
                                z_sun=25.*u.pc,
                                galcen_v_sun=v_sun)

icrs = coord.ICRS(ra=262.806 * u.deg, dec=-39.822 * u.deg,
                  distance=11.5 * u.kpc,
                  pm_ra_cosdec=-2.85 * u.mas/u.yr,
                  pm_dec=2.55 * u.mas/u.yr,
                  radial_velocity=227 * u.km/u.s)
icrs_err = coord.ICRS(ra=0*u.deg, dec=0*u.deg,
                      distance=1.*u.kpc,
                      pm_ra_cosdec=0.1*u.mas/u.yr,
                      pm_dec=0.1*u.mas/u.yr,
                      radial_velocity=1.*u.km/u.s)

icrs_samples = create_error_samples(icrs, icrs_err, n_samples=100)

fsr1758_orbit = compute_orbit(icrs)
fsr1758_orbit_errors = compute_orbit(icrs_samples)

fig, axes = plt.subplots(1, 3, figsize=(3.32*2, 2.3))
fsr1758_orbit.plot(color='C3', zorder=1000, axes=fig.axes)
fsr1758_orbit_errors.plot(color='C0', alpha=0.1, lw=0.3, axes=fig.axes)
fig.align_labels()
plt.savefig("../fsr1758_paper/figures/orbit.pdf", bbox_inches='tight')
plt.close("all")

# Now redo the orbits with lots of samples
icrs_samples = create_error_samples(icrs, icrs_err, n_samples=1000)
fsr1758_orbit_errors = compute_orbit(icrs_samples)

pers = fsr1758_orbit_errors.pericenter(approximate=True)
apos = fsr1758_orbit_errors.apocenter(approximate=True)
eccs = fsr1758_orbit_errors.eccentricity(approximate=True)

# Print the values in LaTeX format.
mc = np.percentile(pers.value, [16, 50, 84])
q = np.diff(mc)
print(f"${mc[1]:0.2f}_{{-{q[0]:0.2f}}}^{{+{q[1]:0.2f}}}$~kpc")

mc = np.percentile(apos.value, [16, 50, 84])
q = np.diff(mc)
print(f"${mc[1]:0.2f}_{{-{q[0]:0.2f}}}^{{+{q[1]:0.2f}}}$~kpc")

mc = np.percentile(eccs.value, [16, 50, 84])
q = np.diff(mc)
print(f"${mc[1]:0.2f}_{{-{q[0]:0.2f}}}^{{+{q[1]:0.2f}}}$")
