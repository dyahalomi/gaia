# gaia_epochfit/gaia_epochs.py
from __future__ import annotations
import numpy as np
import astromet
import gaiascanlaw  # imported globally, as requested


def nominal_times_and_angles(
    ra_deg: float,
    dec_deg: float,
    *,
    tstart: float | None = None,
    tend: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (times_jyr, angles_rad) for Gaia-like visits using gaiascanlaw.
    - times are decimal Julian years
    - angles are radians
    """
    ts0 = tstart if tstart is not None else getattr(gaiascanlaw, "tstart", 2014.6)
    te0 = tend   if tend   is not None else getattr(gaiascanlaw, "tdr4",   2020.1)
    times, angles = gaiascanlaw.scanlaw(ra_deg, dec_deg, tstart=ts0, tend=te0)
    return np.asarray(times, float), np.asarray(angles, float)


def simulate_gaia_al_obs(
    ts_jyr: np.ndarray,
    phis_rad: np.ndarray,
    racs_true_mas: np.ndarray,
    decs_true_mas: np.ndarray,
    mag: float,
):
    """
    Simulate along-scan CCD-bundle measurements using astromet’s error model,
    matching the astromet.py tutorial (angles in degrees for mock_obs).
    """
    al_err = astromet.sigma_ast(mag)
    phis_deg = np.degrees(phis_rad)  # astromet.mock_obs expects degrees
    return astromet.mock_obs(ts_jyr, phis_deg, racs_true_mas, decs_true_mas, err=al_err)


def gaia_like_fit(
    t_obs: np.ndarray,
    x_obs: np.ndarray,
    phi_obs_deg: np.ndarray,
    mag: float,
    ra_ref_deg: float,
    dec_ref_deg: float,
):
    """
    Convenience wrapper around astromet.gaia_fit (Gaia-space fit).
    """
    al_err = astromet.sigma_ast(mag)
    return astromet.gaia_fit(t_obs, x_obs, phi_obs_deg, al_err, ra_ref_deg, dec_ref_deg)
