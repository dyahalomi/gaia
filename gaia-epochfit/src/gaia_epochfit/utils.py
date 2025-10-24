from __future__ import annotations
import numpy as np
from dataclasses import dataclass

JYEAR = 365.25

@dataclass
class StellarParams:
    ra_deg: float
    dec_deg: float
    pmra_masyr: float
    pmdec_masyr: float
    parallax_mas: float
    epoch_jyr: float = 2016.0

@dataclass
class PlanetParams:
    period_yr: float
    a_au: float
    e: float
    i_deg: float
    omega_deg: float  # argument of periastron of planet
    Omega_deg: float  # longitude of ascending node
    tp_jyr: float     # time of periastron passage
    m_jup: float      # planet mass; used to set reflex amplitude

# Basic astro constants
M_SUN = 1.0  # in Msun (scale for host input)
M_JUP_TO_MSUN = 1.0/1047.3486
AU_PER_PC = 206265.0
MAS_PER_AS = 1000.0

def thiele_innes_constants(a_as, i, omega, Omega):
    # All angles in radians; a_as in arcsec
    cO, sO = np.cos(Omega), np.sin(Omega)
    co, so = np.cos(omega), np.sin(omega)
    ci, si = np.cos(i), np.sin(i)
    A = a_as*(cO*co - sO*so*ci)
    B = a_as*(sO*co + cO*so*ci)
    F = a_as*(-cO*so - sO*co*ci)
    G = a_as*(-sO*so + cO*co*ci)
    return A, B, F, G

def kepler_E(M, e, tol=1e-10, max_iter=50):
    E = M.copy()
    for _ in range(max_iter):
        dE = (M - (E - e*np.sin(E))) / (1 - e*np.cos(E))
        E = E + dE
        if np.all(np.abs(dE) < tol):
            break
    return E

def reflex_offsets_mas(times_jyr: np.ndarray,
                        star_mass_msun: float,
                        star_parallax_mas: float,
                        planets: list[PlanetParams]) -> tuple[np.ndarray, np.ndarray]:
    """Compute total sky-plane reflex motion of the *star* in mas, summed over planets.
    Returns (dRAcosDec_mas, dDec_mas) arrays matching times.
    """
    dRAc = np.zeros_like(times_jyr, dtype=float)
    dDec = np.zeros_like(times_jyr, dtype=float)

    d_pc = AU_PER_PC / (star_parallax_mas / MAS_PER_AS)  # distance [pc]

    for p in planets:
        n = 2*np.pi/p.period_yr
        M = n*(times_jyr - p.tp_jyr)
        M = np.mod(M, 2*np.pi)
        E = kepler_E(M, p.e)
        # Relative orbit in orbital plane (a=1 in arcsec later scaled)
        x = np.cos(E) - p.e
        y = np.sqrt(1 - p.e**2) * np.sin(E)
        # Semimajor axis of star's reflex in AU
        a_star_AU = p.a_au * (p.m_jup*M_JUP_TO_MSUN) / (star_mass_msun + p.m_jup*M_JUP_TO_MSUN)
        # Convert AU to arcsec at distance d_pc
        a_star_as = a_star_AU / d_pc
        # Project to sky with Thiele–Innes
        i = np.deg2rad(p.i_deg)
        omega = np.deg2rad(p.omega_deg)
        Omega = np.deg2rad(p.Omega_deg)
        A,B,F,G = thiele_innes_constants(a_star_as, i, omega, Omega)
        X = A*x + F*y
        Y = B*x + G*y
        dRAc += X*MAS_PER_AS
        dDec += Y*MAS_PER_AS

    return dRAc, dDec
