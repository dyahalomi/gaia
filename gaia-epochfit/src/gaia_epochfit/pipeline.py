from __future__ import annotations
import numpy as np
import astromet
from .resolve import resolve_star
from .utils import StellarParams, PlanetParams, reflex_offsets_mas
from . import gaia_epochs

class SimulationAndFit:
    def __init__(self, backend: str = 'exoplanet'):
        self.backend_name = backend.lower()

    def _backend(self, n_planets: int, star_mass_msun: float, parallax_mas: float):
        if self.backend_name == 'exoplanet':
            from .backends.exoplanet_backend import ExoplanetBackend
            return ExoplanetBackend(n_planets)
        elif self.backend_name == 'orvara':
            from .backends.orvara_backend import OrvaraBackend
            return OrvaraBackend(star_mass_msun, parallax_mas)
        elif self.backend_name == 'orbitize':
            from .backends.orbitize_backend import OrbitizeBackend
            return OrbitizeBackend()
        elif self.backend_name == 'octofitter':
            from .backends.octofitter_backend import OctofitterBackend
            return OctofitterBackend()
        else:
            raise ValueError(f"Unknown backend: {self.backend_name}")

    def simulate(self, target_name: str, star_mass_msun: float,
             planets: list[PlanetParams],
             epoch_range=(2014.5, 2025.5),
             parallax_mas: float | None = None,
             mag: float | None = None):

        # Resolve star for ra/dec/parallax
        tgt = resolve_star(target_name)

        # --- Parallax handling ---
        if tgt.plx_mas is None:
            if parallax_mas is not None:
                tgt.plx_mas = parallax_mas
            else:
                raise ValueError(
                    "No parallax from SIMBAD; please pass parallax explicitly as parallax_mas="
                )

        # --- Magnitude handling ---
        if tgt.gaia_mag is None:
            if mag is not None:
                tgt.gaia_mag = mag
            else:
                raise ValueError(
                    "No Gaia mag from SIMBAD; please pass magnitude explicitly as mag ="
                )


        # Get scanning-law epoch times and angles
        ts, phis = gaia_epochs.nominal_times_and_angles(tgt.ra_deg, tgt.dec_deg, tstart=epoch_range[0], tend=epoch_range[1])
        # Restrict to epoch_range
        mask = (ts >= epoch_range[0]) & (ts <= epoch_range[1])
        ts, phis = ts[mask], phis[mask]
        # Compose planetary reflex offsets
        dRAc_mas, dDec_mas = reflex_offsets_mas(ts, star_mass_msun, tgt.plx_mas, planets)
        # Build absolute track by adding 5p astrometric terms around reference epoch
        p = astromet.params()
        p.ra = tgt.ra_deg
        p.dec = tgt.dec_deg
        p.pmrac = tgt.pmra_masyr 
        p.pmdec = tgt.pmdec_masyr
        p.parallax = tgt.plx_mas
        p.G = tgt.gaia_mag
        base_rac, base_dec = astromet.track(ts, p)  # mas offsets for center-of-mass 5p model
        rac_true = base_rac + dRAc_mas
        dec_true = base_dec + dDec_mas
        # Simulate Gaia along-scan measurements
        t_obs, x_obs, phi_obs, rac_obs, dec_obs = gaia_epochs.simulate_gaia_al_obs(ts, phis, rac_true, dec_true, tgt.gaia_mag)
        return dict(resolved=tgt, ts=ts, phis=phis, rac_true=rac_true, dec_true=dec_true,
                    t_obs=t_obs, x_obs=x_obs, phi_obs=phi_obs, rac_obs=rac_obs, dec_obs=dec_obs)

    def fit(self, sim, star_mass_msun: float, backend_kwargs: dict | None = None):
        backend_kwargs = backend_kwargs or {}
        n_planets = backend_kwargs.pop('n_planets', None)
        if n_planets is None:
            n_planets = backend_kwargs.get('nplanets', 1)
        be = self._backend(n_planets, star_mass_msun, sim['resolved'].plx_mas)
        if self.backend_name == 'exoplanet':
            al_err = astromet.sigma_ast(backend_kwargs.get('mag', 12.0))
            idata = be.fit(sim['t_obs'], sim['phi_obs'], sim['x_obs'], al_err,
                           star_mass_msun, sim['resolved'].plx_mas,
                           draws=backend_kwargs.get('draws', 1000),
                           tune=backend_kwargs.get('tune', 1000))
            return idata
        elif self.backend_name == 'orvara':
            return be.fit(sim['t_obs'], sim['x_obs'], sim['phi_obs'],
                          sim['resolved'].ra_deg, sim['resolved'].dec_deg, n_planets)
        else:
            raise NotImplementedError(f"Backend {self.backend_name} not yet wired for fit()")
