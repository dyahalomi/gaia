from __future__ import annotations
import numpy as np
import pymc as pm
import aesara.tensor as at
import exoplanet as xo

class ExoplanetBackend:
    """PyMC/Exoplanet model for absolute astrometry epochs (AL scans)."""
    def __init__(self, n_planets: int):
        self.np = n_planets

    def _model(self, t_obs, phi_obs, x_obs, x_err, star_mass, parallax_mas):
        with pm.Model() as model:
            # Priors
            pm.Normal("M_star", mu=star_mass, sigma=0.1*star_mass)
            pm.Normal("parallax", mu=parallax_mas, sigma=0.1*parallax_mas)
            # Planet parameters
            periods = pm.Lognormal("period", mu=np.log(5.0), sigma=1.0, shape=self.np)
            ecc = xo.distributions.UnitUniform("e", shape=self.np)
            omega = pm.Uniform("omega", lower=0, upper=2*np.pi, shape=self.np)
            Omega = pm.Uniform("Omega", lower=0, upper=2*np.pi, shape=self.np)
            inc = pm.Uniform("inc", lower=0, upper=np.pi, shape=self.np)
            tp = pm.Uniform("tp", lower=t_obs.min()-periods.max(), upper=t_obs.max(), shape=self.np)
            m_ratio = pm.Lognormal("q", mu=np.log(1e-3), sigma=2.0, shape=self.np)

            offsets_ra = 0.0
            offsets_dec = 0.0
            for k in range(self.np):
                orbit = xo.orbits.KeplerianOrbit(
                    period=periods[k], t0=tp[k], ecc=ecc[k], omega=omega[k],
                    Omega=Omega[k], incl=inc[k]
                )
                # Convert AU->mas via parallax: mas per AU at distance d
                mas_per_au = 1000.0 / (1.0/(parallax_mas/1000.0))
                a_star = (m_ratio[k] / (1+m_ratio[k])) * orbit.a  # AU
                ra_of_t, dec_of_t = orbit.get_relpos(t_obs)
                offsets_ra = offsets_ra + a_star * mas_per_au * ra_of_t
                offsets_dec = offsets_dec + a_star * mas_per_au * dec_of_t

            phi = at.as_tensor_variable(phi_obs)
            x_model = offsets_ra*at.sin(phi) + offsets_dec*at.cos(phi)
            pm.Normal("obs", mu=x_model, sigma=x_err, observed=x_obs)
        return model

    def fit(self, t_obs, phi_obs, x_obs, al_err_mas, star_mass, parallax_mas, draws=1000, tune=1000):
        with self._model(t_obs, phi_obs, x_obs, al_err_mas, star_mass, parallax_mas) as model:
            idata = pm.sample(draws=draws, tune=tune, chains=2, cores=2, target_accept=0.9)
        return idata
