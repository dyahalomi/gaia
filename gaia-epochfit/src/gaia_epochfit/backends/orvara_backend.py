from __future__ import annotations
import numpy as np

class OrvaraBackend:
    """Thin wrapper that prepares absolute astrometry epochs for Orvara.
    Requires orvara to be installed.
    """
    def __init__(self, star_mass_msun: float, parallax_mas: float):
        try:
            import orvara  # type: ignore
        except Exception as e:
            raise ImportError("Install orvara to use Orvara backend") from e
        self.star_mass = star_mass_msun
        self.parallax_mas = parallax_mas

    def fit(self, t_obs, x_obs, phi_obs, ra_ref_deg, dec_ref_deg, n_planets: int):
        import orvara
        from orvara import drivers
        data = drivers.Data()
        data.absolute_astrometry = {
            'gaia_epochs': dict(time=np.array(t_obs), angle=np.array(phi_obs), x=np.array(x_obs),
                                ra_ref=ra_ref_deg, dec_ref=dec_ref_deg)
        }
        cfg = drivers.Config()
        cfg.n_planets = n_planets
        cfg.parallax = self.parallax_mas/1000.0  # arcsec
        cfg.mass_star = self.star_mass
        sampler = drivers.Sampler(data, cfg)
        sampler.run()
        return sampler.get_results()
