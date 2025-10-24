from gaia_epochfit.pipeline import SimulationAndFit
from gaia_epochfit.utils import PlanetParams

pipe = SimulationAndFit(backend='exoplanet')
sim = pipe.simulate(
    target_name='Beta Pic', star_mass_msun=1.8, mag=8.0,
    planets=[PlanetParams(period_yr=20, a_au=9.0, e=0.05, i_deg=88, omega_deg=90, Omega_deg=30, tp_jyr=2019.5, m_jup=10.0)],
)
idata = pipe.fit(sim, star_mass_msun=1.8, backend_kwargs=dict(n_planets=1, mag=8.0, draws=1000, tune=1000))
print(idata)
