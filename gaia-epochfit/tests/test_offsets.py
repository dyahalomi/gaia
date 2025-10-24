import numpy as np
from gaia_epochfit.utils import PlanetParams, reflex_offsets_mas

def test_reflex_zero_planet():
    t = np.linspace(2014, 2016, 10)
    dRAc, dDec = reflex_offsets_mas(t, 1.0, 10.0, [])
    assert np.allclose(dRAc, 0)
    assert np.allclose(dDec, 0)

def test_reflex_periodicity():
    t = np.linspace(2014, 2016, 200)
    p = PlanetParams(period_yr=1.0, a_au=1.0, e=0.0, i_deg=60, omega_deg=45, Omega_deg=120, tp_jyr=2014.0, m_jup=1.0)
    dRAc, dDec = reflex_offsets_mas(t, 1.0, 50.0, [p])
    assert np.isfinite(dRAc).all() and np.isfinite(dDec).all()
    assert abs(np.mean(dRAc)) < 1e-2
    assert abs(np.mean(dDec)) < 1e-2
