from __future__ import annotations
from dataclasses import dataclass
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

@dataclass
class ResolvedTarget:
    name: str
    ra_deg: float
    dec_deg: float
    pmra_masyr: float | None
    pmdec_masyr: float | None
    plx_mas: float | None
    gaia_mag: float | None

import numpy as np
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u

def resolve_star(name: str) -> ResolvedTarget:
    """Resolve a star name via SIMBAD and return a ResolvedTarget with RA/Dec/PM/Parallax."""
    sim = Simbad()
    sim.add_votable_fields('ra(d)', 'dec(d)', 'pmra', 'pmdec', 'plx', 'flux(G)')
    res = sim.query_object(name)
    if res is None or len(res) == 0:
        raise ValueError(f"SIMBAD could not resolve: {name}")

    cols = set(res.colnames)

    # RA/Dec
    if {'ra_d', 'dec_d'} <= cols:
        ra_deg = float(res['ra_d'][0])
        dec_deg = float(res['dec_d'][0])
    else:
        c = SkyCoord(res['ra'][0], res['dec'][0], unit=(u.hourangle, u.deg))
        ra_deg, dec_deg = c.ra.deg, c.dec.deg

    # Proper motions
    pmra  = float(res['pmra'][0])  if 'pmra'  in cols and res['pmra'][0]  is not None else None
    pmdec = float(res['pmdec'][0]) if 'pmdec' in cols and res['pmdec'][0] is not None else None

    # Parallax
    plx = float(res['plx_value'][0]) if 'plx_value' in cols and res['plx_value'][0] is not None else None

    # Gaia magnitude — fully coerced to float
    g_mag = None
    if 'G' in cols and res['G'][0] is not None:
        val = res['G'][0]
        if not np.ma.is_masked(val):
            g_mag = float(np.asarray(val, dtype=float))

    return ResolvedTarget(
        name=name,
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        pmra_masyr=pmra,
        pmdec_masyr=pmdec,
        plx_mas=plx,
        gaia_mag=g_mag,
    )

