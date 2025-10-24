# gaia_epochfit/visualize.py
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import astromet

import matplotlib 
matplotlib.rc('xtick', labelsize=23) 
matplotlib.rc('ytick', labelsize=23) 
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

def plot_scan_geometry(params, sim: dict, mag: float | None = None, rng_seed: int | None = 42,
                       ax: plt.Axes | None = None, smooth_points: int = 10000):
    """
    Replicates the astromet.py tutorial visualization:
      - black: true center-of-light (COL) sky-plane track (ΔRA*cosDec vs ΔDec) in mas
      - orange: along-scan 1D measurements (segments of ±σ_AL) at each transit

    Uses per-transit quantities from `sim`: ts, phis (deg), rac_true, dec_true.
    If `mag` is None, tries sim['resolved'].gaia_mag.
    """
    ts       = np.asarray(sim["ts"], float)                   # transit times (yr)
    phis_deg = np.asarray(sim["phis"], float)                 # scan angles (deg)
    racs     = np.asarray(sim["rac_true"], float)             # true ΔRA* [mas]
    decs     = np.asarray(sim["dec_true"], float)             # true ΔDec [mas]

    # choose magnitude (sets AL noise)
    if mag is None:
        tgt = sim.get("resolved", None)
        mag = getattr(tgt, "gaia_mag", None) if tgt is not None else None
    if mag is None:
        raise ValueError("No magnitude available; pass mag= or ensure resolved.gaia_mag is set.")
    mag = float(np.asarray(mag, dtype=float))

    # AL 1σ error from astromet
    al_error = astromet.sigma_ast(mag)  # mas

    # simulated noisy AL centers at each transit (tutorial style)
    if rng_seed is not None:
        rng = np.random.default_rng(rng_seed)
    else:
        rng = np.random.default_rng()
    
    errs = al_error * rng.standard_normal(ts.size)  # mas
    radphis = np.deg2rad(phis_deg)
    obs_rac = racs + errs * np.sin(radphis)
    obs_dec = decs + errs * np.cos(radphis)

    plotts=np.linspace(np.min(ts),np.max(ts),1000)
    plotracs,plotdecs=astromet.track(plotts,params)

    # make plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 7))

    else:
        fig = ax.figure

    # draw the 1D along-scan segments (length 2*al_error) at each transit

    x1 = obs_rac - al_error * np.sin(radphis)
    x2 = obs_rac + al_error * np.sin(radphis)
    y1 = obs_dec - al_error * np.cos(radphis)
    y2 = obs_dec + al_error * np.cos(radphis)

    for i in range(ts.size):
        ax.plot([x1[i], x2[i]], [y1[i], y2[i]], c="orange", lw=0.8, alpha=1.)

    # true track in black
    ax.plot(plotracs, plotdecs, c="k", lw=1, zorder=-1000)

    ax.set_xlabel(r"$RA\cos(\mathrm{Dec})$ [mas]", fontsize=27)
    ax.set_ylabel(r"$Dec$ [mas]", fontsize=27)
    ax.set_title("True track (black) and 1D Gaia AL segments (orange)", fontsize=27)

    fig.tight_layout()

    return fig, ax, radphis
