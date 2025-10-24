# gaia-epochfit

Simulate Gaia-like epoch astrometry and fit N-planet orbits. See `examples/quickstart.py` and `notebooks/Tutorial.ipynb` for end-to-end demos.

## Install (dev)
```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

## Minimal run
```bash
gepochfit simulate-and-fit "Beta Pic"   --nplanets 1   --pl period=20,a=9.0,e=0.05,i=88,omega=90,Omega=30,tp=2019.5,mjup=10   --mag 8.0 --duration 2014,2025   --backend exoplanet --draws 1000 --tune 1000
```
