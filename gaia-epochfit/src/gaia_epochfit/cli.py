import argparse
from .pipeline import SimulationAndFit
from .utils import PlanetParams

def parse_planet_arg(s: str) -> PlanetParams:
    # period=,a=,e=,i=,omega=,Omega=,tp=,mjup=
    kv = dict(x.split('=') for x in s.split(','))
    return PlanetParams(
        period_yr=float(kv['period']), a_au=float(kv['a']), e=float(kv['e']),
        i_deg=float(kv['i']), omega_deg=float(kv['omega']), Omega_deg=float(kv['Omega']),
        tp_jyr=float(kv['tp']), m_jup=float(kv['mjup'])
    )

def main():
    ap = argparse.ArgumentParser(prog='gepochfit')
    sub = ap.add_subparsers(dest='cmd', required=True)

    sp = sub.add_parser('simulate-and-fit')
    sp.add_argument('star_name')
    sp.add_argument('--nplanets', type=int, default=1)
    sp.add_argument('--pl', action='append', required=True,
                    help='planet spec: period=,a=,e=,i=,omega=,Omega=,tp=,mjup= (repeat)')
    sp.add_argument('--mag', type=float, default=12.0)
    sp.add_argument('--mass', type=float, default=1.0, help='stellar mass [Msun]')
    sp.add_argument('--backend', choices=['exoplanet','orvara','orbitize','octofitter'], default='exoplanet')
    sp.add_argument('--draws', type=int, default=1000)
    sp.add_argument('--tune', type=int, default=1000)

    args = ap.parse_args()

    planets = [parse_planet_arg(p) for p in args.pl]
    pipe = SimulationAndFit(backend=args.backend)
    sim = pipe.simulate(args.star_name, star_mass_msun=args.mass, mag=args.mag, planets=planets)
    res = pipe.fit(sim, star_mass_msun=args.mass, backend_kwargs=dict(n_planets=args.nplanets, mag=args.mag, draws=args.draws, tune=args.tune))
    print("FIT COMPLETE:")
    print(res)

if __name__ == '__main__':
    main()
