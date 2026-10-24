"""Command line: ``python -m netsunhinged [serve|analyze] ...``"""

import argparse
import json
import sys

import numpy as np


def cmd_analyze(args):
    from . import deform, load, rigidity
    from .model import MechanicalModel, ModelOptions, preset

    rows = []
    for path in args.cif:
        st = load(path)
        if args.options:
            with open(args.options, encoding="utf-8") as fh:
                opts = ModelOptions.from_dict(json.load(fh))
        else:
            opts = preset(st, args.preset)
        m = MechanicalModel(st, opts)
        rig = rigidity.analyze(m)
        print(f"\n=== {path}  ({st.formula(st.framework)}, {st.N} atoms, {st.n_guest} guest atoms)")
        for w in st.warnings:
            print("  warning:", w)
        s = m.summary()
        print(f"  model: {s['n_bodies']} rigid bodies, {s['n_hinges']} hinges, "
              f"{s['n_ball_joints']} ball joints, {s['n_flex_angle_atoms']} pivot atoms")
        print(f"  infinitesimal: {rig['n_macroscopic']} cell-changing mechanism(s), "
              f"{rig['n_internal']} internal mechanism(s)")
        for k, md in enumerate(rig["macroscopic_modes"]):
            print(f"    mode {k + 1}: {md['kind']}: {md['description']}")
        verdict = "-"
        if not args.no_finite:
            regimes = list(deform.REGIMES) + (list(deform.EXTENDED) if args.extended else [])
            _, res, verdict, dt = deform.run_all(
                m, regimes=regimes, step=args.step, metric=args.metric,
                thr=args.thr, thr2=args.thr2, modes=rig["macroscopic_modes"])
            print(f"  finite deformation ({dt:.1f} s):")
            key = "F" if args.metric == "relative" else "E"
            for name, r in res.items():
                mx = max(p[key] for p in r["points"])
                print(f"    {name:8s} {r['classification']:14s} max {key} = {mx:9.2e}  {r['detail']}")
            print(f"  VERDICT: {verdict.upper()}")
        rows.append((path, rig["n_macroscopic"], verdict))
    if len(rows) > 1:
        print("\nsummary:")
        for path, nm, v in rows:
            print(f"  {v:14s} {nm} mode(s)  {path}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="netsunhinged", description=__doc__)
    sub = p.add_subparsers(dest="cmd")
    s = sub.add_parser("serve", help="start the interactive web app")
    s.add_argument("--port", type=int, default=8765)
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--no-browser", action="store_true")
    a = sub.add_parser("analyze", help="analyse CIF files from the command line")
    a.add_argument("cif", nargs="+")
    a.add_argument("--preset", default="auto",
                   choices=["auto", "paper", "zeolite", "rigid", "bondbend", "donor_pivots"])
    a.add_argument("--options", help="JSON file with model options (as exported by the app)")
    a.add_argument("--step", type=float, default=0.02)
    a.add_argument("--metric", default="relative", choices=["relative", "absolute"])
    a.add_argument("--thr", type=float, default=0.01)
    a.add_argument("--thr2", type=float, default=0.05)
    a.add_argument("--extended", action="store_true", help="also run fully relaxed regimes")
    a.add_argument("--no-finite", action="store_true", help="only infinitesimal analysis")
    args = p.parse_args(argv)
    if args.cmd == "analyze":
        cmd_analyze(args)
    else:
        from .server import serve
        if args.cmd is None:
            serve()
        else:
            serve(args.host, args.port, not args.no_browser)


if __name__ == "__main__":
    sys.exit(main())
