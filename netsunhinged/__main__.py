"""Command line: ``python -m netsunhinged [serve|analyze] ...``"""

import argparse
import json
import sys

import numpy as np


def cmd_analyze(args):
    from . import deform, load, rigidity
    from .cif import cif_block_names
    from .model import MechanicalModel, ModelOptions, preset

    # every structure of every file (a CIF can hold several data blocks)
    jobs = []
    for path in args.cif:
        with open(path, encoding="utf-8", errors="replace") as fh:
            names = cif_block_names(fh.read())
        jobs += [(path if len(names) == 1 else f"{path} [{name}]", path, k) for k, name in enumerate(names)]
    if not jobs:
        sys.exit("no structure (data block with _atom_site_fract_x) found")
    rows = []
    for label, path, block in jobs:
        st = load(path, block)
        if args.options:
            with open(args.options, encoding="utf-8") as fh:
                opts = ModelOptions.from_dict(json.load(fh))
        else:
            opts = preset(st, args.preset)
        m = MechanicalModel(st, opts)
        rig = rigidity.analyze(m)
        print(f"\n=== {label}  ({st.formula(st.framework)}, {st.N} atoms, {st.n_guest} guest atoms)")
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
                m, regimes=regimes, step=args.step, thr=args.thr, thr2=args.thr2, modes=rig["macroscopic_modes"])
            print(f"  finite deformation ({dt:.1f} s):")
            for name, r in res.items():
                mf = max(p["F"] for p in r["points"])
                me = max(p["E"] for p in r["points"])
                print(f"    {name:8s} {r['classification']:14s} max F = {mf:9.2e}  max E = {me:9.2e} A^2"
                      f"  {r['detail']}")
            print(f"  VERDICT: {verdict.upper()}")
        rows.append((label, rig["n_macroscopic"], verdict))
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
    a.add_argument("--thr", type=float, default=0.01, help="flexible threshold on F")
    a.add_argument("--thr2", type=float, default=0.05, help="semi-flexible threshold on F")
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
