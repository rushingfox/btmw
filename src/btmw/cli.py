"""``btmw`` command-line interface.

Each subcommand corresponds to one or more paper figures and delegates
to the relevant module under ``btmw.figures``.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from . import __version__


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _cmd_projection_map(args: argparse.Namespace) -> int:
    from .figures import projection_map
    out = projection_map.plot(
        sim_label=args.sim,
        width_kpc=args.width_kpc,
        output=args.output,
        use_tex=not args.no_tex,
    )
    print(f"wrote {out}")
    return 0


def _cmd_radial_density(args: argparse.Namespace) -> int:
    from .figures import radial_density
    # use_tex=None → auto-select from backend (HBT→True, VR→False).
    # --no-tex forces use_tex=False regardless of backend.
    use_tex = False if args.no_tex else None
    out = radial_density.plot(
        host=args.host,
        backend=args.center_from,
        refresh=args.refresh,
        output=args.output,
        use_tex=use_tex,
    )
    print(f"wrote {out}")
    return 0


def _cmd_hmf(args: argparse.Namespace) -> int:
    from .figures import hmf
    if args.resolution_study:
        out = hmf.plot_resolution_study(refresh=args.refresh, output=args.output, use_tex=not args.no_tex)
        print(f"wrote {out}")
        return 0
    if args.compare_vr:
        out = hmf.plot_compare_vr(refresh=args.refresh, use_tex=not args.no_tex)
        print(f"wrote {out}")
    else:
        out = hmf.plot(refresh=args.refresh, output=args.output, use_tex=not args.no_tex)
        print(f"wrote {out}")
    return 0


def _cmd_hvf(args: argparse.Namespace) -> int:
    from .figures import hvf
    if args.resolution_study:
        out = hvf.plot_resolution_study(refresh=args.refresh, output=args.output, use_tex=not args.no_tex)
        print(f"wrote {out}")
        return 0
    if args.compare_vr:
        out = hvf.plot_compare_vr(refresh=args.refresh, use_tex=not args.no_tex)
        print(f"wrote {out}")
    else:
        out = hvf.plot(refresh=args.refresh, output=args.output, use_tex=not args.no_tex)
        print(f"wrote {out}")
    return 0


def _cmd_hrf(args: argparse.Namespace) -> int:
    from .figures import hrf
    bins = args.bin if args.bin else list(hrf.MASS_BINS)
    if args.resolution_study:
        for mass_bin in bins:
            out = hrf.plot_resolution_study(
                mass_bin=mass_bin,
                refresh=args.refresh,
                use_tex=not args.no_tex,
            )
            print(f"wrote {out}")
        return 0
    if args.compare_vr:
        for mass_bin in bins:
            out = hrf.plot_compare_vr(
                mass_bin=mass_bin,
                refresh=args.refresh,
                use_tex=not args.no_tex,
            )
            print(f"wrote {out}")
    else:
        for mass_bin in bins:
            out = hrf.plot(
                mass_bin=mass_bin,
                refresh=args.refresh,
                output=args.output,
                use_tex=not args.no_tex,
            )
            print(f"wrote {out}")
    return 0


def _cmd_cmf(args: argparse.Namespace) -> int:
    from .figures import cmf
    out = cmf.plot(refresh=args.refresh, output=args.output, use_tex=not args.no_tex)
    print(f"wrote {out}")
    return 0


def _cmd_rvsv(args: argparse.Namespace) -> int:
    from .figures import rvsv
    if args.resolution_study:
        out = rvsv.plot_resolution_study(refresh=args.refresh, output=args.output, use_tex=not args.no_tex)
        print(f"wrote {out}")
        return 0
    hosts = [args.host] if args.host else ["m12i", "m12f"]
    for host in hosts:
        out = rvsv.plot(host=host, refresh=args.refresh, use_tex=not args.no_tex)
        print(f"wrote {out}")
    return 0


def _cmd_mvsv(args: argparse.Namespace) -> int:
    from .figures import mvsv
    hosts = [args.host] if args.host else ["m12i", "m12f"]
    for host in hosts:
        out = mvsv.plot(
            host=host, refresh=args.refresh, output=args.output,
            use_tex=not args.no_tex,
        )
        print(f"wrote {out}")
    return 0


def _cmd_list_sims(args: argparse.Namespace) -> int:
    """Print every simulation declared in ``configs/simulations.yaml``."""
    from .config import load_simulations

    cfg = load_simulations()
    sim_root = cfg.working_root if cfg.working_root is not None else "(not set)"
    print(f"BTMW_SIM_ROOT: {sim_root}")
    print(f"{'label':<26}  {'host':<5}  {'model':<8}  {'res':<8}  snap   has_vr")
    for label, sim in cfg.simulations.items():
        print(
            f"{label:<26}  {sim.host:<5}  {sim.model:<8}  "
            f"{sim.resolution:<8}  {sim.snapshot_num:>4}   {sim.has_vr}"
        )
    return 0


# ---------------------------------------------------------------------------
# Parser construction
# ---------------------------------------------------------------------------

def _add_common_options(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--no-tex",
        action="store_true",
        help="Disable LaTeX rendering (for environments without TeX).",
    )
    p.add_argument(
        "--refresh",
        action="store_true",
        help="Re-extract cache from raw data instead of using data/cache/.",
    )
    p.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output figure path. Defaults to the paper filename.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="btmw",
        description=(
            "Analysis code for Blue-Tilt Milky-Way simulations "
            "(arXiv:2412.16072)."
        ),
    )
    parser.add_argument("--version", action="version", version=f"btmw {__version__}")

    sub = parser.add_subparsers(dest="command", required=True, metavar="<command>")

    # ----- list-sims (utility) -----
    p = sub.add_parser("list-sims", help="List the simulations declared in configs/.")
    p.set_defaults(func=_cmd_list_sims)

    # ----- fig 4 -----
    p = sub.add_parser("projection-map", help="Fig 4: matter projection maps.")
    p.add_argument("--sim", required=True, help="Simulation label, e.g. m12i_cdmo.")
    p.add_argument("--width-kpc", type=float, default=200.0)
    p.add_argument(
        "--no-tex",
        action="store_true",
        help="Disable LaTeX rendering (for environments without TeX).",
    )
    p.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output figure path. Defaults to the paper filename.",
    )
    p.set_defaults(func=_cmd_projection_map)

    # ----- fig 5 -----
    p = sub.add_parser(
        "radial-density",
        help="Fig 5: main-halo radial density profile (m12i / m12f).",
    )
    p.add_argument("--host", choices=["m12i", "m12f"], required=True)
    p.add_argument(
        "--center-from",
        choices=["hbt", "vr"],
        default=None,
        help=(
            "Which data source to use: 'hbt' (HBT particle-based, v1 original) or "
            "'vr' (VelocIRaptor profile-based, PRD version). "
            "Default: 'hbt' for m12i, 'vr' for m12f (paper default)."
        ),
    )
    _add_common_options(p)
    p.set_defaults(func=_cmd_radial_density)

    # ----- fig 6 / 11 / 15 -----
    p = sub.add_parser("hmf", help="Fig 6 (HMF), Fig 11 (HMF VR-SOAP), Fig 15 (HMF resolution).")
    p.add_argument("--compare-vr", action="store_true", help="Fig 11 (PRD version).")
    p.add_argument("--resolution-study", action="store_true", help="Fig 15.")
    _add_common_options(p)
    p.set_defaults(func=_cmd_hmf)

    # ----- fig 7 / 12 / 16 -----
    p = sub.add_parser("hvf", help="Fig 7 (HVF), Fig 12 (HVF VR-SOAP), Fig 16 (HVF resolution).")
    p.add_argument("--compare-vr", action="store_true")
    p.add_argument("--resolution-study", action="store_true")
    _add_common_options(p)
    p.set_defaults(func=_cmd_hvf)

    # ----- fig 8 / 13 / 18 -----
    p = sub.add_parser("hrf", help="Fig 8 / 13 / 18: scaled subhalo radial number-density profile.")
    p.add_argument(
        "--bin",
        action="append",
        type=int,
        default=None,
        choices=[6, 7, 8, 9],
        help="Mass bin (paper uses 6, 7, 8, 9). May be repeated.",
    )
    p.add_argument("--compare-vr", action="store_true")
    p.add_argument("--resolution-study", action="store_true")
    _add_common_options(p)
    p.set_defaults(func=_cmd_hrf)

    # ----- fig 9 -----
    p = sub.add_parser("cmf", help="Fig 9: scaled cumulative subhalo mass function (CMF).")
    _add_common_options(p)
    p.set_defaults(func=_cmd_cmf)

    # ----- fig 10 / 17 -----
    p = sub.add_parser("rvsv", help="Fig 10 (RvsV), Fig 17 (RvsV resolution).")
    p.add_argument("--host", choices=["m12i", "m12f"], default=None)
    p.add_argument("--resolution-study", action="store_true")
    _add_common_options(p)
    p.set_defaults(func=_cmd_rvsv)

    # ----- fig 14 -----
    p = sub.add_parser("mvsv", help="Fig 14: Mvir vs Vmax (VELOCIraptor, both hosts by default).")
    p.add_argument("--host", choices=["m12i", "m12f"], default=None)
    _add_common_options(p)
    p.set_defaults(func=_cmd_mvsv)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
