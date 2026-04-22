#!/usr/bin/env python3
"""Gaudi build/run CLI.

Usage:
    python gaudi.py test [--debug]
    python gaudi.py run <target> [--debug] [--shift-fraction F] [-- EXTRA ...]
    python gaudi.py build <target> [--debug]
    python gaudi.py configure <profile>
    python gaudi.py list
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ── Build profiles ──────────────────────────────────────────────────────────

PROFILES = {
    "test": {
        "dir": "build_test",
        "cmake_args": ["-DGAUDI_WITH_GL=OFF", "-DBUILD_TESTING=ON"],
    },
    "gl": {
        "dir": "build_gl_viewers",
        "cmake_args": ["-DGAUDI_WITH_GL=ON"],
    },
    "asan": {
        "dir": "build_asan",
        "cmake_args": ["-DGAUDI_WITH_GL=ON", "-DCMAKE_BUILD_TYPE=asan"],
    },
}

GENERATOR = "Visual Studio 17 2022"

# ── Target resolution ───────────────────────────────────────────────────────

TEST_ALIASES = {"tests", "test", "gaudi_tests"}


def _enabled_gl_targets():
    """Parse CMakeLists.txt for uncommented add_subdirectory(projects/X)."""
    cml = ROOT / "CMakeLists.txt"
    targets = []
    for line in cml.read_text().splitlines():
        m = re.match(r"\s*add_subdirectory\(projects/(\S+)\)", line)
        if m:
            targets.append(m.group(1))
    return targets


def _all_gl_targets():
    """Every projects/ dir that has a CMakeLists.txt."""
    return sorted(
        d.name
        for d in (ROOT / "projects").iterdir()
        if d.is_dir() and (d / "CMakeLists.txt").exists()
    )


def resolve_target(name):
    """Return (profile_key, cmake_target, exe_path_template)."""
    if name in TEST_ALIASES:
        return "test", "gaudi_tests", "tests/{config}/gaudi_tests.exe"

    gl = _enabled_gl_targets()
    if name in gl:
        return "gl", name, f"projects/{name}/{{config}}/{name}.exe"

    all_gl = _all_gl_targets()
    if name in all_gl:
        print(
            f"Warning: '{name}' exists in projects/ but is commented out in "
            f"CMakeLists.txt. It may fail to build.",
            file=sys.stderr,
        )
        return "gl", name, f"projects/{name}/{{config}}/{name}.exe"

    sys.exit(f"Error: unknown target '{name}'. Run `python gaudi.py list`.")


# ── Actions ─────────────────────────────────────────────────────────────────


def configure(profile_key):
    prof = PROFILES[profile_key]
    build_dir = ROOT / prof["dir"]
    cmd = [
        "cmake",
        "-S", str(ROOT),
        "-B", str(build_dir),
        "-G", GENERATOR,
        *prof["cmake_args"],
    ]
    print(f">> {' '.join(cmd)}")
    return subprocess.run(cmd).returncode


def ensure_configured(profile_key):
    prof = PROFILES[profile_key]
    build_dir = ROOT / prof["dir"]
    if not (build_dir / "CMakeCache.txt").exists():
        print(f"Build directory '{prof['dir']}' not configured. Running cmake...")
        rc = configure(profile_key)
        if rc != 0:
            sys.exit(f"cmake configure failed (exit {rc})")


def build(profile_key, target, config):
    ensure_configured(profile_key)
    build_dir = ROOT / PROFILES[profile_key]["dir"]
    cmd = [
        "cmake",
        "--build", str(build_dir),
        "--target", target,
        "--config", config,
    ]
    print(f">> {' '.join(cmd)}")
    return subprocess.run(cmd).returncode


def run_exe(profile_key, exe_template, config, extra_args=None):
    build_dir = ROOT / PROFILES[profile_key]["dir"]
    exe = build_dir / exe_template.format(config=config)
    if not exe.exists():
        sys.exit(f"Error: executable not found at {exe}")
    argv = [str(exe)]
    if extra_args:
        argv.extend(extra_args)
    print(f">> {' '.join(argv)}")
    return subprocess.run(argv).returncode


# ── CLI ─────────────────────────────────────────────────────────────────────


def cmd_test(args):
    config = "Debug" if args.debug else "Release"
    rc = build("test", "gaudi_tests", config)
    if rc != 0:
        sys.exit(rc)
    sys.exit(run_exe("test", "tests/{config}/gaudi_tests.exe", config))


def cmd_run(args):
    config = "Debug" if args.debug else "Release"
    profile, target, exe_tpl = resolve_target(args.target)
    rc = build(profile, target, config)
    if rc != 0:
        sys.exit(rc)
    forward = list(args.forward or [])
    if getattr(args, "shift_fraction", None) is not None:
        forward = ["--shift-fraction", str(args.shift_fraction)] + forward
    sys.exit(run_exe(profile, exe_tpl, config, forward if forward else None))


def cmd_build(args):
    config = "Debug" if args.debug else "Release"
    profile, target, _ = resolve_target(args.target)
    sys.exit(build(profile, target, config))


def cmd_configure(args):
    key = args.profile
    if key not in PROFILES:
        sys.exit(f"Unknown profile '{key}'. Choose from: {', '.join(PROFILES)}")
    sys.exit(configure(key))


def cmd_list(_args):
    enabled = _enabled_gl_targets()
    all_gl = _all_gl_targets()
    disabled = [t for t in all_gl if t not in enabled]

    print("Headless targets:")
    print(f"  tests          (alias: test, gaudi_tests)  [build_test]")
    print()
    print("GL viewer targets (enabled):                  [build_gl_viewers]")
    for t in enabled:
        print(f"  {t}")
    if disabled:
        print()
        print("GL viewer targets (commented out in CMakeLists.txt):")
        for t in disabled:
            print(f"  {t}  (disabled)")
    print()
    print(f"Build profiles: {', '.join(PROFILES)}")


def main():
    parser = argparse.ArgumentParser(
        prog="gaudi",
        description="Build and run Gaudi targets.",
    )
    sub = parser.add_subparsers(dest="command")

    p_test = sub.add_parser("test", help="Build and run headless tests")
    p_test.add_argument("--debug", action="store_true", help="Debug config")

    p_run = sub.add_parser(
        "run",
        help="Build and run a target",
        epilog=(
            "Example:  python gaudi.py run spectral_modes_test --shift-fraction 0.35\n"
            "          (F=0 low spectrum, F=1 high, 0<F<1 interior via shift-invert)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_run.add_argument("target", help="Target name (e.g. aabb_test)")
    p_run.add_argument("--debug", action="store_true", help="Debug config")
    p_run.add_argument(
        "--shift-fraction",
        type=float,
        default=None,
        metavar="F",
        help=(
            "spectral_modes_test: pass --shift-fraction F clamped to [0,1]: 0=low "
            "eigenmodes, 1=high |lambda|, in-between=shift-invert sigma sweep"
        ),
    )
    p_run.add_argument(
        "forward",
        nargs=argparse.REMAINDER,
        help="Extra args for the executable (often: -- --mid --modes 16)",
    )

    p_build = sub.add_parser("build", help="Build a target (no run)")
    p_build.add_argument("target", help="Target name")
    p_build.add_argument("--debug", action="store_true", help="Debug config")

    p_cfg = sub.add_parser("configure", help="(Re)generate cmake for a profile")
    p_cfg.add_argument("profile", help="Profile: test, gl, asan")

    sub.add_parser("list", help="Show available targets")

    args = parser.parse_args()

    dispatch = {
        "test": cmd_test,
        "run": cmd_run,
        "build": cmd_build,
        "configure": cmd_configure,
        "list": cmd_list,
    }

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    dispatch[args.command](args)


if __name__ == "__main__":
    main()
