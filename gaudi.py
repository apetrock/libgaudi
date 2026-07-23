#!/usr/bin/env python3
"""Gaudi build/run CLI.

Usage:
    python gaudi.py test [--debug] [--asan]
    python gaudi.py run <target> [--debug] [--asan] [--no-record] [--record-path PATH]
                          [--shift-fraction F] [-- EXTRA ...]
    python gaudi.py build <target> [--debug] [--asan]
    python gaudi.py configure <profile> [--debug] [--asan]
    python gaudi.py list

Vermeer runs record to dump/<target>.mp4 by default (requires ffmpeg).
Use --no-record to disable, or --record-path to override the output stem.
"""

import argparse
import os
import re
import shutil
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
    "vermeer": {
        "dir": "build_vermeer",
        "cmake_args": ["-DGAUDI_WITH_VERMEER=ON", "-DBUILD_TESTING=ON"],
    },
    "asan": {
        "dir": "build_asan",
        "cmake_args": ["-DGAUDI_WITH_GL=ON", "-DCMAKE_BUILD_TYPE=asan"],
    },
}

# ── Target resolution ───────────────────────────────────────────────────────

TEST_ALIASES = {"tests", "test", "gaudi_tests"}
_SUBDIR_RE = re.compile(r"^\s*add_subdirectory\(projects/(\S+)\)")


def _cmake_generator():
    if sys.platform == "win32":
        return "Visual Studio 17 2022"
    if shutil.which("ninja"):
        return "Ninja"
    return "Unix Makefiles"


def _default_multi_config():
    return sys.platform == "win32"


def _read_cache_value(build_dir, key):
    cache = build_dir / "CMakeCache.txt"
    if not cache.exists():
        return None
    prefix = f"{key}:"
    for line in cache.read_text().splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return None


def is_multi_config(build_dir):
    gen = _read_cache_value(build_dir, "CMAKE_GENERATOR")
    if gen:
        return "Visual Studio" in gen or "Xcode" in gen
    return _default_multi_config()


def build_config(debug):
    return "Debug" if debug else "Release"


def profile_dir(profile_key, asan=False):
    base = PROFILES[profile_key]["dir"]
    return f"{base}_asan" if asan else base


def profile_cmake_args(profile_key, asan=False):
    args = list(PROFILES[profile_key]["cmake_args"])
    if asan and not any(a.startswith("-DCMAKE_BUILD_TYPE=") for a in args):
        args.append("-DCMAKE_BUILD_TYPE=asan")
    return args


def profile_build_type(profile_key, asan=False):
    for arg in profile_cmake_args(profile_key, asan):
        if arg.startswith("-DCMAKE_BUILD_TYPE="):
            return arg.split("=", 1)[1]
    return None


def effective_config(profile_key, debug, asan=False):
    if asan:
        return "asan"
    pinned = profile_build_type(profile_key, asan=False)
    if pinned:
        return pinned
    return build_config(debug)


def _enabled_targets_by_profile():
    """Parse CMakeLists.txt for enabled GL and Vermeer project targets."""
    gl, vermeer = [], []
    block = None
    for line in (ROOT / "CMakeLists.txt").read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if stripped.startswith("if(GAUDI_WITH_GL)"):
            block = "gl"
            continue
        if stripped.startswith("if(GAUDI_WITH_VERMEER)"):
            block = "vermeer"
            continue
        if block and (stripped.startswith("else()") or stripped.startswith("endif()")):
            block = None
            continue
        m = _SUBDIR_RE.match(line)
        if not m:
            continue
        name = m.group(1)
        if block == "gl":
            gl.append(name)
        elif block == "vermeer":
            vermeer.append(name)
    return gl, vermeer


def _all_project_targets():
    """Every projects/ dir that has a CMakeLists.txt."""
    return sorted(
        d.name
        for d in (ROOT / "projects").iterdir()
        if d.is_dir() and (d / "CMakeLists.txt").exists()
    )


def exe_path(profile_key, target, config, asan=False):
    """Return executable path for a cmake target."""
    build_dir = ROOT / profile_dir(profile_key, asan)
    multi = is_multi_config(build_dir)
    exe_name = "gaudi_tests" if target == "gaudi_tests" else target

    if multi:
        subdir = "tests" if target == "gaudi_tests" else f"projects/{target}"
        return build_dir / subdir / config / f"{exe_name}.exe"
    if target == "gaudi_tests":
        return build_dir / "tests" / exe_name
    return build_dir / "projects" / target / exe_name


def resolve_target(name):
    """Return (profile_key, cmake_target)."""
    if name in TEST_ALIASES:
        return "test", "gaudi_tests"

    gl_enabled, vermeer_enabled = _enabled_targets_by_profile()

    if name in vermeer_enabled:
        return "vermeer", name

    if name in gl_enabled:
        return "gl", name

    all_projects = _all_project_targets()
    if name in all_projects:
        if name.startswith("vermeer_"):
            print(
                f"Warning: '{name}' exists in projects/ but is commented out in "
                f"CMakeLists.txt (GAUDI_WITH_VERMEER block). It may fail to build.",
                file=sys.stderr,
            )
            return "vermeer", name
        print(
            f"Warning: '{name}' exists in projects/ but is commented out in "
            f"CMakeLists.txt (GAUDI_WITH_GL block). It may fail to build.",
            file=sys.stderr,
        )
        return "gl", name

    sys.exit(f"Error: unknown target '{name}'. Run `python gaudi.py list`.")


# ── Actions ─────────────────────────────────────────────────────────────────


def configure(profile_key, config=None, asan=False):
    build_dir = ROOT / profile_dir(profile_key, asan)
    generator = _cmake_generator()
    cmd = [
        "cmake",
        "-S", str(ROOT),
        "-B", str(build_dir),
        "-G", generator,
        *profile_cmake_args(profile_key, asan),
    ]
    if config and not is_multi_config(build_dir) and not profile_build_type(profile_key, asan):
        cmd.append(f"-DCMAKE_BUILD_TYPE={config}")
    print(f">> {' '.join(cmd)}")
    return subprocess.run(cmd).returncode


def ensure_configured(profile_key, config, asan=False):
    prof_dir = profile_dir(profile_key, asan)
    build_dir = ROOT / prof_dir
    if not (build_dir / "CMakeCache.txt").exists():
        print(f"Build directory '{prof_dir}' not configured. Running cmake...")
        rc = configure(profile_key, config, asan)
        if rc != 0:
            sys.exit(f"cmake configure failed (exit {rc})")
        return

    if not is_multi_config(build_dir) and not profile_build_type(profile_key, asan):
        cached_type = _read_cache_value(build_dir, "CMAKE_BUILD_TYPE")
        if cached_type and cached_type != config:
            print(
                f"Build type mismatch ({cached_type} vs {config}) in "
                f"'{prof_dir}', reconfiguring..."
            )
            rc = configure(profile_key, config, asan)
            if rc != 0:
                sys.exit(f"cmake configure failed (exit {rc})")


def build(profile_key, target, config, asan=False):
    ensure_configured(profile_key, config, asan)
    build_dir = ROOT / profile_dir(profile_key, asan)
    cmd = ["cmake", "--build", str(build_dir), "--target", target]
    if is_multi_config(build_dir):
        cmd.extend(["--config", config])
    print(f">> {' '.join(cmd)}")
    return subprocess.run(cmd).returncode


def _vermeer_record_env(env, target, record_path=None):
    """Enable FFmpeg capture for Vermeer demos into dump/<target>.mp4."""
    stem = Path(record_path) if record_path else (ROOT / "dump" / target)
    if not stem.is_absolute():
        stem = ROOT / stem
    stem.parent.mkdir(parents=True, exist_ok=True)
    # Respect a pre-set VERMEER_RECORD_PATH unless --record-path was given.
    if record_path is not None or "VERMEER_RECORD_PATH" not in env:
        env["VERMEER_RECORD_PATH"] = str(stem)
    env["VERMEER_RECORD"] = "1"
    print(f">> VERMEER_RECORD=1 VERMEER_RECORD_PATH={env['VERMEER_RECORD_PATH']}")
    print(f">>   output: {env['VERMEER_RECORD_PATH']}.mp4 (sim frames only)")
    return env


def run_exe(
    profile_key,
    target,
    config,
    extra_args=None,
    asan=False,
    record=False,
    record_path=None,
):
    exe = exe_path(profile_key, target, config, asan)
    if not exe.exists():
        sys.exit(f"Error: executable not found at {exe}")
    argv = [str(exe)]
    if extra_args:
        argv.extend(extra_args)
    env = None
    if asan or record:
        env = os.environ.copy()
    if asan:
        env["ASAN_OPTIONS"] = "symbolize=1:halt_on_error=1"
        env["UBSAN_OPTIONS"] = "print_stacktrace=1:halt_on_error=1"
        if "ASAN_SYMBOLIZER_PATH" not in env:
            symbolizer = shutil.which("llvm-symbolizer")
            if symbolizer:
                env["ASAN_SYMBOLIZER_PATH"] = symbolizer
        print(
            f">> ASAN_OPTIONS={env['ASAN_OPTIONS']} "
            f"UBSAN_OPTIONS={env['UBSAN_OPTIONS']}"
        )
    if record:
        _vermeer_record_env(env, target, record_path)
    print(f">> {' '.join(argv)}")
    return subprocess.run(argv, env=env).returncode


# ── CLI ─────────────────────────────────────────────────────────────────────


def _warn_debug_with_asan(debug, asan):
    if debug and asan:
        print(
            "Warning: --asan pins CMAKE_BUILD_TYPE=asan; ignoring --debug.",
            file=sys.stderr,
        )


def cmd_test(args):
    _warn_debug_with_asan(args.debug, args.asan)
    config = effective_config("test", args.debug, args.asan)
    rc = build("test", "gaudi_tests", config, args.asan)
    if rc != 0:
        sys.exit(rc)
    sys.exit(run_exe("test", "gaudi_tests", config, asan=args.asan))


def cmd_run(args):
    profile, target = resolve_target(args.target)
    _warn_debug_with_asan(args.debug, args.asan)
    config = effective_config(profile, args.debug, args.asan)
    rc = build(profile, target, config, args.asan)
    if rc != 0:
        sys.exit(rc)
    forward = list(args.forward or [])
    if getattr(args, "shift_fraction", None) is not None:
        forward = ["--shift-fraction", str(args.shift_fraction)] + forward

    record = profile == "vermeer" and not getattr(args, "no_record", False)
    record_path = getattr(args, "record_path", None)
    if record_path is not None and profile != "vermeer":
        print(
            "Warning: --record-path ignored (recording is Vermeer-only).",
            file=sys.stderr,
        )
        record_path = None
    if record and not shutil.which("ffmpeg"):
        print(
            "Warning: ffmpeg not found on PATH; recording may fail to start.",
            file=sys.stderr,
        )

    sys.exit(
        run_exe(
            profile,
            target,
            config,
            forward if forward else None,
            args.asan,
            record=record,
            record_path=record_path,
        )
    )


def cmd_build(args):
    profile, target = resolve_target(args.target)
    _warn_debug_with_asan(args.debug, args.asan)
    config = effective_config(profile, args.debug, args.asan)
    sys.exit(build(profile, target, config, args.asan))


def cmd_configure(args):
    key = args.profile
    if key not in PROFILES:
        sys.exit(f"Unknown profile '{key}'. Choose from: {', '.join(PROFILES)}")
    _warn_debug_with_asan(getattr(args, "debug", False), args.asan)
    config = effective_config(key, getattr(args, "debug", False), args.asan)
    sys.exit(configure(key, config, args.asan))


def cmd_list(_args):
    gl_enabled, vermeer_enabled = _enabled_targets_by_profile()
    all_projects = _all_project_targets()
    gl_disabled = [t for t in all_projects if t not in gl_enabled and not t.startswith("vermeer_")]
    vermeer_disabled = [
        t for t in all_projects if t.startswith("vermeer_") and t not in vermeer_enabled
    ]

    print("Headless targets:")
    print(f"  tests          (alias: test, gaudi_tests)  [{profile_dir('test')}]")
    print(f"                 with --asan:                 [{profile_dir('test', True)}]")
    print()
    print(f"GL viewer targets (enabled):                  [{profile_dir('gl')}]")
    for t in gl_enabled:
        print(f"  {t}")
    if gl_disabled:
        print()
        print("GL viewer targets (commented out in CMakeLists.txt):")
        for t in gl_disabled:
            print(f"  {t}  (disabled)")
    print()
    print(f"Vermeer / WebGPU targets (enabled):           [{profile_dir('vermeer')}]")
    print(f"                 with --asan:                 [{profile_dir('vermeer', True)}]")
    for t in vermeer_enabled:
        print(f"  {t}")
    if vermeer_disabled:
        print()
        print("Vermeer targets (commented out in CMakeLists.txt):")
        for t in vermeer_disabled:
            print(f"  {t}  (disabled)")
    print()
    print(f"Build profiles: {', '.join(PROFILES)}")
    print("Use --asan on test/build/run/configure for ASan+UBSan (separate *_asan dirs).")
    print("Vermeer run defaults to FFmpeg dump/<target>.mp4; use --no-record to disable.")
    print(f"CMake generator: {_cmake_generator()}")


def main():
    parser = argparse.ArgumentParser(
        prog="gaudi",
        description="Build and run Gaudi targets.",
    )
    sub = parser.add_subparsers(dest="command")

    p_test = sub.add_parser("test", help="Build and run headless tests")
    p_test.add_argument("--debug", action="store_true", help="Debug config")
    p_test.add_argument(
        "--asan",
        action="store_true",
        help="ASan+UBSan build (CMAKE_BUILD_TYPE=asan, uses build_test_asan)",
    )

    p_run = sub.add_parser(
        "run",
        help="Build and run a target",
        epilog=(
            "Examples:\n"
            "  python gaudi.py run vermeer_dipole_tunneling_demo\n"
            "          (Vermeer: records to dump/<target>.mp4 by default)\n"
            "  python gaudi.py run vermeer_dipole_tunneling_demo --no-record\n"
            "  python gaudi.py run vermeer_dipole_tunneling_demo --record-path dump/my_run\n"
            "  python gaudi.py run spectral_modes_test --shift-fraction 0.35\n"
            "          (F=0 low spectrum, F=1 high, 0<F<1 interior via shift-invert)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_run.add_argument("target", help="Target name (e.g. vermeer_dipole_tunneling_demo)")
    p_run.add_argument("--debug", action="store_true", help="Debug config")
    p_run.add_argument(
        "--asan",
        action="store_true",
        help="ASan+UBSan build (CMAKE_BUILD_TYPE=asan, uses <profile>_asan dir)",
    )
    p_run.add_argument(
        "--no-record",
        action="store_true",
        help="Disable default FFmpeg recording for Vermeer targets",
    )
    p_run.add_argument(
        "--record-path",
        default=None,
        metavar="PATH",
        help=(
            "Vermeer recording stem (default: dump/<target> → dump/<target>.mp4)"
        ),
    )
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
    p_build.add_argument(
        "--asan",
        action="store_true",
        help="ASan+UBSan build (CMAKE_BUILD_TYPE=asan, uses <profile>_asan dir)",
    )

    p_cfg = sub.add_parser("configure", help="(Re)generate cmake for a profile")
    p_cfg.add_argument("profile", help=f"Profile: {', '.join(PROFILES)}")
    p_cfg.add_argument("--debug", action="store_true", help="Debug build type")
    p_cfg.add_argument(
        "--asan",
        action="store_true",
        help="ASan+UBSan build (CMAKE_BUILD_TYPE=asan, uses <profile>_asan dir)",
    )

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
