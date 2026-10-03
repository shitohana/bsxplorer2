"""Validate distributions from isolated installs outside the source checkout."""

import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import tomllib

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "python"


def run(*command: str, cwd: Path, env: dict[str, str]) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, env=env, check=True)


def validate(
    wheel: Path, workspace: Path, requirements: Path, env: dict[str, str]
) -> None:
    workspace.mkdir()
    venv = workspace / ".venv"
    run("uv", "venv", "--python", sys.executable, str(venv), cwd=workspace, env=env)
    interpreter = str(venv / "bin" / "python")
    run(
        "uv",
        "pip",
        "install",
        "--python",
        interpreter,
        "--require-hashes",
        "-r",
        str(requirements),
        cwd=workspace,
        env=env,
    )
    run(
        "uv",
        "pip",
        "install",
        "--python",
        interpreter,
        "--no-deps",
        str(wheel),
        cwd=workspace,
        env=env,
    )
    shutil.copytree(
        PROJECT / "tests",
        workspace / "tests",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    # -I and the temporary working directory prevent checkout imports from
    # concealing missing files or a stale extension in a distribution.
    probe = """
from pathlib import Path
import bsx2
from bsx2 import _bsx2
root = Path(bsx2.__file__).parent
assert '.venv' in root.parts, root
assert (root / 'py.typed').is_file()
for name in ['__init__.pyi', '_bsx2.pyi', 'io.pyi', 'types.pyi']:
    assert (root / name).is_file(), name
assert _bsx2.BsxBatch.empty().is_empty()
print('Installed package:', root)
"""
    run(interpreter, "-I", "-c", probe, cwd=workspace, env=env)
    run(
        interpreter,
        "-I",
        "-m",
        "pytest",
        "tests",
        "-q",
        "--import-mode=importlib",
        cwd=workspace,
        env=env,
    )
    run(
        interpreter,
        "-I",
        "-m",
        "mypy.stubtest",
        "bsx2._bsx2",
        "--concise",
        cwd=workspace,
        env=env,
    )
    run(
        interpreter,
        "-I",
        "-m",
        "mypy.stubtest",
        "bsx2.io",
        "bsx2.types",
        "--concise",
        cwd=workspace,
        env=env,
    )
    run(
        str(venv / "bin" / "ty"),
        "check",
        "tests",
        "--python",
        interpreter,
        cwd=workspace,
        env=env,
    )


def main() -> None:
    tag = f"cp{sys.version_info.major}{sys.version_info.minor}"
    wheels = sorted((ROOT / "dist").glob(f"bsx2-*-{tag}-{tag}-*.whl"))
    sdists = sorted((ROOT / "dist").glob("bsx2-*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise RuntimeError(
            "Expected one matching wheel and one sdist in dist/; remove stale bsx2 artifacts"
        )
    env = os.environ.copy()
    for name in ["PYTHONPATH", "MYPYPATH", "VIRTUAL_ENV", "CONDA_PREFIX"]:
        env.pop(name, None)
    with tempfile.TemporaryDirectory(prefix="bsx2-package-") as temporary:
        workspace = Path(temporary)
        requirements = workspace / "requirements.txt"
        run(
            "uv",
            "export",
            "--quiet",
            "--project",
            str(PROJECT),
            "--locked",
            "--no-emit-project",
            "--output-file",
            str(requirements),
            cwd=workspace,
            env=env,
        )
        validate(wheels[0], workspace / "wheel", requirements, env)
        source = workspace / "source"
        source.mkdir()
        with tarfile.open(sdists[0]) as archive:
            archive.extractall(source, filter="data")
        extracted = next(source.iterdir())
        if not (extracted / "Cargo.lock").is_file():
            raise RuntimeError(
                "Source distribution is missing the workspace Cargo.lock"
            )
        rebuilt = workspace / "rebuilt"
        # Reuse dependency compilation, but build the extracted core and wrapper.
        # No source paths in the original checkout are used by this manifest.
        source_env = env | {
            "CARGO_TARGET_DIR": str(ROOT / "target"),
            "CARGO_NET_OFFLINE": "true",
        }
        before = tomllib.loads((extracted / "Cargo.lock").read_text())
        # Recompile our crates from the archive while retaining dependency
        # compilation. A reused checkout binary alone cannot verify the sdist.
        run(
            "cargo",
            "clean",
            "--manifest-path",
            str(extracted / "python" / "Cargo.toml"),
            "--profile",
            "pyrelease",
            "-p",
            "bsxplorer2",
            "-p",
            "bsxplorer2-py",
            cwd=extracted,
            env=source_env,
        )
        run(
            "uv",
            "build",
            "--wheel",
            "--no-build-isolation",
            "--python",
            sys.executable,
            "--out-dir",
            str(rebuilt),
            str(extracted),
            cwd=extracted,
            env=source_env,
        )
        # Maturin removes the CLI workspace member from sdists. Cargo can prune
        # its unused lock entries, but must preserve every registry version.
        after = tomllib.loads((extracted / "Cargo.lock").read_text())

        def registry_packages(lock):
            return {
                (p["name"], p["version"], p.get("source"), p.get("checksum"))
                for p in lock["package"]
                if "source" in p
            }

        if not registry_packages(after) <= registry_packages(before):
            raise RuntimeError("Source build changed locked dependency versions")
        validate(next(rebuilt.glob("*.whl")), workspace / "sdist", requirements, env)
    print("Wheel and source distribution validation passed.", flush=True)


if __name__ == "__main__":
    main()
