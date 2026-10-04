"""Check the wheel and source distribution inside the devenv environment."""

import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(*args, cwd, env):
    print("+", *args, flush=True)
    subprocess.run(args, cwd=cwd, env=env, check=True)


def registry_packages(source):
    lock = tomllib.loads((source / "Cargo.lock").read_text())
    return {
        (p["name"], p["version"], p["source"], p.get("checksum"))
        for p in lock["package"]
        if "source" in p
    }


def validate(wheel, work, requirements, env):
    venv = work / ".venv"
    run("uv", "venv", "--python", sys.executable, str(venv), cwd=work, env=env)
    python = str(venv / "bin/python")
    run(
        "uv",
        "pip",
        "install",
        "--python",
        python,
        "--require-hashes",
        "-r",
        str(requirements),
        cwd=work,
        env=env,
    )
    run(
        "uv",
        "pip",
        "install",
        "--python",
        python,
        "--no-deps",
        str(wheel),
        cwd=work,
        env=env,
    )
    run(
        python,
        "-I",
        "-c",
        """
from pathlib import Path
import sys
import bsx2
from bsx2 import _bsx2
root = Path(bsx2.__file__).parent
assert root.is_relative_to(Path(sys.prefix)), root
assert (root / 'py.typed').is_file()
for name in ('__init__.pyi', '_bsx2.pyi', 'io.pyi', 'types.pyi'):
    assert (root / name).is_file(), name
assert _bsx2.BsxBatch.empty().is_empty()
print('Installed package:', root)
""",
        cwd=work,
        env=env,
    )
    shutil.copytree(
        ROOT / "python/tests",
        work / "tests",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    run(python, "-I", "-m", "pytest", "tests", "-q", cwd=work, env=env)
    for modules in [("bsx2._bsx2",), ("bsx2.io", "bsx2.types")]:
        run(
            python,
            "-I",
            "-m",
            "mypy.stubtest",
            *modules,
            "--concise",
            cwd=work,
            env=env,
        )


def main():
    (wheel,) = (ROOT / "dist").glob("bsx2-*.whl")
    (sdist,) = (ROOT / "dist").glob("bsx2-*.tar.gz")
    env = dict(os.environ)
    for name in ("PYTHONPATH", "MYPYPATH", "VIRTUAL_ENV", "UV_PROJECT_ENVIRONMENT"):
        env.pop(name, None)
    with tempfile.TemporaryDirectory(prefix="bsx2-package-") as directory:
        work = Path(directory)
        requirements = work / "requirements.txt"
        run(
            "uv",
            "export",
            "--quiet",
            "--project",
            str(ROOT / "python"),
            "--locked",
            "--no-emit-project",
            "--output-file",
            str(requirements),
            cwd=work,
            env=env,
        )
        validate(wheel, work, requirements, env)
        with tarfile.open(sdist) as archive:
            archive.extractall(work / "source", filter="data")
        (source,) = (work / "source").iterdir()
        before = registry_packages(source)
        # This cache contains only archive-source builds, so a checkout build
        # cannot conceal missing sources. Cargo owns its incremental reuse.
        run(
            "uv",
            "build",
            "--wheel",
            "--no-build-isolation",
            "--python",
            sys.executable,
            "--out-dir",
            str(work / "rebuilt"),
            str(source),
            cwd=source,
            env=env
            | {
                "CARGO_TARGET_DIR": str(ROOT / "target/package-check"),
                "CARGO_NET_OFFLINE": "true",
            },
        )
        assert registry_packages(source) <= before, (
            "Source build changed locked dependencies"
        )
        (rebuilt,) = (work / "rebuilt").glob("*.whl")
        shutil.rmtree(work / ".venv")
        shutil.rmtree(work / "tests")
        validate(rebuilt, work, requirements, env)
    print("Wheel and source distribution validation passed.")


if __name__ == "__main__":
    main()
