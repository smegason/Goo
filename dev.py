import subprocess
import argparse
import shutil
import os

ADDON_NAME = 'sample_addon'


def do_build(install_at=None, include_tests=False):
    shutil.rmtree("out", True)
    target_dir = os.path.join("out", ADDON_NAME)
    ignore_files = [".gitignore", "dev.py", "README.md", "CONTRIBUTING.md", "setup.cfg"]

    shutil.copytree(
        "npanel",
        f"{target_dir}/npanel",
        ignore=shutil.ignore_patterns("__pycache__"),
    )

    for item in os.listdir():
        if os.path.isdir(item):
            continue  # we copied directories above
        if item in ignore_files:
            continue
        if include_tests is False and item == "tests.py":
            continue
        if include_tests is False and item.startswith("test_"):
            continue  # we do not include test files
        shutil.copy(item, os.path.join(target_dir, item))

    # CREATE ZIP
    shutil.make_archive(target_dir, "zip", "out", ADDON_NAME)

    if install_at is not None:
        install_location = os.path.join(install_at, ADDON_NAME)
        shutil.rmtree(install_location, ignore_errors=True)
        shutil.copytree(target_dir, install_location)


def run_tests():
    test = subprocess.Popen(
        [
            "blender",
            "--background",
            "-noaudio",
            "--factory-startup",
            "--python-exit-code",
            "1",
            "--python",
            "tests/tests.py",
        ]
    )
    test.wait()

    if test.returncode == 1:
        exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        default="build",
        choices=["build", "test"],
        help="""
    TEST = build with test files and run tests
    BUILD = copy relevant files into ./out/blenderkit.
    """,
    )
    parser.add_argument(
        "--install-at",
        type=str,
        default=None,
        help="If path is specified, then builded addon will \
            be copied to that location.",
    )
    args = parser.parse_args()

    if args.command == "build":
        do_build(args.install_at)
    elif args.command == "test":
        do_build(args.install_at, include_tests=True)
        run_tests()
    elif args.command == "bundle":
        pass  # bundle_dependencies()
    elif args.command == "format":
        pass  # format_code()
    else:
        parser.print_help()
