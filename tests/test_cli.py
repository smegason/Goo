import os
import subprocess
import sys

import pytest


try:
    from pytest_blender.test import pytest_blender_active
except ImportError:
    pytest_blender_active = True


@pytest.mark.skipif(
    pytest_blender_active,
    reason="Requires testing without loading the pytest-blender plugin.",
)
@pytest.mark.skipif(
    not os.environ.get("BLENDER_EXECUTABLE"),
    reason="Environment variable 'BLENDER_EXECUTABLE' must be set.",
)
def test_pytest_blender_cli():
    proc = subprocess.run(
        [
            sys.executable,
            os.path.join("pytest_blender", "__main__.py"),
            "--blender-executable",
            os.environ["BLENDER_EXECUTABLE"],
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert proc.stderr == b""

    stdout = proc.stdout.decode("utf-8")
    assert stdout.count("\n") == 1
    assert stdout.endswith("\n")
    assert stdout.split(os.sep)[-1].startswith("python")
    assert os.path.isfile(stdout.rstrip("\n"))