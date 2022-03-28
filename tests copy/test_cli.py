import os
import subprocess
import sys

import pytest


try:
    from goo.test import goo_active
except ImportError:
    goo_active = True


@pytest.mark.skipif(
    goo_active,
    reason="Requires testing without loading the pytest-blender plugin.",
)
@pytest.mark.skipif(
    not os.environ.get("BLENDER_EXECUTABLE"),
    reason="Environment variable 'BLENDER_EXECUTABLE' must be set.",
)
def test_goo_cli():
    proc = subprocess.run(
        [
            sys.executable,
            os.path.join("goo", "__main__.py"),
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