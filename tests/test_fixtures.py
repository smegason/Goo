"""pytest-blender fixtures tests."""

import os
import re
import shutil
import subprocess

import pytest


try:
    from pytest_blender.test import pytest_blender_unactive
except ImportError:
    pytest_blender_unactive = False


@pytest.mark.skipif(
    pytest_blender_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_blender_executable(blender_executable):
    assert blender_executable
    blender_executable_path = shutil.which(blender_executable)
    if not os.path.isfile(blender_executable_path):
        assert os.path.islink(blender_executable_path)
    else:
        assert os.path.isfile(blender_executable_path)
    assert "blender" in os.path.basename(blender_executable_path).lower()


@pytest.mark.skipif(
    pytest_blender_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_blender_python_executable(blender_python_executable):
    assert blender_python_executable
    if not os.path.isfile(blender_python_executable):
        assert os.path.islink(blender_python_executable)
    else:
        assert os.path.isfile(blender_python_executable)
    assert "python" in os.path.basename(blender_python_executable)


@pytest.mark.skipif(
    pytest_blender_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_blender_version(blender_executable, blender_version):
    stdout = subprocess.check_output([blender_executable, "--version"])
    expected_blender_version = stdout.decode("utf-8").splitlines()[0].split(" ")[1]

    assert expected_blender_version == blender_version
    assert re.match(r"\d+\.\d", blender_version)


@pytest.mark.skipif(
    pytest_blender_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_blender_python_version(blender_python_version, blender_python_executable):
    blender_python_version_stdout = subprocess.check_output(
        [
            blender_python_executable,
            "--version",
        ]
    )
    expected_blender_python_version = (
        blender_python_version_stdout.decode("utf-8").splitlines()[0].split(" ")[1]
    )

    assert blender_python_version == expected_blender_python_version
    assert re.match(r"\d+\.\d\.?\d*", blender_python_version)