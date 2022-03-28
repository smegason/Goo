"""pytest-blender tests configuration."""

import io
import logging
import os
import sys
from contextlib import redirect_stdout

import pytest


goo_logger = logging.getLogger("goo")
goo_logger.setLevel(logging.DEBUG)
goo_logger.addHandler(logging.StreamHandler())

try:
    from Goo.test import goo_active
except ImportError:
    # executing pytest from Python Blender executable, the plugin is active
    goo_active = True

if goo_active:

    @pytest.fixture(scope="session", autouse=True)
    def _register_addons(request, install_addons_from_dir, disable_addons):
        addons_dir = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            "addons",
        )

        f = io.StringIO()
        with redirect_stdout(f):
            addon_module_names = install_addons_from_dir(addons_dir)
        yield
        sys.stdout.write("\n")
        with redirect_stdout(f):
            disable_addons(addon_module_names)