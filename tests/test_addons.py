"""Example addons testing."""

import pytest


try:
    from goo.test import goo_unactive
except ImportError:
    goo_unactive = False


@pytest.mark.skipif(
    goo_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_basic_addon():
    import _bpy
    import addon_utils

    installed_addons = [addon.__name__ for addon in addon_utils.modules()]
    assert "goo_basic" in installed_addons
    assert "__init__" not in installed_addons

    operator_classes = [cls.__name__ for cls in _bpy.types.Operator.__subclasses__()]
    assert "PytestBlenderObjectMoveX" in operator_classes