"""Example addons testing."""

import pytest


try:
    from pytest_blender.test import pytest_blender_unactive
except ImportError:
    pytest_blender_unactive = False


@pytest.mark.skipif(
    pytest_blender_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_basic_addon():
    import _bpy
    import addon_utils

    installed_addons = [addon.__name__ for addon in addon_utils.modules()]
    assert "pytest_blender_basic" in installed_addons
    assert "__init__" not in installed_addons

    operator_classes = [cls.__name__ for cls in _bpy.types.Operator.__subclasses__()]
    assert "PytestBlenderObjectMoveX" in operator_classes