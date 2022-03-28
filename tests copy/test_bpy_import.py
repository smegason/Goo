import pytest


try:
    from goo.test import goo_unactive
except ImportError:
    goo_unactive = False


@pytest.mark.skipif(
    goo_unactive,
    reason="Requires testing loading the pytest-blender plugin.",
)
def test_bpy_import():
    import bpy  # noqa F401