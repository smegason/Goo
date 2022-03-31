# test_goo.py
# this file is read by pytest to test all functions in goo

import bpy

# importing goo
import importlib.util
spec = importlib.util.spec_from_file_location("goo", "scripts/modules/goo/goo.py")
goo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(goo)


def test_sample_func():
    print("test sample function------------")

    version = bpy.app.version_string
    print("Blender version=" + version)
    ret = goo.sample_func()
    print("Goo ret=" + ret)
    assert True


def test_calculate_volume():
    # make cell

    # volume = goo.calculate_volume(obj)
    # assert volume = ??
    print("test calculate volume")

    assert True
