# test_goo.py
# this file is read by pytest to test all functions in goo

import bpy

# importing goo
import importlib.util
spec = importlib.util.spec_from_file_location("goo", "scripts/modules/goo/goo.py")
goo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(goo)

# show Blender version
print("***------- Testing Goo -------***")
version = bpy.app.version_string
print("Blender version=" + version)

print("Enable Extra objects Add-on")
bpy.ops.preferences.addon_enable(module='add_mesh_extra_objects')


def test_sample_func():
    print("test sample function------------")

    ret = goo.sample_func()
    print("Goo ret=" + ret)
    assert True


def test_calculate_volume():
    # make cell

    # volume = goo.calculate_volume(obj)
    # assert volume = ??
    print("test calculate volume")

    cell = goo.Cell(name_string="Cell1_", loc=(0, 0, 0))
    goo.make_cell(cell)

    volume = goo.calculate_volume(cell.get_blender_object())

    print("Volume of cell=")
    print(volume)

    assert(volume > 4 and volume < 4.1)
