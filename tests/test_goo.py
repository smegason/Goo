# test_goo.py
# this file is read by pytest to test all functions in goo

import bpy
import math

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


def test_calculate_volume():
    print("Test calculate_volume")

    # make cell
    cell = goo.Cell(name_string="Cell1_", loc=(0, 0, 0))
    goo.make_cell(cell)

    # get volume
    volume = goo.calculate_volume(cell.get_blender_object())
    print("Volume of cell=")
    print(volume)
    delta = math.abs(4/3 * math.pi * 1**3 - volume)
    print("Delta=")
    print(delta)

    # volume should be 4/3 pi 1^3 = 4.19 but runs low (?)
    assert(delta < 0.2)


def test_get_major_axis():
    print("Test get_major_axis")

    # make cell
    cell = goo.Cell(name_string="Cell1_", loc=(0, 0, 0))
    cell.data['size'] = (2, 1, 1)
    goo.make_cell(cell)

    # get major axis
    axis = goo.get_major_axis(cell.get_blender_object())
    print("Axis of cell=")
    print(axis[0])
    print(axis[1])
    print(axis[2])

    # calculate difference
    delta = math.abs(-1 - axis[0])
    print("Delta=")
    print(delta)

    # axis should align with x-axis
    assert(delta < 0.1)
