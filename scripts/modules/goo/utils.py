from functools import reduce

import bpy, bmesh
from mathutils import *


class BlenderObject:
    def __init__(self, obj):
        self._obj = obj

    @property
    def obj(self):
        return self._obj

    @property
    def name(self):
        return self._obj.name

    @name.setter
    def name(self, name):
        self._obj.name = name

    @property
    def loc(self):
        return self._obj.location

    @loc.setter
    def loc(self, loc):
        self._obj.location = loc


def create_mesh(
    name,
    loc,
    mesh="icosphere",
    size=1,
    rotation=(0, 0, 0),
    scale=(1, 1, 1),
    subdivisions=2,
    **kwargs,
):
    bm = bmesh.new()

    if isinstance(rotation, tuple):
        rotation = Euler(rotation)
    elif isinstance(rotation, Quaternion):
        rotation = rotation.to_euler()

    if mesh == "icosphere":
        bmesh.ops.create_icosphere(bm, subdivisions=subdivisions, radius=size, **kwargs)
    elif mesh == "plane":
        bmesh.ops.create_grid(bm, x_segments=1, y_segments=1, size=size / 2, **kwargs)
    elif mesh == "monkey":
        bmesh.ops.create_monkey(bm, **kwargs)
    else:
        raise ValueError("""mesh must be one of "icosphere", "plane", or "monkey".""")

    me = bpy.data.meshes.new(f"{name}_mesh")
    bm.to_mesh(me)
    bm.free()

    obj = bpy.data.objects.new(name, me)
    obj.location = loc
    obj.rotation_euler = rotation
    obj.scale = scale

    return obj


# Utility functions for settings
def rsetattr(obj, attr, val):
    pre, _, post = attr.rpartition(".")
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)


def rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)

    return reduce(_getattr, [obj] + attr.split("."))
