import bpy
from goo.utils import *


class Force(BlenderObject):
    def __init__(self, obj):
        super(Force, self).__init__(obj)
        self.type = "force"

    def enable(self):
        if not self._obj.field:
            with bpy.context.temp_override(active_object=self._obj):
                bpy.ops.object.forcefield_toggle()
        else:
            self._obj.field.type = "FORCE"

    def disable(self):
        if self._obj.field:
            self._obj.field.type = "NONE"

    def enabled(self):
        return self._obj.field and self._obj.field.type == "FORCE"

    @property
    def strength(self):
        return self._obj.field.strength

    @strength.setter
    def strength(self, strength):
        self._obj.field.strength = strength

    @property
    def falloff(self):
        return self._obj.field.falloff

    @falloff.setter
    def falloff(self, falloff):
        self._obj.field.falloff_power = falloff

    @property
    def min_dist(self):
        if self._obj.field.use_min_distance:
            return self._obj.field.min_dist
        return None

    @min_dist.setter
    def min_dist(self, min_dist):
        if min_dist is None:
            self._obj.field.use_min_distance = False
        else:
            self._obj.field.use_min_distance = True
            self._obj.field.distance_min = min_dist

    @property
    def max_dist(self):
        if self._obj.field.use_max_distance:
            return self._obj.field.max_dist
        return None

    @max_dist.setter
    def max_dist(self, max_dist):
        if max_dist is None:
            self._obj.field.use_max_distance = False
        else:
            self._obj.field.use_max_distance = True
            self._obj.field.distance_max = max_dist


class AdhesionForce(Force):
    def __init__(self, obj):
        super(AdhesionForce, self).__init__(obj)
        self.type = "adhesion"

    @property
    def strength(self):
        return -self._obj.field.strength

    @strength.setter
    def strength(self, strength):
        self._obj.field.strength = -strength


def create_force(name, loc, strength, falloff=0, min_dist=None, max_dist=None) -> Force:
    obj = bpy.data.objects.new(name, None)
    obj.location = loc

    force = Force(obj)
    force.enable()
    force.strength = strength
    force.falloff = falloff
    force.min_dist = min_dist
    force.max_dist = max_dist

    bpy.context.collection.objects.link(obj)
    return force


def create_adhesion(strength, obj=None, name="", loc=(0, 0, 0)) -> AdhesionForce:
    if obj is None:
        assert name != ""
        obj = bpy.data.objects.new(name, None)
        obj.location = loc

    adhesion_force = AdhesionForce(obj)
    adhesion_force.enable()

    adhesion_force.strength = strength
    adhesion_force.min_dist = 0.6
    adhesion_force.max_dist = 1.4
    adhesion_force.falloff = 0.5
    return adhesion_force


def create_motion(name, loc, strength) -> Force:
    obj = bpy.data.objects.new(name, None)
    obj.location = loc

    force = Force(obj)
    force.enable()

    force.strength = strength
    force.falloff = 2
    return force


class ForceCollection:
    def __init__(self, name: str):
        self._col = bpy.data.collections.new(name)
        self._forces = []

    @property
    def name(self) -> str:
        return self._col.name

    @property
    def collection(self) -> bpy.types.Collection:
        return self._col

    def add_force(self, force: Force):
        self._col.objects.link(force.obj)
        self._forces.append(force)

    def remove_force(self, force: Force):
        self._col.objects.unlink(force.obj)
        self._forces.remove(force)

    @property
    def forces(self):
        return self._forces
