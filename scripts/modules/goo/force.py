import bpy
from goo.utils import *


class Force(BlenderObject):
    def __init__(self, obj, type="FORCE"):
        super(Force, self).__init__(obj)
        self.type = type

        # Instantiate force field object
        if not self._obj.field:
            with bpy.context.temp_override(active_object=self._obj, object=self._obj):
                bpy.ops.object.forcefield_toggle()
        self.enable()

    def enable(self):
        self._obj.field.type = self.type

    def disable(self):
        if self._obj.field:
            self._obj.field.type = "NONE"

    def enabled(self):
        return self._obj.field and self._obj.field.type == self.type

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
        super(AdhesionForce, self).__init__(obj, "FORCE")

    @property
    def strength(self):
        return -self._obj.field.strength

    @strength.setter
    def strength(self, strength):
        self._obj.field.strength = -strength


def create_force(
    name, loc, strength, type="FORCE", falloff=0, min_dist=None, max_dist=None
) -> Force:
    obj = bpy.data.objects.new(name, None)
    obj.location = loc
    force = Force(obj)

    force.strength = strength
    force.falloff = falloff
    force.min_dist = min_dist
    force.max_dist = max_dist

    ForceCollection.global_forces().add_force(force)
    return force


def create_turbulence(name, loc, strength) -> Force:
    force = create_force(name, loc, strength, max_dist=10, type="TURBULENCE")

    force.obj.field.size = 7
    force.obj.field.noise = 0
    force.obj.field.seed = 1
    return force


def get_adhesion(strength, obj=None, name=None, loc=(0, 0, 0)) -> AdhesionForce:
    if obj is None:
        obj = bpy.data.objects.new(name, None)
        obj.location = loc
    adhesion_force = AdhesionForce(obj)

    adhesion_force.strength = strength
    adhesion_force.min_dist = 0.6
    adhesion_force.max_dist = 1.4
    adhesion_force.falloff = 0.5
    return adhesion_force


def get_motion(name, loc, strength) -> Force:
    obj = bpy.data.objects.new(name, None)
    obj.location = loc
    force = Force(obj)

    force.strength = strength
    force.falloff = 2
    return force


class ForceCollection:
    _global_forces = None

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

    @staticmethod
    def global_forces():
        if ForceCollection._global_forces is None:
            ForceCollection._global_forces = ForceCollection("globals")
            bpy.context.scene.collection.children.link(
                ForceCollection._global_forces.collection
            )
        return ForceCollection._global_forces
