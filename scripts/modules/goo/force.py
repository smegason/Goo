import bpy
from goo.utils import *


class Force(BlenderObject):
    def __init__(self, obj, type="FORCE"):
        super(Force, self).__init__(obj)
        self.type = type

        # Instantiate force field object
        if not self.obj.field:
            with bpy.context.temp_override(active_object=self.obj, object=self.obj):
                bpy.ops.object.forcefield_toggle()
        self.enable()

    def enable(self):
        self.obj.field.type = self.type

    def disable(self):
        if self.obj.field:
            self.obj.field.type = "NONE"

    def enabled(self):
        return self.obj.field and self.obj.field.type == self.type

    @property
    def strength(self):
        return self.obj.field.strength

    @strength.setter
    def strength(self, strength):
        self.obj.field.strength = strength

    @property
    def falloff(self):
        return self.obj.field.falloff

    @falloff.setter
    def falloff(self, falloff):
        self.obj.field.falloff_power = falloff

    @property
    def min_dist(self):
        if self.obj.field.use_min_distance:
            return self.obj.field.min_dist
        return None

    @min_dist.setter
    def min_dist(self, min_dist):
        if min_dist is None:
            self.obj.field.use_min_distance = False
        else:
            self.obj.field.use_min_distance = True
            self.obj.field.distance_min = min_dist

    @property
    def max_dist(self):
        if self.obj.field.use_max_distance:
            return self.obj.field.max_dist
        return None

    @max_dist.setter
    def max_dist(self, max_dist):
        if max_dist is None:
            self.obj.field.use_max_distance = False
        else:
            self.obj.field.use_max_distance = True
            self.obj.field.distance_max = max_dist


class AdhesionForce(Force):
    def __init__(self, obj):
        super(AdhesionForce, self).__init__(obj, "FORCE")

    @property
    def strength(self):
        return -self.obj.field.strength

    @strength.setter
    def strength(self, strength):
        self.obj.field.strength = -strength


class MotionForce(Force):
    def __init__(self, obj: bpy.types.Object):
        super(MotionForce, self).__init__(obj, "FORCE")
        obj.field.shape = "PLANE"
        obj.field.apply_to_rotation = False

    @property
    def strength(self):
        return -self.obj.field.strength

    @strength.setter
    def strength(self, strength):
        self.obj.field.strength = -strength

    def set_loc(self, new_loc, target_loc):
        dir = target_loc - new_loc
        self.loc = new_loc
        self.obj.rotation_euler = dir.to_track_quat("Z", "X").to_euler()


def create_force(
    name, loc, strength, type="FORCE", falloff=0, min_dist=None, max_dist=None
) -> Force:
    obj = bpy.data.objects.new(name, None)
    obj.location = loc
    force = Force(obj, type)

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


def get_motion(name, loc, strength) -> MotionForce:
    obj = bpy.data.objects.new(name, None)
    obj.location = loc
    force = MotionForce(obj)

    force.strength = strength
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


class Boundary(BlenderObject):
    def setup_physics(self):
        BoundaryCollisionConstructor().construct(self.obj)


def create_boundary(loc, size, mesh="icosphere"):
    obj = create_mesh("Boundary", loc, mesh=mesh, size=size)
    bpy.context.scene.collection.objects.link(obj)

    boundary = Boundary(obj)
    boundary.setup_physics()
    boundary.hide()
    return boundary
