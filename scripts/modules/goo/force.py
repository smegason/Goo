from typing import Optional
from typing_extensions import override

import bpy
from goo.utils import *


class Force(BlenderObject):
    """A force.

    Forces are represented in Blender by force field objects. They interact
    with cells to influence their motion.

    Args:
        obj: Blender object to be used as a representation for the force.

    Attributes:
        type (str): Type of force.
    """

    def __init__(self, obj: bpy.types.Object, type="FORCE"):
        super(Force, self).__init__(obj)
        self.type = type

        # Instantiate force field object
        if not self.obj.field:
            with bpy.context.temp_override(active_object=self.obj, object=self.obj):
                bpy.ops.object.forcefield_toggle()
        self.enable()

    def enable(self):
        """Enables the force."""
        self.obj.field.type = self.type

    def disable(self):
        """Disables the force."""
        if self.obj.field:
            self.obj.field.type = "NONE"

    def enabled(self) -> bool:
        """Checks if the force field is enabled."""
        return self.obj.field and self.obj.field.type == self.type

    @property
    def strength(self) -> int:
        """Strength of the force field."""
        return self.obj.field.strength

    @strength.setter
    def strength(self, strength: int):
        self.obj.field.strength = strength

    @property
    def falloff(self) -> float:
        """Falloff power of the force. Strength at distance :math:`r` is given by
        :math:`\\text{strength} / r ^ \\text{falloff}`."""
        return self.obj.field.falloff

    @falloff.setter
    def falloff(self, falloff):
        self.obj.field.falloff_power = falloff

    @property
    def min_dist(self) -> float:
        """Minimum distance an object must be from a force to be affected."""
        if self.obj.field.use_min_distance:
            return self.obj.field.min_dist
        return None

    @min_dist.setter
    def min_dist(self, min_dist: Optional[float]):
        if min_dist is None:
            self.obj.field.use_min_distance = False
        else:
            self.obj.field.use_min_distance = True
            self.obj.field.distance_min = min_dist

    @property
    def max_dist(self) -> float:
        """Maximum distance an object can be from a force to be affected."""
        if self.obj.field.use_max_distance:
            return self.obj.field.max_dist
        return None

    @max_dist.setter
    def max_dist(self, max_dist: Optional[float]):
        if max_dist is None:
            self.obj.field.use_max_distance = False
        else:
            self.obj.field.use_max_distance = True
            self.obj.field.distance_max = max_dist


class AdhesionForce(Force):
    """An adhesion force."""

    def __init__(self, obj):
        super(AdhesionForce, self).__init__(obj, "FORCE")

    @override
    @property
    def strength(self) -> int:
        return -self.obj.field.strength

    @strength.setter
    def strength(self, strength: int):
        self.obj.field.strength = -strength


class MotionForce(Force):
    """A motion force."""

    def __init__(self, obj: bpy.types.Object):
        super(MotionForce, self).__init__(obj, "FORCE")
        obj.field.shape = "PLANE"
        obj.field.apply_to_rotation = False

    @override
    @property
    def strength(self):
        return -self.obj.field.strength

    @strength.setter
    def strength(self, strength):
        self.obj.field.strength = -strength

    def set_loc(self, new_loc: Vector, target_loc: Vector):
        """Set location of a motion force, towards which a cell will move.

        Args:
            new_loc: Location of the motion force.
            target_loc: Location of the cell upon which the motion acts.
        """
        dir = target_loc - new_loc
        self.loc = new_loc
        self.obj.rotation_euler = dir.to_track_quat("Z", "X").to_euler()


def create_force(
    name: str,
    loc: tuple,
    strength: int,
    type: str = "FORCE",
    falloff: float = 0,
    min_dist: float = None,
    max_dist: float = None,
) -> Force:
    """Creates a new force field.

    Args:
        name: The name of the force field object.
        loc: The location of the force field object.
        strength: The strength of the force field.
        type: The type of the force field.
        falloff: The falloff power of the force field.
        min_dist: The minimum distance for the force field.
        max_dist: The maximum distance for the force field.

    Returns:
        Force: The created force field object.
    """
    obj = bpy.data.objects.new(name, None)
    obj.location = loc
    force = Force(obj, type)

    force.strength = strength
    force.falloff = falloff
    force.min_dist = min_dist
    force.max_dist = max_dist

    ForceCollection.global_forces().add_force(force)
    return force


def create_turbulence(name: str, loc: tuple, strength: int) -> Force:
    """Creates a new turbulence force. See :func:`create_force`."""
    force = create_force(name, loc, strength, max_dist=10, type="TURBULENCE")

    force.obj.field.size = 7
    force.obj.field.noise = 0
    force.obj.field.seed = 1
    return force


def get_adhesion(
    strength: int,
    obj: Optional[bpy.types.Object] = None,
    name: str = None,
    loc: tuple = (0, 0, 0),
) -> AdhesionForce:
    """Creates a new adhesion force.

    Adhesion forces can either be created from cells, in which they are
    homotypic forces, Or they can be created de novo, in which they are
    heterotypic adhesion forces meant to allow two different cell types to
    interact with each other.

    Args:
        strength: Strength of the adhesion force.
        obj: Cell to use as origin of the adhesion force. 
            If None, a new object is created.
        name: Name of the adhesion force.
        loc: Initial location of the adhesion force.
    """
    if obj is None:
        obj = bpy.data.objects.new(name, None)
        obj.location = loc
    adhesion_force = AdhesionForce(obj)

    adhesion_force.strength = strength
    adhesion_force.min_dist = 0.6
    adhesion_force.max_dist = 1.4
    adhesion_force.falloff = 0.5
    return adhesion_force


def get_motion(name: str, loc: tuple, strength: int) -> MotionForce:
    """Creates a new motion force.

    Args:
        name: Name of the motion force.
        loc: Initial location of the motion force.
        strength: Strength of the motion force.
    """
    obj = bpy.data.objects.new(name, None)
    obj.location = loc
    force = MotionForce(obj)

    force.strength = strength
    return force


class ForceCollection:
    """A class representing a collection of forces.

    A force collection is represented by Blender collections of Blender forces.
    """

    _global_forces = None

    def __init__(self, name: str):
        self._col = bpy.data.collections.new(name)
        self._forces = []

    @property
    def name(self) -> str:
        """Name of the collection of forces."""
        return self._col.name

    @property
    def collection(self) -> bpy.types.Collection:
        """The underlying Blender collection."""
        return self._col

    def add_force(self, force: Force):
        """Add a Force to the collection."""
        self._col.objects.link(force.obj)
        self._forces.append(force)

    def remove_force(self, force: Force):
        """Remove a Force from the collection."""
        self._col.objects.unlink(force.obj)
        self._forces.remove(force)

    @property
    def forces(self):
        """List of forces in this collection."""
        return self._forces

    @staticmethod
    def global_forces() -> "ForceCollection":
        """Collection of forces that affects all cells.

        Returns:
            The collection containing global forces.
        """
        if ForceCollection._global_forces is None:
            ForceCollection._global_forces = ForceCollection("globals")
            bpy.context.scene.collection.children.link(
                ForceCollection._global_forces.collection
            )
        return ForceCollection._global_forces


class Boundary(BlenderObject):
    """A boundary for cells."""

    def setup_physics(self):
        """Set up physics for the boundary."""
        BoundaryCollisionConstructor().construct(self.obj)


def create_boundary(loc: tuple, size: float, mesh: str = "icosphere"):
    """Create a boundary.

    Args:
        loc: Center of the boundary.
        size: Radius of the boundary.
        mesh: Shape fo the boundary.
    """

    obj = create_mesh("Boundary", loc, mesh=mesh, size=size)
    bpy.context.scene.collection.objects.link(obj)

    boundary = Boundary(obj)
    boundary.setup_physics()
    boundary.hide()
    return boundary
