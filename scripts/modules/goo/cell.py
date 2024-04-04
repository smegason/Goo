from functools import reduce
from typing import Optional

import numpy as np

import bpy, bmesh
from mathutils import *
from goo.force import *
from goo.utils import *


class Cell(BlenderObject):
    def __init__(self, obj: bpy.types.Object):
        super(Cell, self).__init__(obj)
        self._effectors = bpy.data.collections.new(f"{obj.name}_effectors")
        bpy.context.scene.collection.children.link(self._effectors)

        self.celltype: CellType = None

        self._physics_enabled = False
        self.cloth_settings = {}
        self.collision_settings = {}

        self._adhesion_forces: list[AdhesionForce] = []
        self._motion_force: Force = None

        self._last_division_time = 0

    @property
    def name(self) -> str:
        return self._obj.name

    @name.setter
    def name(self, name: str):
        self._obj.name = name
        self._effectors.name = f"{name}_effectors"

    def copy(self) -> "Cell":
        """Copies cell data and physics modifiers, including homotypic adhesion forces."""
        obj_copy = self._obj.copy()
        obj_copy.data = self._obj.data.copy()
        bpy.context.scene.collection.objects.link(obj_copy)

        cell_copy = Cell(obj_copy)
        return cell_copy

    # TODO: turn into custom properties
    @property
    def last_division_time(self):
        return self._last_division_time

    @last_division_time.setter
    def last_division_time(self, time):
        self._last_division_time = time
        self._obj.data["last_division_time"] = time

    @property
    def obj_eval(self):
        dg = bpy.context.evaluated_depsgraph_get()
        obj_eval = self._obj.evaluated_get(dg)
        return obj_eval

    def get_vertices(self, local_coords=False):
        verts = self.obj_eval.data.vertices
        if local_coords:
            return [v.co for v in verts]
        else:
            return [self.obj_eval.matrix_world @ v.co for v in verts]

    def get_volume(self) -> float:
        """Calculates the volume of the Blender mesh."""
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())
        bm.transform(self.obj_eval.matrix_world)
        volume = bm.calc_volume()
        bm.free()

        return volume

    def get_COM(self, local_coords: bool = False) -> Vector:
        """Calculates the center of mass of a mesh."""
        vert_coords = self.get_vertices(local_coords)
        com = Vector(np.mean(vert_coords, axis=0))
        return com

    def _get_eigenvector(self, n):
        """Returns the nth eigenvector (axis) in object space as a line defined by two Vectors."""
        obj_eval = self.obj_eval
        verts = self.get_vertices()

        # Calculate the eigenvectors and eigenvalues of the covariance matrix
        covariance_matrix = np.cov(verts, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
        eigenvectors = eigenvectors[:, eigenvalues.argsort()[::-1]]

        axis = eigenvectors[:, n]

        # sort indices by distance
        vert_indices = np.argsort(np.asarray(verts).dot(axis))
        first_vertex = verts[vert_indices[0]]
        last_vertex = verts[vert_indices[-1]]

        return Axis(Vector(axis), first_vertex, last_vertex, obj_eval.matrix_world)

    def get_major_axis(self) -> "Axis":
        return self._get_eigenvector(0)

    def get_minor_axis(self) -> "Axis":
        return self._get_eigenvector(1)

    def recenter(self):
        """Recenter origin to COM."""
        com = self.get_COM()
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())
        bmesh.ops.translate(bm, verts=bm.verts, vec=-self.get_COM(local_coords=True))
        bm.to_mesh(self._obj.data)
        bm.free()

        self.loc = com

    def divide(self, divisionLogic):
        mother, daughter = divisionLogic.make_divide(self)
        if mother.celltype:
            mother.celltype.add_cell(daughter)
        return mother, daughter

    def remesh(self, smooth: bool = True):
        # use of object ops is 2x faster than remeshing with modifiers
        self._obj.data.remesh_mode = "VOXEL"
        self._obj.data.remesh_voxel_size = 0.25
        with bpy.context.temp_override(active_object=self._obj):
            bpy.ops.object.voxel_remesh()

        for f in self._obj.data.polygons:
            f.use_smooth = smooth

    # ----- PHYSICS -----
    def physics_enabled(self) -> bool:
        return self._physics_enabled

    def enable_physics(self, cloth=True, collision=True, forces=True):
        if cloth and self.cloth_mod is None:
            self.setup_cloth()
        if collision and self.collision_mod is None:
            self.setup_collision()
        if forces:
            for force in self.adhesion_forces:
                force.enable()
        self._physics_enabled = True

    def disable_physics(self, cloth=True, collision=True, forces=True):
        if cloth and self.cloth_mod is not None:
            store_settings(self.cloth_mod, default_cloth_settings, self.cloth_settings)
            self._obj.modifiers.remove(self.cloth_mod)
        if collision and self.collision_mod is not None:
            store_settings(
                self.collision_mod, default_collision_settings, self.collision_settings
            )
            self._obj.modifiers.remove(self.collision_mod)
        if forces:
            for force in self.adhesion_forces:
                force.disable()
        self._physics_enabled = False

    def setup_cloth(self):
        cloth_mod = self._obj.modifiers.new(name="Cloth", type="CLOTH")
        update_settings(cloth_mod, default_cloth_settings, self.cloth_settings)
        self.cloth_mod.settings.effector_weights.collection = self._effectors
        self.cloth_settings = {}

    def setup_collision(self):
        collision_mod = self._obj.modifiers.new(name="Collision", type="COLLISION")
        update_settings(
            collision_mod, default_collision_settings, self.collision_settings
        )
        self.collision_settings = {}

    def get_modifier(self, type) -> Optional[bpy.types.Modifier]:
        return next((m for m in self._obj.modifiers if m.type == type), None)

    @property
    def cloth_mod(self) -> Optional[bpy.types.ClothModifier]:
        return self.get_modifier("CLOTH")

    @property
    def collision_mod(self) -> Optional[bpy.types.CollisionModifier]:
        return self.get_modifier("COLLISION")

    @property
    def stiffness(self) -> float:
        return self.cloth_mod.settings.tension_stiffness

    @stiffness.setter
    def stiffness(self, stiffness: float):
        self.cloth_mod.settings.tension_stiffness = stiffness
        self.cloth_mod.settings.compression_stiffness = stiffness
        self.cloth_mod.settings.shear_stiffness = stiffness

    @property
    def pressure(self) -> float:
        return self.cloth_mod.settings.uniform_pressure_force

    @pressure.setter
    def pressure(self, pressure: float):
        self.cloth_mod.settings.uniform_pressure_force = pressure

    # ----- FORCES -----
    def link_adhesion_force(self, force: AdhesionForce):
        """Add force that stems from this cell."""
        self._adhesion_forces.append(force)

    def add_effector(self, force: Force | ForceCollection):
        """Add a force or a collection of forces that affect this cell."""
        if isinstance(force, Force):
            self._effectors.objects.link(force.obj)
        elif isinstance(force, ForceCollection):
            self._effectors.children.link(force.collection)

    def remove_effector(self, force: Force | ForceCollection):
        """Removes a force or a collection of forces that affect this cell."""
        if isinstance(force, Force):
            self._effectors.objects.unlink(force.obj)
        elif isinstance(force, ForceCollection):
            self._effectors.children.unlink(force.collection)

    @property
    def adhesion_forces(self) -> list[AdhesionForce]:
        return self._adhesion_forces

    @property
    def motion_force(self) -> Force:
        return self._motion_force

    @motion_force.setter
    def motion_force(self, force: Force):
        if self.motion_force:
            self.remove_effector(self.motion_force)
        self.add_effector(force)
        self._motion_force = force


def create_cell(name, loc, physics_on=True, **kwargs) -> Cell:
    obj = create_mesh(name, loc, mesh="icosphere", **kwargs)
    bpy.context.scene.collection.objects.link(obj)
    cell = Cell(obj)
    cell.remesh()

    if physics_on:
        cell.enable_physics()
    return cell


# --- Physics Modifier Utilities ---
default_stiffness = 15
default_pressure = 5
default_cloth_settings = {
    "settings.quality": 10,
    "settings.air_damping": 10,
    "settings.bending_model": "ANGULAR",
    "settings.mass": 1,
    "settings.time_scale": 1,
    # Cloth > Stiffness
    "settings.tension_stiffness": default_stiffness,
    "settings.compression_stiffness": default_stiffness,
    "settings.shear_stiffness": default_stiffness,
    "settings.bending_stiffness": 1,
    # Cloth > Damping
    "settings.tension_damping": 50,
    "settings.compression_damping": 50,
    "settings.shear_damping": 50,
    "settings.bending_damping": 0.5,
    # Cloth > Pressure
    "settings.use_pressure": True,
    "settings.uniform_pressure_force": default_pressure,
    "settings.use_pressure_volume": True,
    "settings.target_volume": 1,
    "settings.pressure_factor": 2,
    "settings.fluid_density": 1.05,
    # Cloth > Field Weights
    "settings.effector_weights.collection": None,
    # Cloth > Collisions
    "collision_settings.collision_quality": 5,
    "collision_settings.use_collision": True,
    "collision_settings.use_self_collision": True,
    "collision_settings.self_friction": 0,
    "collision_settings.friction": 0,
    "collision_settings.self_distance_min": 0.005,
    "collision_settings.distance_min": 0.005,
    "collision_settings.self_impulse_clamp": 0,
}

default_collision_settings = {
    "settings.use_culling": True,
    "settings.damping": 1,
    "settings.thickness_outer": 0.025,
    "settings.thickness_inner": 0.25,
    "settings.cloth_friction": 0,
}


def store_settings(mod, default_settings, custom_settings):
    for k in default_settings.keys():
        v = rgetattr(mod, k)
        custom_settings[k] = v


def update_settings(mod, default_settings, custom_settings=None):
    for k in default_settings.keys():
        v = custom_settings[k] if custom_settings else default_settings[k]
        rsetattr(mod, k, v)


class Axis:
    def __init__(self, axis, start, end, world_matrix):
        """
        :param axis: Vector of axis
        :param start: Vector of start endpoint in object space
        :param end: Vector of end endpoint in object space
        :param world_matrix: 4x4 matrix of object to world transformation
        :param local_coords: if coordinates given are in local space
        """
        self._axis = axis
        self._start = start
        self._end = end
        self._matrix_world = world_matrix.inverted()

    def axis(self, local_coords=False):
        if local_coords:
            axis = self._axis.copy()
            axis.rotate(self._matrix_world.to_quaternion())
            return axis
        return self._axis

    def endpoints(self, local_coords=False):
        mat = self._matrix_world if local_coords else Matrix.Identity(4)
        return [mat @ self._start, mat @ self._end]

    def length(self, local_coords=False):
        start, end = self.endpoints(local_coords)
        return (end - start).length


class CellType:
    def __init__(self, collection: ForceCollection, physics_on=True):
        self.homo_adhesions = collection
        self._cells = set()

        self.physics_enabled = physics_on
        self.hetero_adhesions = {}

        self.homo_adhesion_strength = 2000
        self.motion_strength = 0

    @property
    def name(self):
        return self.homo_adhesions.name

    def add_cell(self, cell: Cell):
        cell.celltype = self
        self._cells.add(cell)

        if not self.physics_enabled:
            return
        # add homotypic adhesion force
        homo_adhesion = create_adhesion(self.homo_adhesion_strength, obj=cell._obj)
        cell.link_adhesion_force(homo_adhesion)
        self.homo_adhesions.add_force(homo_adhesion)

        cell.add_effector(self.homo_adhesions)

        # add hetero adhesion forces
        for celltype, item in self.hetero_adhesions.items():
            outgoing_forces, incoming_forces, strength = item
            hetero_adhesion = create_adhesion(
                strength, name=cell.name + "_to_" + celltype.name, loc=cell.loc
            )
            outgoing_forces.add_force(hetero_adhesion)
            cell.link_adhesion_force(hetero_adhesion)

            cell.add_effector(incoming_forces)

        # add motion force
        motion = create_motion(
            name=cell.name + "_motion",
            loc=cell.get_COM(),
            strength=self.motion_strength,
        )
        cell.motion_force = motion

    def create_cell(self, name, loc, **kwargs) -> Cell:
        cell = create_cell(name, loc, physics_on=self.physics_enabled, **kwargs)
        self.add_cell(cell)
        return cell

    # TODO: create cells in shape or in specified locations of specified cell types
    def create_cells(celltype, locs=None, shape=None, size=None):
        pass

    @property
    def cells(self) -> list[Cell]:
        return list(self._cells)

    def set_homo_adhesion(self, strength: float):
        self.homo_adhesion_strength = strength

    def set_hetero_adhesion(self, other_celltype: "CellType", strength: float):
        outgoing_forces = ForceCollection(self.name + "_to_" + other_celltype.name)
        incoming_forces = ForceCollection(other_celltype.name + "_to_" + self.name)
        self.hetero_adhesions[other_celltype] = (
            outgoing_forces,
            incoming_forces,
            strength,
        )
        other_celltype.hetero_adhesions[self] = (
            incoming_forces,
            outgoing_forces,
            strength,
        )

    def set_motion(self, strength: float):
        self.motion_strength = strength


def create_celltype(name, physics_on=True) -> CellType:
    collection = ForceCollection(name)
    return CellType(collection, physics_on)
