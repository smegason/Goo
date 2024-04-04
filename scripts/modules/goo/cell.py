from functools import reduce
import numpy as np

import bpy, bmesh
from mathutils import Vector, Matrix, Euler, Quaternion

from goo.force import BlenderObject, create_force, create_adhesion


class Cell(BlenderObject):
    def __init__(self, obj: bpy.types.Object, obj_col: bpy.types.Collection):
        super(Cell, self).__init__(obj)
        self.obj_col = obj_col
        self.celltype = None

        self.physics_enabled = False
        self.cloth_settings = {}
        self.collision_settings = {}
        self.adhesion_forces = []

        self._last_division_time = 0

    @property
    def name(self):
        return self.obj.name

    @name.setter
    def name(self, name):
        self.obj.name = name
        self.obj_col.name = name

    # TODO: fix
    def copy(self):
        """Copies cell data and physics modifiers (except forces)."""
        obj_copy = self.obj.copy()
        obj_copy.data = self.obj.data.copy()

        obj_col = bpy.data.collections.new(self.name)
        bpy.context.scene.collection.children.link(obj_col)
        obj_col.objects.link(obj_copy)

        cell_copy = Cell(obj_copy, obj_col)
        if self.celltype:
            self.celltype.add_cell(cell_copy)
        else:
            bpy.context.collection.objects.link(obj_copy)

        return cell_copy

    @property
    def last_division_time(self):
        return self._last_division_time

    @last_division_time.setter
    def last_division_time(self, time):
        self._last_division_time = time
        self.obj.data["last_division_time"] = time

    @property
    def obj_eval(self):
        dg = bpy.context.evaluated_depsgraph_get()
        obj_eval = self.obj.evaluated_get(dg)
        return obj_eval

    def get_vertices(self, local_coords=False):
        verts = self.obj_eval.data.vertices
        if local_coords:
            return [v.co for v in verts]
        else:
            return [self.obj_eval.matrix_world @ v.co for v in verts]

    def get_volume(self):
        """Calculates the volume of the Blender mesh."""
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())
        bm.transform(self.obj_eval.matrix_world)
        volume = bm.calc_volume()
        bm.free()

        return volume

    def get_COM(self, local_coords=False):
        """Calculates the center of mass of a mesh."""
        vert_coords = self.get_vertices(local_coords)
        com = Vector(np.mean(vert_coords, axis=0))
        return com

    def recenter(self):
        """Recenter origin to COM."""
        com = self.get_COM()
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())
        bmesh.ops.translate(bm, verts=bm.verts, vec=-self.get_COM(local_coords=True))
        bm.to_mesh(self.obj.data)
        bm.free()

        self.loc = com

    def _get_eigenvectors(self):
        """Returns a list eigenvectors in object space sorted by descending eigenvalue."""
        # Calculate the covariance matrix of the vertices
        covariance_matrix = np.cov(self.get_vertices(), rowvar=False)
        # Calculate the eigenvectors and eigenvalues of the covariance matrix
        eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
        eigenvectors = eigenvectors[:, eigenvalues.argsort()[::-1]]
        return eigenvectors

    def _get_eigenvector(self, n):
        """Returns the nth eigenvector in object space as a line defined by two Vectors."""
        obj_eval = self.obj_eval
        verts = self.get_vertices()
        eigenvectors = self._get_eigenvectors()
        axis = eigenvectors[:, n]

        # sort indices by distance
        vert_indices = np.argsort(np.asarray(verts).dot(axis))
        first_vertex = verts[vert_indices[0]]
        last_vertex = verts[vert_indices[-1]]

        return Axis(Vector(axis), first_vertex, last_vertex, obj_eval.matrix_world)

    def get_minor_axis(self):
        return self._get_eigenvector(1)

    def get_major_axis(self):
        return self._get_eigenvector(0)

    def create_division_plane(self):
        return _create_division_plane(
            self.obj.name, self.get_major_axis(), self.get_COM()
        )

    def divide(self, divisionLogic):
        mother, daughter = divisionLogic.make_divide(self)
        if mother.celltype:
            mother.celltype.add_cell(daughter)
        return mother, daughter

    def remesh(self, smooth=True):
        # use of object ops is 2x faster than remeshing with modifiers
        self.obj.data.remesh_mode = "VOXEL"
        self.obj.data.remesh_voxel_size = 0.25
        with bpy.context.temp_override(active_object=self.obj):
            bpy.ops.object.voxel_remesh()

        for f in self.obj.data.polygons:
            f.use_smooth = smooth

    def get_modifier(self, type):
        return next((m for m in self.obj.modifiers if m.type == type), None)

    # ----- PHYSICS -----
    def enable_physics(self, cloth=True, collision=True, forces=True):
        if cloth and self.cloth_mod is None:
            self.setup_cloth()
        if collision and self.collision_mod is None:
            self.setup_collision()
        if forces:
            for force in self.forces:
                force.enable()
        self.physics_enabled = True

    def setup_cloth(self):
        cloth_mod = self.obj.modifiers.new(name="Cloth", type="CLOTH")
        update_settings(cloth_mod, default_cloth_settings, self.cloth_settings)
        self.cloth_mod.settings.effector_weights.collection = self.obj_col
        self.cloth_settings = {}

    def setup_collision(self):
        collision_mod = self.obj.modifiers.new(name="Collision", type="COLLISION")
        update_settings(
            collision_mod, default_collision_settings, self.collision_settings
        )
        self.collision_settings = {}

    def disable_physics(self, cloth=True, collision=True, forces=True):
        if cloth and self.cloth_mod is not None:
            store_settings(self.cloth_mod, default_cloth_settings, self.cloth_settings)
            self.obj.modifiers.remove(self.cloth_mod)
        if collision and self.collision_mod is not None:
            store_settings(
                self.collision_mod, default_collision_settings, self.collision_settings
            )
            self.obj.modifiers.remove(self.collision_mod)
        if forces:
            for force in self.forces:
                force.disable()
        self.physics_enabled = False

    def toggle_physics(self, on):
        if self.cloth_mod:
            self.cloth_mod.show_viewport = on
            self.cloth_mod.show_render = on
            self.collision_mod.collision.use = on

    @property
    def cloth_mod(self):
        return next((m for m in self.obj.modifiers if m.type == "CLOTH"), None)

    @property
    def collision_mod(self):
        return next((m for m in self.obj.modifiers if m.type == "COLLISION"), None)

    @property
    def stiffness(self):
        return self.cloth_mod.settings.tension_stiffness

    @stiffness.setter
    def stiffness(self, stiffness):
        self.cloth_mod.settings.tension_stiffness = stiffness
        self.cloth_mod.settings.compression_stiffness = stiffness
        self.cloth_mod.settings.shear_stiffness = stiffness

    @property
    def pressure(self):
        return self.cloth_mod.settings.uniform_pressure_force

    @pressure.setter
    def pressure(self, pressure):
        self.cloth_mod.settings.uniform_pressure_force = pressure

    # ----- FORCES -----
    def add_force(self, force):
        """Add force towards a specific cell type (homotypic = same cell type)."""
        self.adhesion_forces.append(force)

    def get_force(self, celltype):
        """Find force towards a specific cell type."""
        return self.adhesion_forces

    @property
    def forces(self):
        return self.adhesion_forces


def create_cell(name, loc, obj=None, physics_on=True, **kwargs) -> Cell:
    if obj is None:
        obj = _create_mesh(name, loc, mesh="icosphere", **kwargs)
    obj_col = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(obj_col)
    obj_col.objects.link(obj)

    cell = Cell(obj, obj_col)
    cell.remesh()
    if physics_on:
        cell.enable_physics()

    return cell


# --- Begin: Blender Functions ---


def _create_mesh(
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


def _create_division_plane(name, major_axis, com, collection=None):
    """
    Creates a plane orthogonal to the long axis vector
    and passing through the cell's center of mass.
    """
    # Define new plane
    plane = _create_mesh(
        f"{name}_division_plane",
        loc=com,
        mesh="plane",
        size=major_axis.length() + 1,
        rotation=major_axis.axis().to_track_quat("Z", "Y"),
    )

    # Add thickness to plane
    solid_mod = plane.modifiers.new(name="Solidify", type="SOLIDIFY")
    solid_mod.offset = 0
    solid_mod.thickness = 0.025

    plane.hide_set(True)
    return plane


# Utility functions for settings
def rsetattr(obj, attr, val):
    pre, _, post = attr.rpartition(".")
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)


def rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)

    return reduce(_getattr, [obj] + attr.split("."))


default_cloth_settings = {
    "settings.quality": 10,
    "settings.air_damping": 10,
    "settings.bending_model": "ANGULAR",
    "settings.mass": 1,
    "settings.time_scale": 1,
    # Cloth > Stiffness
    "settings.tension_stiffness": 15,
    "settings.compression_stiffness": 15,
    "settings.shear_stiffness": 15,
    "settings.bending_stiffness": 1,
    # Cloth > Damping
    "settings.tension_damping": 50,
    "settings.compression_damping": 50,
    "settings.shear_damping": 50,
    "settings.bending_damping": 0.5,
    # Cloth > Pressure
    "settings.use_pressure": True,
    "settings.uniform_pressure_force": 5,
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
    def __init__(self, collection: bpy.types.Collection, physics_on=True):
        self.collection = collection
        self._cells = set()

        self.physics_enabled = physics_on
        self.homo_adhesion_strength = 2000

        self.hetero_adhesions = {}

        self.motion_strength = 0

    @property
    def name(self):
        return self.collection.name

    def add_cell(self, cell: Cell):
        cell.celltype = self
        self._cells.add(cell)

        if not self.physics_enabled:
            return
        # add homotypic adhesion force
        homo_adhesion = create_adhesion(self.homo_adhesion_strength, obj=cell.obj)
        try:
            self.collection.objects.link(homo_adhesion.obj)
            cell.obj_col.children.link(self.collection)
        except:
            pass
        cell.add_force(homo_adhesion)

        # TODO: hetero adhesion
        for celltype, item in self.hetero_adhesions.items():
            outgoing_col, incoming_col, strength = item
            hetero_adhesion = create_adhesion(
                strength, name=cell.name + "_to_" + celltype.name, loc=cell.loc
            )
            outgoing_col.objects.link(hetero_adhesion.obj)
            cell.add_force(hetero_adhesion)

            cell.obj_col.children.link(incoming_col)

        motion = create_force(
            cell.name + "_motion", loc=cell.get_COM(), strength=self.motion_strength
        )
        cell.obj_col.objects.link(motion.obj)
        cell.add_force(motion)

    def create_cell(self, name, loc, **kwargs) -> Cell:
        cell = create_cell(name, loc, physics_on=self.physics_enabled, **kwargs)
        self.add_cell(cell)
        return cell

    # TODO: create cells in shape or in specified locations of specified cell types
    def create_cells(celltype, locs=None, shape=None, size=None):
        pass

    @property
    def cells(self):
        return list(self._cells)

    def set_homo_adhesion(self, strength):
        self.homo_adhesion_strength = strength

    def set_hetero_adhesion(self, other_celltype, strength):
        outgoing_col = bpy.data.collections.new(
            self.name + "_to_" + other_celltype.name
        )
        incoming_col = bpy.data.collections.new(
            other_celltype.name + "_to_" + self.name
        )
        self.hetero_adhesions[other_celltype] = outgoing_col, incoming_col, strength
        other_celltype.hetero_adhesions[self] = incoming_col, outgoing_col, strength

    def set_motion(self, strength):
        self.motion_strength = strength


def create_celltype(name, physics_on=True) -> CellType:
    col = bpy.data.collections.new(name)
    return CellType(col, physics_on)
