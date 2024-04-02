import numpy as np
import bpy, bmesh
from mathutils import Vector, Matrix, Euler, Quaternion
from goo.force import BlenderObject, make_homo_adhesion


# TODO: change so that origin updates to cell center
class Cell(BlenderObject):
    def __init__(self, obj, celltype=None):
        super(Cell, self).__init__(obj)
        self._celltype = celltype
        self._last_division_time = 0
        self.adhesion_forces = {}
        self.cloth_settings = {}
        self.cloth_collision_settings = {}
        self.collision_settings = {}
        self.physics_enabled = False

    def copy(self):
        """Copies cell data and physics modifiers (except forces)."""
        obj_copy = self.obj.copy()
        obj_copy.data = self.obj.data.copy()
        bpy.context.collection.objects.link(obj_copy)

        cell_copy = Cell(obj_copy)
        # if self.cloth_mod is not None:
        #     cell_copy.setup_cloth(settings=self.cloth_mod.settings)
        # if self.collision_mod is not None:
        #     cell_copy.setup_cloth(settings=self.collision_mod.settings)
        return cell_copy

    @property
    def celltype(self):
        return self._celltype

    @celltype.setter
    def celltype(self, celltype):
        self._celltype = celltype
        self.obj.data["celltype"] = celltype.name

    @property
    def last_division_time(self):
        return self._last_division_time

    @last_division_time.setter
    def last_division_time(self, time):
        self._last_division_time = time
        self.obj.data["last_division_time"] = time

    @property
    def obj_eval(self):
        # return self.obj
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
        """Calculates the volume of the Blender mesh.

        In order to retrieve the mesh as it is currrently evaluated - including
        the effect of all modifiers - in the simulation, its corresponding evaluated
        ID is obtained from the dependency graph for the current context.

        .. seealso:: [Blender API Documentation >
                    ``evaluated_get(depsgraph)``]
                    (https://docs.blender.org/api/current/bpy.types.ID.html?highlight=evaluated_get#bpy.types.ID.evaluated_get)

        :param bpy.data.objects['name'] obj: The Blender mesh.
        :returns: The volume of the mesh.
        :rtype: float

        .. note:: The function may return a negative volume.
                    See ``calc_volume(signed=True)``.

        """
        obj_eval = self.obj_eval

        # Use BMesh to calculate volume of a mesh
        mesh_from_eval = obj_eval.to_mesh()
        bm = bmesh.new()
        bm.from_mesh(mesh_from_eval)
        bm.transform(obj_eval.matrix_world)
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

        vert_indices = np.argsort(
            np.asarray(verts).dot(axis)
        )  # sort indices by distance
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

    def has_modifier(self, type):
        return len([m for m in self.obj.modifiers if m.type == type]) > 0

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
        if self.cloth_settings and self.cloth_collision_settings:
            _update_cloth_settings(
                cloth_mod, self.cloth_settings, self.cloth_collision_settings
            )
            self.cloth_settings = {}
            self.cloth_collision_settings = {}
        else:
            _setup_cloth_defaults(cloth_mod)

    def setup_collision(self):
        collision_mod = self.obj.modifiers.new(name="Collision", type="COLLISION")
        if self.collision_settings:
            _update_collision_settings(collision_mod, self.collision_settings)
            self.collision_settings = {}
        else:
            _setup_collision_defaults(collision_mod)

    def disable_physics(self, cloth=True, collision=True, forces=True):
        if cloth and self.cloth_mod is not None:
            _store_cloth_settings(
                self.cloth_mod, self.cloth_settings, self.cloth_collision_settings
            )
            self.obj.modifiers.remove(self.cloth_mod)
        if collision and self.collision_mod is not None:
            _store_collision_settings(self.collision_mod, self.collision_settings)
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
    def add_force(self, celltype, force):
        """Add force towards a specific cell type (homotypic = same cell type)."""
        self.adhesion_forces[celltype] = force

    def get_force(self, celltype):
        """Find force towards a specific cell type."""
        return self.adhesion_forces[celltype]

    @property
    def forces(self):
        return list(self.adhesion_forces.values())


def create_cell(name, loc, physics_on=True, smooth=True, **kwargs):
    obj = _create_mesh(name, loc, mesh="icosphere", **kwargs)
    cell = Cell(obj)
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
    # use of bmesh is about 10x faster than bpy.ops
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

    bpy.context.collection.objects.link(obj)
    return obj


def _create_division_plane(name, major_axis, com):
    """
    Creates a plane orthogonal to the long axis vector
    and passing through the cell's center of mass.
    """
    # TODO: not sure about length of the plane, maybe short axis?
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

    # Hide plane
    plane.hide_set(True)
    return plane


cloth_setting_attrs = [
    "quality",
    "air_damping",
    "bending_model",
    "mass",
    "time_scale",
    "tension_stiffness",
    "compression_stiffness",
    "shear_stiffness",
    "bending_stiffness",
    "tension_damping",
    "compression_damping",
    "shear_damping",
    "bending_damping",
    "use_pressure",
    "uniform_pressure_force",
    "use_pressure_volume",
    "target_volume",
    "pressure_factor",
    "fluid_density",
]
cloth_collision_setting_attrs = [
    "collision_quality",
    "use_collision",
    "use_self_collision",
    "self_friction",
    "friction",
    "self_distance_min",
    "distance_min",
    "self_impulse_clamp",
]


def _setup_cloth_defaults(cloth_mod, stiffness=15, pressure=5):
    cloth_mod.settings.quality = 10
    cloth_mod.settings.air_damping = 10
    cloth_mod.settings.bending_model = "ANGULAR"
    cloth_mod.settings.mass = 1
    cloth_mod.settings.time_scale = 2

    # Cloth > Stiffness
    cloth_mod.settings.tension_stiffness = stiffness
    cloth_mod.settings.compression_stiffness = stiffness
    cloth_mod.settings.shear_stiffness = stiffness
    cloth_mod.settings.bending_stiffness = 1
    # Cloth > Damping
    cloth_mod.settings.tension_damping = 50
    cloth_mod.settings.compression_damping = 50
    cloth_mod.settings.shear_damping = 50
    cloth_mod.settings.bending_damping = 0.5
    # Cloth > Pressure
    cloth_mod.settings.use_pressure = True
    cloth_mod.settings.uniform_pressure_force = pressure
    cloth_mod.settings.use_pressure_volume = True
    cloth_mod.settings.target_volume = 1
    cloth_mod.settings.pressure_factor = 2
    cloth_mod.settings.fluid_density = 1.05
    # Cloth > Collisions
    cloth_mod.collision_settings.collision_quality = 5
    cloth_mod.collision_settings.use_collision = True
    cloth_mod.collision_settings.use_self_collision = True
    cloth_mod.collision_settings.self_friction = 0
    cloth_mod.collision_settings.friction = 0
    cloth_mod.collision_settings.self_distance_min = 0.005
    cloth_mod.collision_settings.distance_min = 0.005
    cloth_mod.collision_settings.self_impulse_clamp = 0


def _store_cloth_settings(cloth_mod, cloth_settings, cloth_collision_settings):
    for attr in cloth_setting_attrs:
        cloth_settings[attr] = getattr(cloth_mod.settings, attr)
    for attr in cloth_collision_setting_attrs:
        cloth_collision_settings[attr] = getattr(cloth_mod.collision_settings, attr)


def _update_cloth_settings(cloth_mod, cloth_settings, cloth_collision_settings):
    for attr in cloth_setting_attrs:
        setattr(cloth_mod.settings, attr, cloth_settings[attr])
    for attr in cloth_collision_setting_attrs:
        setattr(cloth_mod.collision_settings, attr, cloth_collision_settings[attr])


collision_setting_attrs = [
    "use_culling",
    "damping",
    "thickness_outer",
    "thickness_inner",
    "cloth_friction",
]


def _setup_collision_defaults(collision_mod):
    collision_mod.settings.use_culling = True
    collision_mod.settings.damping = 1
    collision_mod.settings.thickness_outer = 0.025
    collision_mod.settings.thickness_inner = 0.25
    collision_mod.settings.cloth_friction = 0


def _store_collision_settings(collision_mod, collision_settings):
    for attr in collision_setting_attrs:
        collision_settings[attr] = getattr(collision_mod.settings, attr)


def _update_collision_settings(collision_mod, collision_settings):
    for attr in collision_setting_attrs:
        setattr(collision_mod.settings, attr, collision_settings[attr])


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
    def __init__(self, name, physics_on=True):
        self.name = name
        self._cells = set()

        self.physics_on = physics_on
        self.homo_adhesion_strength = 2000

    def add_cell(self, cell):
        cell.celltype = self
        self._cells.add(cell)
        if self.physics_on:
            homo_adhesion = make_homo_adhesion(cell.obj, self._homo_adhesion_strength)
            cell.add_force(self.name, homo_adhesion)

    def create_cell(self, name, loc, **kwargs):
        cell = create_cell(name, loc, physics_on=self.physics_on, **kwargs)
        if self.physics_on:
            homo_adhesion = make_homo_adhesion(cell.obj, self._homo_adhesion_strength)
            cell.add_force(self.name, homo_adhesion)
        self.add_cell(cell)
        return cell

    def create_cells(celltype, locs=None, shape=None, size=None):
        # TODO: create cells in shape or in specified locations of specified cell types
        pass

    @property
    def cells(self):
        return list(self._cells)

    @property
    def homo_adhesion_strength(self):
        return self._homo_adhesion_strength

    @homo_adhesion_strength.setter
    def homo_adhesion_strength(self, strength):
        self._homo_adhesion_strength = strength
        if not self.physics_on:
            return
        for cell in self._cells:
            homo_adhesion = cell.get_force(self.name)
            if homo_adhesion:
                homo_adhesion.strength = strength
            else:
                homo_adhesion = make_homo_adhesion(
                    cell.obj, self._homo_adhesion_strength
                )
                cell.add_force(self.name, homo_adhesion)
