import numpy as np
import bpy, bmesh
from mathutils import Vector, Matrix, Euler, Quaternion


class Cell:
    def __init__(self, obj, celltype=None):
        self.obj = obj
        self._celltype = celltype
        self._last_division_time = 0

    @property
    def name(self):
        return self.obj.name

    @name.setter
    def name(self, name):
        self.obj.name = name

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
        return self.obj
        # TODO: this is expensive! make it so that evaluated_depsgraph_get is called just once per cycle!
        # dg = bpy.context.evaluated_depsgraph_get()
        # obj_eval = self.obj.evaluated_get(dg)
        # return obj_eval

    def volume(self):
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

    def COM(self, global_coords=True):
        """Calculates the center of mass of a mesh."""
        obj_eval = self.obj_eval
        vert_coords = _get_vertex_coords(obj_eval)
        com = Vector(np.mean(vert_coords, axis=0))
        if not global_coords:
            com = self.obj.matrix_world.inverted() @ com
        return Vector(com)

    def _eigenvector(self, n):
        """Returns the nth eigenvector in object space as a line defined by two Vectors."""
        obj_eval = self.obj_eval
        vert_coords = _get_vertex_coords(obj_eval)
        eigenvectors = _eigenvectors(vert_coords)

        axis = eigenvectors[:, n]

        vertices = obj_eval.data.vertices
        vert_indices = np.argsort(vert_coords.dot(axis))  # sort indices by distance
        first_vertex = vertices[vert_indices[0]].co
        last_vertex = vertices[vert_indices[-1]].co

        return Axis(Vector(axis), first_vertex, last_vertex, obj_eval.matrix_world)

    def minor_axis(self):
        return self._eigenvector(1)

    def major_axis(self):
        return self._eigenvector(0)

    def create_division_plane(self):
        return _create_division_plane(self.obj.name, self.major_axis(), self.COM())

    def divide(self, divisionLogic):
        mother, daughter = divisionLogic.make_divide(self)
        if mother.celltype:
            mother.celltype.add_cell(daughter)
        return mother, daughter

    def remesh(self):
        # use of object ops is 2x faster than remeshing with modifiers
        self.obj.data.remesh_mode = "VOXEL"
        self.obj.data.remesh_voxel_size = 0.25
        with bpy.context.temp_override(active_object=self.obj):
            bpy.ops.object.voxel_remesh()

    def copy(self):
        obj_copy = self.obj.copy()
        obj_copy.data = self.obj.data.copy()
        bpy.context.collection.objects.link(obj_copy)
        return Cell(obj_copy)


def create_cell(name, loc, **kwargs):
    obj = _create_mesh(name, loc, mesh="icosphere", **kwargs)
    cell = Cell(obj)
    cell.remesh()
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


def _get_vertex_coords(obj_eval):
    """Returns a list of mesh vertex coordinates in global coordinates."""
    vertices = obj_eval.data.vertices
    vert_coords = np.asarray([obj_eval.matrix_world @ v.co for v in vertices])
    return vert_coords


def _eigenvectors(vert_coords):
    """Returns a list eigenvectors in object space sorted by descending eigenvalue."""
    # Calculate the covariance matrix of the vertices
    covariance_matrix = np.cov(vert_coords, rowvar=False)
    # Calculate the eigenvectors and eigenvalues of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
    eigenvectors = eigenvectors[:, eigenvalues.argsort()[::-1]]
    return eigenvectors


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
        size=major_axis.length(global_coords=True) + 1,
        rotation=major_axis.axis().to_track_quat("Z", "Y"),
    )

    # Add thickness to plane
    solid_mod = plane.modifiers.new(name="Solidify", type="SOLIDIFY")
    solid_mod.offset = 0
    solid_mod.thickness = 0.025

    # Hide plane
    plane.hide_set(True)
    return plane


# --- End: Blender Functions ---


class Axis:
    # TODO: global vs. local coordinate parameters are probably wrong
    def __init__(self, axis, start, end, world_matrix):
        """
        :param axis: Vector of axis
        :param start: Vector of start endpoint in object space
        :param end: Vector of end endpoint in object space
        :param world_matrix: 4x4 matrix of object to world transformation
        """
        self._axis = axis
        self._start = start
        self._end = end
        self._matrix_world = world_matrix

    def axis(self, global_coords=False):
        if global_coords:
            return self._matrix_world @ self._axis
        return self._axis

    def endpoints(self, global_coords=False):
        if global_coords:
            return [self._matrix_world @ self._start, self._matrix_world @ self._end]
        return [self._start, self._end]

    def length(self, global_coords=False):
        start, end = self.endpoints(global_coords)
        return (end - start).length


class CellType:
    def __init__(self, name):
        self.name = name
        self._cells = set()

    def add_cell(self, cell):
        cell.celltype = self
        self._cells.add(cell)

    def create_cell(self, name, loc, **kwargs):
        cell = create_cell(name, loc, **kwargs)
        self.add_cell(cell)
        return cell

    def create_cells(celltype, locs=None, shape=None, size=None):
        # TODO: create cells in shape or in specified locations of specified cell types
        pass

    @property
    def cells(self):
        return list(self._cells)
