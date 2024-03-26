import numpy as np
import bpy, bmesh
from mathutils import Vector, Matrix, Euler


class Cell:
    def __init__(self, obj, celltype=None):
        self.obj = obj
        self.celltype = celltype

    @property
    def name(self):
        return self.obj.name

    @name.setter
    def name(self, name):
        self.obj.name = name

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
        # print(f"Volume of {obj.name}: {abs(volume)}")  # Output the result
        bm.free()

        return volume

    def COM(self):
        """Calculates the center of mass of a mesh.

        This function fetch the evaluated object's dependency graph,
        retrieves the coordinates of each vertex then computes the center of mass
        as the mean position among the set of vertices.

        :param bpy.data.objects['name'] obj: The Blender mesh.
        :returns: The coordinates of the center of mass of the mesh as a tuple(x, y, z).
        :rtype: tuple
        """
        obj_eval = self.obj_eval
        vert_coords = _get_vertex_coords(obj_eval)
        com = np.mean(vert_coords, axis=0)
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
        return _create_division_plane(
            self.obj.name,
            self.major_axis(),
            self.COM(),
        )

    def divide(self, dividerHandler):
        return dividerHandler.make_divide(self)


def create_cell(name, loc, **kwargs):
    obj = _create_mesh(name, loc, mesh="icosphere", **kwargs)
    cell = Cell(obj)
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

    if mesh == "icosphere":
        matrix_world = Matrix.LocRotScale(loc, rotation, scale)
        bmesh.ops.create_icosphere(
            bm, subdivisions=subdivisions, radius=size, matrix=matrix_world, **kwargs
        )
    elif mesh == "plane":
        matrix_world = Matrix.LocRotScale(loc, rotation, scale)
        bmesh.ops.create_grid(
            bm, x_segments=1, y_segments=1, size=size / 2, matrix=matrix_world, **kwargs
        )
    elif mesh == "monkey":
        matrix_world = Matrix.LocRotScale(loc, rotation, size * Vector(scale))
        bmesh.ops.create_monkey(bm, matrix=matrix_world, **kwargs)
    else:
        raise ValueError("""mesh must be one of "icosphere", "plane", or "monkey".""")

    me = bpy.data.meshes.new(f"{name}_mesh")
    bm.to_mesh(me)
    bm.free()

    obj = bpy.data.objects.new(name, me)
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
    # Sort the eigenvectors by descending eigenvalues
    eigenvectors = eigenvectors[:, eigenvalues.argsort()[::-1]]
    return eigenvectors


def _create_division_plane(name, long_axis, com):
    """
    Creates a plane orthogonal to the long axis vector
    and passing through the cell's center of mass.

    :param long_axis: The long axis vector.
    :type long_axis: numpy.ndarray
    :param com: The center of mass of the mesh.
    :type com: numpy.ndarray
    """
    # Define a new plane object
    plane = _create_mesh(
        f"{name}_division_plane",
        loc=com,
        mesh="plane",
        size=long_axis.axis().length + 1,
        rotation=long_axis.axis().to_track_quat("Z", "Y"),
    )

    # Add solidify modifier to the plane, add thickness to the plane
    solid_mod = plane.modifiers.new(name="Solidify", type="SOLIDIFY")
    solid_mod.offset = 0
    solid_mod.thickness = 0.025

    # Hide the plane
    plane.hide_set(True)
    return plane


# --- End: Blender Functions ---


class Axis:
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
    def __init__(self):
        pass

    def add_cell(self, cell):
        cell.celltype = self
        self.cells.append(cell)

    def create_cell(celltype, loc):
        cell = Cell(loc)
        celltype.add_cell(cell)

    def create_cells(celltype, locs=None, shape=None, size=None):
        # create cells in shape or in specified locations of specified cell types
        pass
