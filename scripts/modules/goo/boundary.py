from typing import Optional
from typing_extensions import override
import bpy
import numpy as np
import math
from goo.utils import *


class Boundary(BlenderObject):
    """A boundary for cells."""

    def __init__(self, obj: bpy.types.Object, mat=None):
        super(Boundary, self).__init__(obj)
        self._mat = mat
        self.obj.data.materials.append(mat)

    def setup_physics(self):
        """Set up physics for the boundary."""
        BoundaryCollisionConstructor().construct(self.obj)
    
    @property
    def size(self) -> int:
        """Size of the boundary."""
        return self.obj.size
    
    @size.setter
    def size(self, size: int):
        self.obj.size = size


class Voxel:
    def __init__(self, obj: bpy.types.Object, loc: tuple, mat=None):
        self.obj = obj
        self.location = loc
        self._mat = mat
        self.obj.data.materials.append(mat)


def create_boundary(loc: tuple, size: float, mesh: str = "icosphere"):
    """Create a boundary.

    Args:
        loc: Center of the boundary.
        size: Radius of the boundary.
        mesh: Shape of the boundary.
    """
    if mesh not in ["icosphere", "cube"]:
        raise ValueError(f"Unsupported mesh type: {mesh}."
                         "Supported types are 'icosphere' and 'cube'.")

    obj = create_mesh("Boundary", loc, mesh=mesh, size=size, subdivisions=4)
    bpy.context.scene.collection.objects.link(obj)

    boundary = Boundary(obj)
    boundary.setup_physics()
    boundary.hide()
    return boundary


def create_grid(loc: tuple, size: tuple, voxel_size: float):
    """Create a 3D grid.

    Args:
        loc: Center of the grid.
        size: Dimensions of the grid (width, height, depth).
        voxel_size: Size of a voxel in the grid.
    """

    if not (isinstance(loc, tuple) and len(loc) == 3):
        raise TypeError("loc must be a tuple of length 3 (x, y, z coordinates)")
    
    if not (isinstance(size, tuple) and len(size) == 3):
        raise TypeError("size must be a tuple of length 3 (width, height, depth)")
    
    min_x = loc[0] - size[0] / 2
    min_y = loc[1] - size[1] / 2
    min_z = loc[2] - size[2] / 2
    max_x = loc[0] + size[0] / 2
    max_y = loc[1] + size[1] / 2
    max_z = loc[2] + size[2] / 2

    voxels = []

    # Calculate the number of voxels needed in each direction
    num_x = math.ceil((max_x - min_x) / voxel_size)
    num_y = math.ceil((max_y - min_y) / voxel_size)
    num_z = math.ceil((max_z - min_z) / voxel_size)

    # Create voxels
    for i in range(num_x):
        for j in range(num_y):
            for k in range(num_z):
                # Calculate the position of each voxel
                x = min_x + i * voxel_size + voxel_size / 2
                y = min_y + j * voxel_size + voxel_size / 2
                z = min_z + k * voxel_size + voxel_size / 2
                voxel_location = (x, y, z)

                # Create an empty cube at the voxel position
                obj = create_mesh(f"voxel_{i}{j}{k}", 
                                  voxel_location, 
                                  mesh="cube", 
                                  size=voxel_size, 
                                  subdivisions=1)
                bpy.context.scene.collection.objects.link(obj)
                
                # Add concentrations vector
                concentration_vector = np.random.rand(5)
                for idx, conc in enumerate(concentration_vector):
                    obj[f"conc_{idx}"] = float(conc)

                # Store the voxel information
                voxel = Voxel(obj, voxel_location, concentration_vector)
                voxels.append(voxel)

    return voxels