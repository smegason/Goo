from typing import Optional
from typing_extensions import override

import bpy
import numpy as np
import math
from goo.utils import *


class Boundary(BlenderObject):
    """A boundary for cells."""

    def __init__(self, obj: bpy.types.Object):
        super(Boundary, self).__init__(obj)

    def setup_physics(self):
        """Set up physics for the boundary."""
        BoundaryCollisionConstructor().construct(self.obj)

    def create_voxels(self, voxel_size):
        """Instantiate grid inside the boundary."""
        bounding_box = self.obj.bound_box
        min_x = min([v[0] for v in bounding_box])
        min_y = min([v[1] for v in bounding_box])
        min_z = min([v[2] for v in bounding_box])
        max_x = max([v[0] for v in bounding_box])
        max_y = max([v[1] for v in bounding_box])
        max_z = max([v[2] for v in bounding_box])
        
        voxels = []

        # Calculate the number of voxels needed in each direction
        num_x = math.ceil((max_x - min_x) / voxel_size)
        num_y = math.ceil((max_y - min_y) / voxel_size)
        num_z = math.ceil((max_z - min_z) / voxel_size)
        
        # Create empty cubes as voxels
        for i in range(num_x):
            for j in range(num_y):
                for k in range(num_z):
                    # Calculate the position of each voxel
                    x = min_x + i * voxel_size
                    y = min_y + j * voxel_size
                    z = min_z + k * voxel_size
                    location = (x, y, z)

                    # Create an empty cube at the voxel position
                    bpy.ops.object.empty_add(type='CUBE', 
                                             location=location, 
                                             radius=voxel_size / 2)
                    voxel = bpy.context.object
                    voxels.append(voxel.location)

                    # Add custom properties
                    concentration_vector = np.random.rand(5)
                    for idx, conc in enumerate(concentration_vector):
                        voxel[f"conc_{idx}"] = float(conc)

    @property
    def size(self) -> int:
        """Size of the boundary."""
        return self.obj.size
    
    @size.setter
    def size(self, size: int):
        self.obj.size = size


def create_boundary(loc: tuple, size: float, mesh: str = "icosphere"):
    """Create a boundary.

    Args:
        loc: Center of the boundary.
        size: Radius of the boundary.
        mesh: Shape fo the boundary.
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


def voxelize_boundary(boundary, voxel_size):
    """Voxelize a boundary.

    Args:
        boundary: A Boundary object.
        voxel_size: Size of a voxel in the boundary grid.
    """
    voxels = boundary.create_voxels(voxel_size)
    return voxels 