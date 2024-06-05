from typing import Optional
from typing_extensions import override

import numpy as np
import math
from scipy.ndimage import convolve
from goo.utils import *


class Molecule(): 
    """A molecule involved in the reaction-diffusion system.

    Args:
        concentration (float): The initial concentration of the molecule.
        diffusion_rate (float): The diffusion rate of the molecule.

    Attributes:
        concentration (float): The concentration of the molecule.
        diffusion_rate (float): The diffusion rate of the molecule.
    """

    def __init__(self, name: str, conc: float, D: float, mat=None):
        self.conc = conc
        self.D = D
        self._mat = mat

    def __repr__(self):
        str = (f"Molecule(concentration={self.conc}"
               "diffusion_rate={self.diffusion_rate})")
        return str

    @property
    def name(self) -> str:
        """Name of the cell. Also defines the name of related forces and
        collections of effectors.
        """
        return self.name
        
    @name.setter
    def name(self, name: str):
        self.name = name


class Voxel(BlenderObject):

    def __init__(self, obj: bpy.types.Object, loc: tuple, conc: tuple, mat=None):
        super(Voxel, self).__init__(obj)
        self.location = loc
        self.conc = conc
        self._mat = mat
        if self._mat:
            self.obj.data.materials.append(mat)   

    def recolor(self, color: tuple[float, float, float]):
        """Recolors the material of the cell.

        This function changes the diffuse color of the cell's material to the
        specified color while preserving the alpha value. If the material uses
        nodes, it also updates the 'Base Color' input of any nodes that have it.

        Args:
            color: A tuple (r, g, b) representing the new color to apply.
        """
        r, g, b = color
        self._mat.diffuse_color = (r, g, b, 1)

        if self._mat.use_nodes:
            for node in self._mat.node_tree.nodes:
                if "Base Color" in node.inputs:
                    _, _, _, a = node.inputs["Base Color"].default_value
                    node.inputs["Base Color"].default_value = r, g, b, a


class ReactionDiffusionSystem:
    """A reaction-diffusion system simulation.

    Args:
        width (int): The width of the 2D grid.
        height (int): The height of the 2D grid.
        molecule_a (Molecule): The first type of molecule in the system.
        diffusion_rate_a (float): The diffusion rate of molecule A.

    Attributes:
        width (int): The width of the boundary.
        height (int): The height of the boundary.
        molecule_a (Molecule): The first type of molecule in the system.
        diffusion_rate_a (float): The diffusion rate of molecule A.
        grid_a (numpy.ndarray): The concentration grid. 
    """
    color = (0.07, 0.21, 0.3)

    def __init__(self, loc: tuple, size: tuple, voxel_size: int, mol_a: Molecule):
        self.size = size
        self.loc = loc
        self.voxel_size = voxel_size
        self.mol_a = mol_a
        self.grid = np.zeros((size[0], size[1], size[2]))
        self._voxel_data = []
        self._voxels = set()

    @property
    def voxels(self) -> list[Voxel]:
        """The list of cells associated with this cell type."""
        return list(self._voxels)
    
    def initialize(self, initial_concentration: Optional[float] = None):
        num_x, num_y, num_z = self.size

        min_x = self.loc[0] - self.size[0] / 2
        min_y = self.loc[1] - self.size[1] / 2
        min_z = self.loc[2] - self.size[2] / 2

        # Create voxels
        for i in range(num_x):
            for j in range(num_y):
                for k in range(num_z):
                    # Calculate the position of each voxel
                    x = min_x + i * self.voxel_size + self.voxel_size / 2
                    y = min_y + j * self.voxel_size + self.voxel_size / 2
                    z = min_z + k * self.voxel_size + self.voxel_size / 2
                    voxel_location = (x, y, z)

                    # linear gradient along x-axis
                    if initial_concentration is not None:
                        concentration = initial_concentration * (i / (num_x - 1))
                    else:
                        # Otherwise, use default concentration
                        concentration = self.default_concentration

                    # Store the voxel information with the calculated concentration
                    self._voxel_data.append((voxel_location, concentration))

    def update(self):
        # Define the 3D diffusion kernel
        diffusion_kernel = np.array([[[0, 0.125, 0],
                                    [0.125, 0.25, 0.125],
                                    [0, 0.125, 0]],

                                    [[0.125, 0.25, 0.125],
                                     [0.25, 1, 0.25],
                                     [0.125, 0.25, 0.125]],

                                    [[0, 0.125, 0],
                                     [0.125, 0.25, 0.125],
                                     [0, 0.125, 0]]])

        # Perform diffusion using convolution
        da = self.mol_a.D * convolve(self.grid, diffusion_kernel, mode='wrap')

        # Update the concentration grid
        self.grid += da

        # Normalize values to a reasonable range
        max_val = np.max(self.grid)
        if max_val > 0:
            self.grid /= max_val

    def is_converged(self, threshold=0.1):
        """Check if the concentration grid has converged."""
        change = np.max(np.abs(self.grid - np.mean(self.grid)))
        return change < threshold
    
    def toggle_voxel_grid(self):
        """Create a 3D grid in Blender.

        Args:
            loc: Center of the grid.
            size: Dimensions of the grid (width, height, depth).
            voxel_size: Size of a voxel in the grid.
        """

        for idx, loc in enumerate(self._voxel_data):      
            # Create an empty cube at the voxel position
            name = f"voxel_{idx}"
            obj = create_mesh(name, 
                              loc[0], 
                              mesh="cube", 
                              size=self.voxel_size, 
                              subdivisions=1)
            bpy.context.scene.collection.objects.link(obj)  

            # Add concentrations vector
            concentration_vector = np.random.rand(5)
            for idx, conc in enumerate(concentration_vector):
                obj[f"conc_{idx}"] = float(conc)

            color = self.__class__.color

            mat = create_material(f"{name}_material", color=color) if color else None
            # Store the voxel information
            voxel = Voxel(obj=obj, loc=loc, conc=concentration_vector, mat=mat)
            self._voxels.add(voxel)