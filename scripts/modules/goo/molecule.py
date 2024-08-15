from typing import Optional, Tuple
from typing_extensions import override
import numpy as np
import math

from scipy.ndimage import laplace
from scipy.integrate import solve_ivp
from scipy.spatial import KDTree

from goo.utils import *


class Molecule:
    """A molecule involved in the diffusion system.

    Args:
        name (str): The name of the molecule.
        conc (float): The initial concentration of the molecule.
        D (float): The diffusion rate of the molecule.
        gradient (str, optional): The gradient of the molecule. Defaults to None.
    """

    def __init__(
        self, name: str, conc: float, D: float, gradient: Optional[str] = None
    ):
        self.name = name
        self.D = D
        self.conc = conc
        self.gradient = gradient

    def __str__(self):
        return self.name


class DiffusionSystem:
    """A diffusion system that simulates the diffusion of molecules in a 3D grid.
    It also handles the secretion and sensing of molecular signals by cells.

    Args:
        molecules (list[Molecule]): The list of molecules in the system.
        grid_size (Tuple[int, int, int], optional): The size of the 3D grid. Defaults to (25, 25, 25).
        grid_center (Tuple[int, int, int], optional): The center of the 3D grid. Defaults to (0, 0, 0).
        time_step (float, optional): The time step of the simulation. Defaults to 0.1.
        total_time (int, optional): The total time of the simulation. Defaults to 10.
        element_size (Tuple[float, float, float], optional): The size of each element in the grid. Defaults to (1.0, 1.0, 1.0).
    """

    def __init__(
        self,
        molecules: list[Molecule],
        grid_size: Tuple[int, int, int] = (25, 25, 25),
        grid_center: Tuple[int, int, int] = (0, 0, 0),
        time_step: float = 0.1,
        total_time: int = 10,
        element_size=(1.0, 1.0, 1.0),
    ) -> None:
        self._molecules = molecules
        self._grid_size = grid_size
        self._time_step = time_step
        self._total_time = total_time
        self._grid_center = grid_center
        self._grid_concentrations = None
        self._kd_tree = None
        self._element_size = element_size
        self._laplacian_kernel = (
            np.array(
                [
                    [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
                    [[1, 1, 1], [1, -26, 1], [1, 1, 1]],
                    [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
                ]
            )
            / 26.0
        )

    def _build_kdtree(self):
        """Initialize the grid and build its corresponding KD-Tree."""
        x = (
            np.linspace(
                -(self._grid_size[0] - 1) / 2,
                (self._grid_size[0] - 1) / 2,
                self._grid_size[0],
            )
            * self._element_size[0]
        )

        y = (
            np.linspace(
                -(self._grid_size[1] - 1) / 2,
                (self._grid_size[1] - 1) / 2,
                self._grid_size[1],
            )
            * self._element_size[1]
        )

        z = (
            np.linspace(
                -(self._grid_size[2] - 1) / 2,
                (self._grid_size[2] - 1) / 2,
                self._grid_size[2],
            )
            * self._element_size[2]
        )

        grid_points = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)
        grid_points += np.array(self._grid_center)
        self._kd_tree = KDTree(grid_points)

        # initialize the grid to store concentrations for each molecule
        self._grid_concentrations = np.zeros((len(self._molecules), *self._grid_size))

        for idx, mol in enumerate(self._molecules):
            match mol.gradient:
                case None:
                    continue
                case "random":
                    variation = 0.1  # 10% variation
                    noise = np.random.normal(
                        0, variation * mol.conc, size=self._grid_size
                    )
                    rand_conc = mol.conc + noise
                    # conc is always positive
                    rand_conc = np.clip(rand_conc, 0, None)
                    self._grid_concentrations[idx] = rand_conc
                case "center":
                    center_index = self._get_nearest_grid_index(self._grid_center)
                    self._grid_concentrations[idx][center_index] += mol.conc
                case "linear":
                    x_grad = np.linspace(0, mol.conc, self._grid_size[0]).reshape(
                        -1, 1, 1
                    )

                    self._grid_concentrations[idx] = np.tile(
                        x_grad, (1, self._grid_size[1], self._grid_size[2])
                    )

        # add 8 empty points at the corners of the grid in Blender
        corner_offsets = [
            (-1, -1, -1),
            (-1, -1, 1),
            (-1, 1, -1),
            (-1, 1, 1),
            (1, -1, -1),
            (1, -1, 1),
            (1, 1, -1),
            (1, 1, 1),
        ]
        for offset in corner_offsets:
            corner_position = (
                self._grid_center[0]
                + offset[0] * (self._grid_size[0] - 1) / 2 * self._element_size[0],
                self._grid_center[1]
                + offset[1] * (self._grid_size[1] - 1) / 2 * self._element_size[1],
                self._grid_center[2]
                + offset[2] * (self._grid_size[2] - 1) / 2 * self._element_size[2],
            )
            bpy.ops.object.empty_add(type="PLAIN_AXES", location=corner_position)

        return self._kd_tree

    def _get_nearest_grid_index(self, point):
        """Get the nearest grid index for the given point."""
        return self._kd_tree.query(point)[1]

    def update_concentration(self, mol_idx, index, value):
        """Update the concentration value of a molecule at a given voxel in the grid."""
        # Convert flat index to 3D index
        z_index, y_index, x_index = np.unravel_index(index, self._grid_size)
        self._grid_concentrations[mol_idx, z_index, y_index, x_index] += value
        return self._grid_concentrations[mol_idx, z_index, y_index, x_index]

    def get_concentration(self, mol_idx, index):
        """Get the concentration value of a molecule at a given voxel in the grid."""
        z_index, y_index, x_index = np.unravel_index(index, self._grid_size)
        return self._grid_concentrations[mol_idx, z_index, y_index, x_index]

    def get_all_concentrations(self, center, radius):
        """Get concentrations of all molecules in a sphere with given center and radius."""
        distances, indices = self._kd_tree.query(
            center, k=500, distance_upper_bound=1.5 * radius, p=2
        )
        signaling_concs = {}
        for mol_idx, molecule in enumerate(self._molecules):
            signaling_concs[molecule] = 0
            for cell_distance, index in zip(distances, indices):
                if np.isinf(cell_distance):
                    continue
                elif cell_distance >= radius:
                    signaling_concs[molecule] += self.get_concentration(mol_idx, index)
        return signaling_concs

    @property
    def molecules(self) -> list[Molecule]:
        """The list of molecules in the system."""
        return self._molecules

    @molecules.setter
    def molecules(self, molecules: list[Molecule]):
        self._molecules = molecules

    @property
    def gradient(self) -> list[(Molecule, str)]:
        """The gradient of the molecules in the system."""
        return self._gradient

    @gradient.setter
    def gradient(self, gradient: list[(Molecule, str)]):
        self._gradient = gradient

    @property
    def grid_size(self) -> tuple:
        """The size of the 3D grid."""
        return self._grid_size

    @grid_size.setter
    def grid_size(self, grid_size: tuple):
        self._grid_size = grid_size

    @property
    def time_step(self) -> float:
        """The time step of the simulation."""
        return self._time_step

    @time_step.setter
    def time_step(self, time_step: float):
        self._time_step = time_step

    @property
    def total_time(self) -> int:
        """The total time of the simulation."""
        return self._total_time

    @total_time.setter
    def total_time(self, total_time: int):
        self._total_time = total_time

    def diffuse(self, t, y):
        concentrations = []
        offset = 0
        for i, molecule in enumerate(self.molecules):
            C = y[offset : offset + self.grid_size_prod].reshape(self.grid_size)
            dC_dt = self.diffusion_rates[i] * laplace(C)
            concentrations.append(dC_dt.ravel())
            offset += self.grid_size_prod
        return np.concatenate(concentrations)

    def run(self):
        return
        # Initialize the combined state for all molecules
        y0 = np.concatenate([molecule.conc for molecule in self.molecules])
        t_span = (0, self.total_time)
        t_eval = np.arange(0, self.total_time, self.time_step)

        # Solve the diffusion system
        solution = solve_ivp(self.diffuse(), t_span, y0, method="RK45", t_eval=t_eval)

        # Reshape the solution for each molecule
        results = []
        offset = 0
        for molecule in self.molecules:
            concentration_sol = solution.y[offset : offset + np.prod(self.grid_size)]
            concentration_sol = concentration_sol.reshape(
                self.grid_size + (len(t_eval),)
            )
            results.append(concentration_sol)
            offset += np.prod(self.grid_size)
        return results
