from typing import Optional
from typing_extensions import override
import numpy as np
import math

from scipy.ndimage import laplace
from scipy.integrate import solve_ivp
from scipy.spatial import KDTree

from goo.utils import *
from goo.cell import Cell


class Molecule:
    """A molecule involved in the diffusion system. 

    Args:
        name (str): The name of the molecule.
        conc (float): The initial concentration of the molecule.
        D (float): The diffusion rate of the molecule.
    """
    def __init__(
        self,
        name: str,
        conc: float,
        D: float,
    ):
        self._name = name
        self._D = D
        self._conc = conc
        
    def __repr__(self):
        str = (f"Molecule(concentration={self.conc}"
               f"diffusion_rate={self.D})")
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

    @property
    def conc(self) -> float:
        """The concentration of the molecule."""
        return self.conc
    
    @conc.setter
    def conc(self, conc: float):
        self.conc = conc

    @property
    def D(self) -> float:
        """The diffusion rate of the molecule."""
        return self.D
    
    @D.setter
    def D(self, D: float):
        self.D = D

        '''if conc == 'random':
            self.concentration = np.random.rand(*grid_size)
        elif conc == 'linear_gradient':
            self.concentration = np.linspace(0, 1, grid_size[0]).reshape(-1, 1, 1)
            self.concentration = np.tile(self.concentration, (1, grid_size[1], grid_size[2]))
        elif conc == 'localized' and location is not None:
            self.concentration = np.zeros(grid_size)
            self.concentration[location] = 1.0'''


class DiffusionSystem:
    def __init__(
        self, 
        molecules: list[Molecule], 
        cells: list[Cell] = None,
        grid_size: tuple = (25, 25, 25), 
        grid_center: tuple = (0, 0, 0),
        time_step: float = 0.1, 
        total_time: int = 10
    ) -> None: 
        self._molecules = molecules
        self._grid_size = grid_size
        self._time_step = time_step
        self._total_time = total_time
        self._cells = cells
        self._grid_center = grid_center
        self._grid_concentrations = None
        self._kd_tree = None
        self._laplacian_kernel = np.array([
            [[1,  1,  1], [1,  1,  1], [1,  1,  1]],
            [[1,  1,  1], [1, -26, 1], [1,  1,  1]],
            [[1,  1,  1], [1,  1, 1], [1,  1,  1]]
        ]) / 26.0
        
    def _build_kdtree(self):
        """Initialize the grid and build the KD-Tree."""
        x = np.linspace(-(self._grid_size[0] - 1) / 2, (self._grid_size[0] - 1) / 2, self._grid_size[0])
        y = np.linspace(-(self._grid_size[1] - 1) / 2, (self._grid_size[1] - 1) / 2, self._grid_size[1])
        z = np.linspace(-(self._grid_size[2] - 1) / 2, (self._grid_size[2] - 1) / 2, self._grid_size[2])
        grid_points = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)
        grid_points += np.array(self._grid_center)
        self._kd_tree = KDTree(grid_points)
        self._grid_concentrations = np.zeros(self._grid_size)
        print("self._grid_concentrations.shape", self._grid_concentrations.shape)
        return self._kd_tree
    
    def update_concentration(self, index, value):
        """Update the concentration value at the given grid index."""
        # Convert flat index to 3D index
        z_index, y_index, x_index = np.unravel_index(index, self._grid_size)
        self._grid_concentrations[z_index, y_index, x_index] += value
        return self._grid_concentrations[z_index, y_index, x_index]

    @property
    def molecules(self) -> list[Molecule]:
        """The list of molecules in the system."""
        return self.molecules
    
    @molecules.setter
    def molecules(self, molecules: list[Molecule]):
        self.molecules = molecules

    @property
    def grid_size(self) -> tuple:
        """The size of the 3D grid."""
        return self.grid_size
    
    @grid_size.setter
    def grid_size(self, grid_size: tuple):
        self.grid_size = grid_size

    @property
    def time_step(self) -> float:
        """The time step of the simulation."""
        return self.time_step
    
    @time_step.setter
    def time_step(self, time_step: float):
        self.time_step = time_step

    @property
    def total_time(self) -> int:
        """The total time of the simulation."""
        return self.total_time
    
    @total_time.setter
    def total_time(self, total_time: int):
        self.total_time = total_time

    def diffusion_system(self, t, y):
        concentrations = []
        offset = 0
        for i, molecule in enumerate(self.molecules):
            C = y[offset:offset + self.grid_size_prod].reshape(self.grid_size)
            dC_dt = self.diffusion_rates[i] * laplace(C)
            concentrations.append(dC_dt.ravel())
            offset += self.grid_size_prod
        return np.concatenate(concentrations)

    def run_simulation(self):
        # Initialize the combined state for all molecules
        y0 = np.concatenate([molecule.concentration.flatten() for molecule in self.molecules])
        t_span = (0, self.total_time)
        t_eval = np.arange(0, self.total_time, self.time_step)
        
        # Solve the diffusion system
        solution = solve_ivp(self.diffusion_system, t_span, y0, method='RK45', t_eval=t_eval)

        # Reshape the solution for each molecule
        results = []
        offset = 0
        for molecule in self.molecules:
            concentration_sol = solution.y[offset:offset + np.prod(self.grid_size)].reshape(self.grid_size + (len(t_eval),))
            results.append(concentration_sol)
            offset += np.prod(self.grid_size)
        return results