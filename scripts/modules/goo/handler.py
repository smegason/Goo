from typing import Callable, Any, List
from typing_extensions import override

from enum import Enum, Flag, auto
from datetime import datetime
import json
import os

import numpy as np
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.ndimage import laplace

import bpy
import bmesh
from mathutils import Vector
from goo.cell import Cell
from goo.molecule import DiffusionSystem


class Handler:
    def setup(
        self, 
        get_cells: Callable[[], list[Cell]], 
        get_diffsystems: Callable[[], list[DiffusionSystem]],
        dt: float
    ):
        """Set up the handler.

        Args:
            get_cells: A function that, when called, 
                retrieves the list of cells that may divide.
            dt: The time step for the simulation.
        """
        self.get_cells = get_cells
        self.get_diff_systems = get_diffsystems
        self.dt = dt

    def run(self, scene: bpy.types.Scene, depsgraph: bpy.types.Depsgraph):
        """Run the handler.

        This is the function that gets passed to Blender, to be called
        upon specified events (e.g. post-frame change).

        Args:
            scene: The Blender scene.
            depsgraph: The dependency graph.
        """
        raise NotImplementedError("Subclasses must implement run() method.")


class RemeshHandler(Handler):
    """Handler for remeshing cells at given frequencies.

    Attributes:
        freq (int): Number of frames between remeshes.
        smooth_factor (float): Factor to pass to `bmesh.ops.smooth_vert`. 
            Disabled if set to 0.
        voxel_size (float): Factor to pass to `voxel_remesh()`. Disabled if set to 0.
        sphere_factor (float): Factor to pass to Cast to sphere modifier. 
            Disabled if set to 0.
    """

    def __init__(self, freq=1, voxel_size=None, smooth_factor=0.1, sphere_factor=0):
        self.freq = freq
        self.voxel_size = voxel_size
        self.smooth_factor = smooth_factor
        self.sphere_factor = sphere_factor

    @override
    def run(self, scene, depsgraph):
        if scene.frame_current % self.freq != 0:
            return
        for cell in self.get_cells():
            if not cell.physics_enabled:
                continue

            # Update mesh and disable physics
            bm = bmesh.new()
            bm.from_mesh(cell.obj_eval.to_mesh())
            cell.disable_physics()
            if self.smooth_factor:
                bmesh.ops.smooth_vert(
                    bm,
                    verts=bm.verts,
                    factor=self.smooth_factor,
                )
            bm.to_mesh(cell.obj.data)
            bm.free()
            cell.recenter()

            if self.voxel_size is not None:
                cell.remesh(self.voxel_size)
                cell.recenter()
            else: 
                cell.remesh()
                cell.recenter()
            
            # Recenter and re-enable physics
            cell.enable_physics()
            cell.cloth_mod.point_cache.frame_start = scene.frame_current


class DiffusionHandler(Handler):
    """Handler for simulating diffusion of a substance in the grid in the scene.

    Args:
        diffusionSystem: The reaction-diffusion system to simulate.
    """

    def __init__(self, diffusionSystem: DiffusionSystem) -> None:
        self.diffusionSystem = diffusionSystem
        self.kd_tree = None

    def build_kd_tree(self):
        """Build the KD-Tree from the grid coordinates if not already built."""
        self.kd_tree = self.diffusionSystem._build_kdtree()

    def read_molecular_signal(
        self,
        cell: Cell,
        cell_distances: np.ndarray,
        indices: np.ndarray,
        radius: float
    ) -> None: 
        """Read the concentration of molecules in the cell."""
        for mol_idx, molecule in enumerate(self.diffusionSystem.molecules):
            valid_indices = ~np.isinf(cell_distances) & (cell_distances >= radius)
            print(f"Valid indices: {valid_indices}")    
            if valid_indices.any():
                index = indices[valid_indices][0]
                total_conc = self.diffusionSystem.get_concentration(mol_idx, index)
                print(f"Conc of cell {cell.name} for {molecule._name}: {total_conc}")
                cell.molecules_conc.update({molecule._name: total_conc})
        print(f"Molecular concentrations: {cell.molecules_conc}")

    def update_molecular_signal(
        self,
        cell: Cell,
        cell_distances: np.ndarray,
        indices: np.ndarray,
        radius: float
    ) -> None: 
        """Update the concentration of molecules in the cell."""
                        
        k = 0.1
        for mol_idx, molecule in enumerate(self.diffusionSystem._molecules):
            valid_indices = ~np.isinf(cell_distances) & (cell_distances >= radius)
            valid_distances = cell_distances[valid_indices]
            valid_indices = indices[valid_indices]

            for cell_distance, index in zip(valid_distances, valid_indices):
                add_conc = k * (cell_distance / radius)
                self.diffusionSystem.update_concentration(mol_idx, index, add_conc)

    def diffuse(self, mol_idx: int):
        """Update the concentration of molecules based on diffusion."""
        conc = self.diffusionSystem._grid_concentrations[mol_idx]
        laplacian = laplace(conc, mode='wrap')
        diff_coeff = self.diffusionSystem._molecules[mol_idx]._D
        conc += self.diffusionSystem._time_step * diff_coeff * laplacian
        conc = np.clip(conc, 0, None)
        self.diffusionSystem._grid_concentrations[mol_idx] = conc

    def simulate_diffusion(self):
        """Run the diffusion simulation over the total time."""
        tot_time = self.diffusionSystem._total_time
        t_step = self.diffusionSystem._time_step
        num_steps = int(tot_time / t_step)
        for _ in range(num_steps):
            for mol_idx in range(len(self.diffusionSystem._molecules)):
                self.diffuse(mol_idx)

    @override
    def run(self, scene, depsgraph) -> None:
        if self.kd_tree is None:
            self.build_kd_tree()
        
        print("Current frame", scene.frame_current)

        # diffuse molecules on grid
        # self.simulate_diffusion()

        cells = self.get_cells()
        
        for cell in cells:
            radius = cell.get_radius()
            com = cell.COM()
            scaling_factor = 1 / self.diffusionSystem._element_size[0]
            cell_distances, indices = self.kd_tree.query(com, 
                                                         k=1000 * scaling_factor**2, 
                                                         distance_upper_bound=1.25*radius, 
                                                         p=2)

            if len(cell_distances) > 0 and not np.all(np.isinf(cell_distances)):
                # self.update_molecular_signal(cell, cell_distances, indices, radius)
                self.read_molecular_signal(cell, cell_distances, indices, radius)
            else:
                # If no valid distances, set concentration to 0
                for molecule in self.diffusionSystem._molecules:
                    cell.molecules_conc.update({molecule._name: 0})
                    print(cell.molecules_conc)
                    print(f"Total conc of cell {cell.name} for {molecule._name}: 0")


class AdhesionLocationHandler(Handler):
    """Handler for updating cell-associated adhesion locations every frame."""

    @override
    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            cell_size = cell.major_axis().length() / 2

            for force in cell.adhesion_forces:
                if not force.enabled():
                    continue
                force.loc = cell.COM()
                force.min_dist = cell_size - 0.4
                force.max_dist = cell_size + 0.4


"""Possible types of growth."""
Growth = Enum("Growth", ["LINEAR", "EXPONENTIAL", "LOGISTIC"])


class GrowthPIDHandler(Handler):
    """Handler for simulating cell growth based off of internal pressure.

    Growth is determined by a PID controller, in which changes to a cell's
    internal pressure governs how much it grows in the next frame.

    Attributes:
        growth_type (Growth): Type of growth exhibited by cells.
        growth_rate (float): Rate of growth of cells.
        initial_pressure (float): Initial pressure of cells.
        target_volume (float): Target volume of cells.
        Kp (float): P variable of the PID controller.
        Ki (float): I variable of the PID controller.
        Kd (float): D variable of the PID controller.
    """

    def __init__(
        self,
        growth_type: Growth = Growth.LINEAR,
        growth_rate: float = 1,
        initial_pressure=0.01,
        target_volume=30,
        Kp=0.05,
        Ki=0.00001,
        Kd=0.5,
    ):
        self.growth_type = growth_type
        self.growth_rate = growth_rate  # in cubic microns per frame
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.PID_scale = 60
        self.initial_pressure = initial_pressure
        self.target_volume = target_volume

    @override
    def setup(self, 
              get_cells: Callable[[], list[Cell]], 
              get_diffsystems: Callable[[], list[DiffusionSystem]], 
              dt):
        super(GrowthPIDHandler, self).setup(get_cells, get_diffsystems, dt)
        for cell in self.get_cells():
            self.initialize_PID(cell)

    def initialize_PID(self, cell: Cell):
        """Initialize PID controller for a cell.

        Args:
            cell: Cell to initialize PID controller.
        """
        cell["Kp"] = self.Kp
        cell["Ki"] = self.Ki
        cell["Kd"] = self.Kd
        cell["PID_scale"] = self.PID_scale
        cell["growth_rate"] = self.growth_rate

        cell["integral"] = 0
        cell["previous_error"] = 0

        cell["previous_pressure"] = self.initial_pressure
        cell["next_volume"] = cell.volume()
        cell["target_volume"] = self.target_volume

    @override
    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            if "target_volume" not in cell:
                self.initialize_PID(cell)
            if "divided" in cell and cell["divided"]:
                # if divided, reset certain values
                cell["previous_pressure"] = self.initial_pressure
                cell["next_volume"] = cell.volume()
            if not cell.physics_enabled:
                continue

            cell["volume"] = cell.volume()

            match self.growth_type:
                case Growth.LINEAR:
                    cell["next_volume"] += cell["growth_rate"] * self.dt
                case Growth.EXPONENTIAL:
                    cell["next_volume"] *= 1 + cell["growth_rate"] * self.dt
                case Growth.LOGISTIC:
                    cell["next_volume"] = cell["next_volume"] * (
                        1
                        + cell["growth_rate"]
                        * (1 - cell["next_volume"] / cell["target_volume"])
                        * self.dt
                    )
                case _:
                    raise ValueError(
                        "Growth type must be one of LINEAR, EXPONENTIAL, or LOGISTIC."
                    )
            cell["next_volume"] = min(cell["next_volume"], cell["target_volume"])
            volume_deviation = 1 - cell["volume"] / cell["next_volume"]

            # Update pressure based on PID output
            error = volume_deviation
            integral = cell["integral"] + error
            derivative = error - cell["previous_error"]
            pid = cell["Kp"] * error + cell["Ki"] * integral + cell["Kd"] * derivative

            cell.pressure = cell["previous_pressure"] + pid * cell["PID_scale"]

            # Update previous error and pressure for the next iteration
            cell["previous_error"] = error
            cell["integral"] = integral
            cell["previous_pressure"] = cell.pressure


"""Possible distributions of random motion."""
ForceDist = Enum("ForceDist", ["CONSTANT", "UNIFORM", "GAUSSIAN"])


class RandomMotionHandler(Handler):
    """Handler for simulating random cell motion.

    At every frame, the direction of motion is randomized, and the strength
    of the motion force is randomly selected from a specified distribution.

    Attributes:
        distribution (ForceDist): Distribution of random strength of motion force.
        max_strength (int): Maximum strength motion force.
    """

    def __init__(
        self,
        distribution: ForceDist = ForceDist.UNIFORM,
        max_strength: int = 0,
    ):
        self.distribution = distribution
        self.max_strength = max_strength

    @override
    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            if not cell.physics_enabled:
                continue
            if not cell.motion_force.enabled:
                cell.motion_force.enable()

            dir = Vector(np.random.uniform(low=-1, high=1, size=(3,)))
            match self.distribution:
                case ForceDist.CONSTANT:
                    strength = self.max_strength
                case ForceDist.UNIFORM:
                    strength = np.random.random_sample() * self.max_strength
                case ForceDist.GAUSSIAN:
                    strength = np.random.normal() * self.max_strength
                case _:
                    raise ValueError(
                        "Motion noise distribution must be one of UNIFORM or GAUSSIAN."
                    )
            cell.motion_force.strength = strength
            cell.move_towards(dir)


"""Possible properties by which cells are colored."""
Colorizer = Enum("Colorizer", ["PRESSURE", "VOLUME", "RANDOM"])


class ColorizeHandler(Handler):
    """Handler for coloring cells based off of a specified property.

    Cells are colored on a blue-red spectrum, based on the relative value
    of the specified property to all other cells. For example, the cell with
    the highest pressure is colored red, while the cell with an average
    pressure is colored purple.

    Attributes:
        colorizer (Colorizer): the property by which cells are colored.
    """

    def __init__(self, colorizer: Colorizer = Colorizer.PRESSURE):
        self.colorizer = colorizer

    @override
    def run(self, scene, depsgraph):
        red = Vector((1.0, 0.0, 0.0))
        blue = Vector((0.0, 0.0, 1.0))

        match self.colorizer:
            case Colorizer.PRESSURE:
                ps = np.array([cell.pressure for cell in self.get_cells()])
                ps = (ps - np.min(ps)) / max(np.max(ps) - np.min(ps), 1)
            case Colorizer.VOLUME:
                ps = np.array([cell.volume() for cell in self.get_cells()])
                ps = (ps - np.min(ps)) / max(np.max(ps) - np.min(ps), 1)
            case Colorizer.RANDOM:
                ps = np.random.rand(len(self.get_cells()))
            case _:
                raise ValueError(
                    "Colorizer must be one of PRESSURE, VOLUME, or RANDOM."
                )

        for cell, p in zip(self.get_cells(), ps):
            color = blue.lerp(red, p)
            cell.recolor(tuple(color))


@staticmethod
def _get_divisions(cells: list[Cell]):
    """Calculate a list of cells that have divided in the past frame.

    Each element of the list contains a tuple of three names: that of the mother
    cell, and then the two daughter cells.

    Returns:
        List of tuples of mother and daughter cell names.
    """
    divisions = set()
    for cell in cells:
        if "divided" in cell and cell["divided"]:
            divisions.add(
                (cell.name[:-2], cell.name[:-2] + ".0", cell.name[:-2] + ".1")
            )
    return list(divisions)


@staticmethod
def _contact_area(cell1: Cell, cell2: Cell, threshold=0.1):
    """Calculate the contact areas between two cells.

    Contact is defined as two faces that are within a set threshold distance
    from each other.

    Args:
        cell1: First cell to calculate contact.
        cell2: Second cell to calculate contact.
        threshold: Maximum distance between two faces of either cell to consider
            as contact.

    Returns:
        A tuple containing for elements:
            - Total area of cell1 in contact with cell2
            - Total area of cell2 in contact with cell1
            - Ratio of area of cell1 in contact with cell2
            - Ratio of area of cell2 in contact with cell1
    """
    faces1 = cell1.obj_eval.data.polygons
    faces2 = cell2.obj_eval.data.polygons

    centers1 = [cell1.obj_eval.matrix_world @ f.center for f in faces1]
    centers2 = [cell2.obj_eval.matrix_world @ f.center for f in faces2]

    dists = np.array(cdist(centers1, centers2, "euclidean"))

    contact_faces1 = np.any(dists < threshold, axis=1)
    contact_faces2 = np.any(dists < threshold, axis=0)

    areas1 = np.array([f.area for f in faces1])
    areas2 = np.array([f.area for f in faces2])

    contact_areas1 = np.sum(areas1[contact_faces1])
    contact_areas2 = np.sum(areas2[contact_faces2])

    ratio1 = contact_areas1 / np.sum(areas1)
    ratio2 = contact_areas2 / np.sum(areas2)

    return contact_areas1, contact_areas2, ratio1, ratio2


@staticmethod
def _contact_areas(cells, threshold=4):
    """Calculate the pairwise contact areas between a list of cells.

    Contact is calculated heuristically by first screening cells that are within
    a certain threshold distance between each other.

    Args:
        cells: The list of cells to calculate contact areas over.
        threshold: The maximum distance between cells to consider them for contact.

    Returns:
        A list of tuples containing pairwise contact areas and contact ratios.
            See :func:`_contact_area`.
    """
    coms = [cell.COM() for cell in cells]
    dists = squareform(pdist(coms, "euclidean"))

    mask = dists < threshold
    mask = np.triu(mask, k=1)

    pairs = np.where(mask)

    areas = {cell.name: [] for cell in cells}
    ratios = {cell.name: [] for cell in cells}
    for i, j in zip(pairs[0], pairs[1]):
        contact_area_i, contact_area_j, ratio_i, ratio_j = _contact_area(
            cells[i], cells[j]
        )
        areas[cells[i].name].append((cells[j].name, contact_area_i))
        areas[cells[j].name].append((cells[i].name, contact_area_j))
        ratios[cells[i].name].append((cells[j].name, ratio_i))
        ratios[cells[j].name].append((cells[i].name, ratio_j))

    return areas, ratios


class _all:
    def __get__(self, instance, cls):
        return ~cls(0)


class DataFlag(Flag):
    """Enum of data flags used by the :func:`DataExporter` handler.

    Attributes:
        TIMES: time elapsed since beginning of simulation.
        DIVISIONS: list of cells that have divided and their daughter cells.
        MOTION_PATH: list of the current position of each cell.
        FORCE_PATH: list of the current positions of the associated 
            motion force of each cell.
        VOLUMES: list of the current volumes of each cell.
        PRESSURES: list of the current pressures of each cell.
        CONTACT_AREAS: list of contact areas between each pair of cells.
        CONCENTRATIONS: concentrations of each molecule in the grid system.
    """

    TIMES = auto()
    DIVISIONS = auto()
    MOTION_PATH = auto()
    FORCE_PATH = auto()
    VOLUMES = auto()
    PRESSURES = auto()
    CONTACT_AREAS = auto()
    GRID = auto()
    CELL_CONCENTRATIONS = auto()

    ALL = _all()


class DataExporter(Handler):
    """Handler for the reporting and saving of data generated during
    the simulation.

    Attributes:
        path (str): Path to save .json file of calculated metrics. If empty, statistics
            are printed instead.
        options: (DataFlag): Flags of which metrics to calculated and save/print. 
            Flags can be combined by binary OR operation, 
            i.e. `DataFlag.TIMES | DataFlag.DIVISIONS`.
    """

    def __init__(self, path="", options: DataFlag = DataFlag.ALL):
        self.path = path
        self.options = options

    @override
    def setup(self, 
              get_cells: Callable[[], list[Cell]], 
              get_diffsystems: Callable[[], list[DiffusionSystem]],
              dt):
        super(DataExporter, self).setup(get_cells, get_diffsystems, dt)
        self.time_start = datetime.now()
        out = {"seed": bpy.context.scene["seed"], "frames": []}

        if self.path:
            with open(self.path, "w") as f:
                f.write(json.dumps(out))
        else:
            print(out)
        self.run(bpy.context.scene, bpy.context.evaluated_depsgraph_get())

    @override
    def run(self, scene, depsgraph):
        frame_out = {"frame": scene.frame_current}

        if self.options & DataFlag.TIMES:
            frame_out["time"] = (datetime.now() - self.time_start).total_seconds()
        if self.options & DataFlag.DIVISIONS:
            frame_out["divisions"] = _get_divisions(self.get_cells())

        frame_out["cells"] = []
        for cell in self.get_cells():
            cell_out = {"name": cell.name}
            frame_out["cells"].append(cell_out)

            if self.options & DataFlag.MOTION_PATH:
                cell_out["loc"] = tuple(cell.loc)
            if self.options & DataFlag.FORCE_PATH:
                cell_out["motion_loc"] = tuple(cell.motion_force.loc)
            if self.options & DataFlag.VOLUMES:
                cell_out["volume"] = cell.volume()
            if self.options & DataFlag.PRESSURES and cell.physics_enabled:
                cell_out["pressure"] = cell.pressure
            if self.options & DataFlag.CELL_CONCENTRATIONS:
                cell_out["concentrations"] = cell.molecules_conc

        if self.options & DataFlag.CONTACT_AREAS:
            areas, ratios = _contact_areas(self.get_cells())
            frame_out["contact_areas"] = areas
            frame_out["contact_ratios"] = ratios

        if self.options & DataFlag.GRID:
            for diff_system in self.get_diff_systems():
                grid_conc = diff_system._grid_concentrations
                mol = diff_system._molecules[0]
                if mol._name not in frame_out:
                    frame_out[mol._name] = {
                        "concentrations": self._convert_numpy_to_list(grid_conc)
                    }

        if self.path:
            with open(self.path, "r") as f:
                out = json.load(f)
                out["frames"].append(frame_out)
            with open(self.path, "w") as f:
                f.write(json.dumps(out))
        else:
            print(frame_out)

    def _convert_numpy_to_list(self, obj):
        """Convert numpy arrays to lists for JSON serialization."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: self._convert_numpy_to_list(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_numpy_to_list(i) for i in obj]
        else:
            return obj
        

class SliceExporter(Handler):
    """Handler to save Z slices of the simulation at each frame."""

    def __init__(
        self,
        output_dir: str = "", 
        z_range: tuple[float, float] = (5, -5), 
        z_step: float = 0.2, 
        microscope_dt: int = 10
    ):
        self.output_dir = output_dir
        self.z_range = z_range
        self.z_step = z_step
        self.microscope_dt = microscope_dt

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def run(self, scene: bpy.types.Scene, depsgraph: bpy.types.Depsgraph):
        # export slices every microscope_dt
        if scene.frame_current % self.microscope_dt != 0:
            return
        
        time_step = scene.frame_current // self.microscope_dt
        time_step_dir = os.path.join(self.output_dir, f"T{time_step:03d}")
        if not os.path.exists(time_step_dir):
            os.makedirs(time_step_dir)

        # Check if a camera named 'SliceExporterCamera' already exists
        camera = bpy.data.objects.get('SliceExporterCamera')

        if camera is None:
            # No camera exists, so create and center it
            camera = create_and_center_camera(location=(0, 0, 5), target=(0, 0, 0))
        
        # Set the created or existing camera as the active camera
        bpy.context.scene.camera = camera

        z_start, z_end = self.z_range
        num_slices = int(abs(z_end - z_start) / self.z_step) + 1
        print("Number of slices:", num_slices)

        # Iterate over z-axis
        for i in range(num_slices):
            # current slice position
            slice_z = z_start + i * self.z_step if z_start < z_end else z_start - i * self.z_step
            camera.location.z = slice_z
            print("Camera location:", camera.location)

            # thin slice
            camera.data.clip_start = 0  # front of the camera
            camera.data.clip_end = self.z_step  # end of the slice being photographed

            slice_filename = f"slice_z{i:03d}.png"
            filepath = os.path.join(time_step_dir, slice_filename)
            bpy.context.scene.render.filepath = filepath
            
            # OpenGL render from the active camera's perspective
            bpy.ops.render.opengl(write_still=True, view_context=False)

            # add a step to convert from RGB  to 8-bit grayscale


@staticmethod
def create_and_center_camera(location=(0, 0, 10), target=(0, 0, 0)):
    """Create a new camera object, set its parameters, 
    center it to the specified location, and orient it towards the target.

    Args:
        location: The location to place the camera (default is (0, 0, 10)).
        target: The point the camera should face (default is the origin in (0, 0, 0)).
    """
    # new camera data block
    cam_data = bpy.data.cameras.new(name="SliceExporterCamera")
    cam_object = bpy.data.objects.new("SliceExporterCamera", cam_data)
    bpy.context.collection.objects.link(cam_object)
    cam_object.location = location
    direction = Vector(target) - cam_object.location
    cam_object.rotation_euler = direction.to_track_quat('Z', 'Y').to_euler()
    # acquisition parameters
    cam_object.data.type = 'ORTHO'
    cam_object.data.ortho_scale = 20

    return cam_object