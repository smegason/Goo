from typing import Callable, Union
from typing_extensions import override, Optional
from abc import ABC, abstractmethod

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
from goo.gene import Gene
from goo.molecule import Molecule, DiffusionSystem


class Handler(ABC):
    def setup(
        self,
        get_cells: Callable[[], list[Cell]],
        get_diffsystem: Callable[[], DiffusionSystem],
        dt: float,
    ) -> None:
        """Set up the handler.

        Args:
            get_cells: A function that, when called,
                retrieves the list of cells that may divide.
            dt: The time step for the simulation.
        """
        self.get_cells = get_cells
        self.get_diffsystem = get_diffsystem
        self.dt = dt

    @abstractmethod
    def run(self, scene: bpy.types.Scene, depsgraph: bpy.types.Depsgraph) -> None:
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

    @override
    def setup(
        self,
        get_cells: Callable[[], list[Cell]],
        get_diffsystems: Callable[[], list[DiffusionSystem]],
        dt,
    ):
        """Build the KD-Tree from the grid coordinates if not already built."""
        super(DiffusionHandler, self).setup(get_cells, get_diffsystems, dt)
        self.get_diffsystem().build_kdtree()

    def run(self, scene, depsgraph) -> None:
        self.get_diffsystem().simulate_diffusion()


class NetworkHandler(Handler):
    """Handler for gene regulatory networks."""

    def run(self, scene, despgraph):
        for cell in self.get_cells():
            cell.step_grn(self.get_diffsystem())


class RecenterHandler(Handler):
    """Handler for updating cell origin and location of
    cell-associated adhesion locations every frame."""

    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            cell.recenter()

            cell_size = cell.major_axis().length() / 2
            for force in cell.adhesion_forces:
                if not force.enabled():
                    continue
                force.min_dist = cell_size - 0.4
                force.max_dist = cell_size + 0.4

            if cell.motion_force:
                cell.move()


class GrowthPIDHandler(Handler):
    @override
    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            cell.step_growth()


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
            cell.move(dir)


"""Possible properties by which cells are colored."""
Colorizer = Enum("Colorizer", ["PRESSURE", "VOLUME", "RANDOM", "GENE"])


class ColorizeHandler(Handler):
    """Handler for coloring cells based off of a specified property.

    Cells are colored on a blue-red spectrum, based on the relative value
    of the specified property to all other cells. For example, the cell with
    the highest pressure is colored red, while the cell with an average
    pressure is colored purple.

    Attributes:
        colorizer (Colorizer): the property by which cells are colored.
        gene (str): optional, the gene off of which cell color is based.
    """

    def __init__(
        self,
        colorizer: Colorizer = Colorizer.PRESSURE,
        gene: Union[Gene, str] = None,
        range: Optional[tuple] = None,
    ):
        self.colorizer = colorizer
        self.gene = gene
        self.range = range

    def _scale(self, values):
        if self.range is None:
            # Scaled based on min and max
            return (values - np.min(values)) / max(np.max(values) - np.min(values), 1)
        if self.range is not None:
            # Truncate numbers into specific range, then scale based on max of range.
            min, max = self.range
            values = np.minimum(values, max)
            values = np.maximum(values, min)
            return (values - min) / (max - min)

    def run(self, scene, depsgraph):
        match self.colorizer:
            case Colorizer.PRESSURE:
                ps = np.array([cell.pressure for cell in self.get_cells()])
                values = self._scale(ps)
            case Colorizer.VOLUME:
                vs = np.array([cell.volume() for cell in self.get_cells()])
                values = self._scale(vs)
            case Colorizer.GENE:
                gs = np.array(
                    [cell.metabolites[self.gene] for cell in self.get_cells()]
                )
                values = self._scale(gs)
            case Colorizer.RANDOM:
                values = np.random.rand(len(self.get_cells()))
            case _:
                raise ValueError(
                    "Colorizer must be one of PRESSURE, VOLUME, GENE, or RANDOM."
                )

        red = Vector((1.0, 0.0, 0.0))
        blue = Vector((0.0, 0.0, 1.0))

        for cell, p in zip(self.get_cells(), values):
            color = blue.lerp(red, p)
            cell.recolor(tuple(color))


def _get_divisions(cells: list[Cell]) -> list[tuple[str, str, str]]:
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
def _contact_area(
    cell1: Cell, cell2: Cell, threshold=0.1
) -> tuple[float, float, float, float]:
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
def _contact_areas(cells: list[Cell], threshold=4) -> tuple[dict, dict]:
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


@staticmethod
def _shape_features(cells: list[Cell]) -> tuple[float, float, float, float]:
    """Calculate a set of shape features of a cell.

    Inlcudes the aspect ratio, sphericity

    Args:
        cell: A cell.

    Returns:
        Shape features (aspect ratio, sphericity, compactness, sav_ratio).
    """

    aspect_ratios = []
    sphericities = []
    compactnesses = []
    sav_ratios = []

    for cell in cells:
        aspect_ratio = cell.aspect_ratio()
        sphericity = cell.sphericity()
        compactness = cell.compactness()
        sav_ratio = cell.sav_ratio()

        aspect_ratios.append(aspect_ratio)
        sphericities.append(sphericity)
        compactnesses.append(compactness)
        sav_ratios.append(sav_ratio)

    return (aspect_ratios, sphericities, compactnesses, sav_ratios)


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
    SHAPE_FEATURES = auto()
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
    def setup(
        self,
        get_cells: Callable[[], list[Cell]],
        get_diffsystems: Callable[[], list[DiffusionSystem]],
        dt,
    ):
        super(DataExporter, self).setup(get_cells, get_diffsystems, dt)
        self.time_start = datetime.now()
        out = {"seed": bpy.context.scene["seed"], "frames": []}

        if self.path:
            with open(self.path, "w") as f:
                f.write(json.dumps(out))
        else:
            print(out)
        self.run(bpy.context.scene, bpy.context.evaluated_depsgraph_get())

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

        if self.options & DataFlag.SHAPE_FEATURES:
            aspect_ratios, sphericities, compactnesses, sav_ratios = _shape_features(
                self.get_cells()
            )
            frame_out["aspect_ratios"] = aspect_ratios
            frame_out["aspect_ratios"] = aspect_ratios
            frame_out["compactnesses"] = compactnesses
            frame_out["sav_ratios"] = sav_ratios

        if self.options & DataFlag.CONTACT_AREAS:
            areas, ratios = _contact_areas(self.get_cells())
            frame_out["contact_areas"] = areas
            frame_out["contact_ratios"] = ratios

        if self.options & DataFlag.GRID:
            for diff_system in self.get_diffsystem():
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
