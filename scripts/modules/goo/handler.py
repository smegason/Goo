from typing import Callable
from enum import Enum, Flag, auto
from datetime import datetime

import numpy as np
from scipy.spatial.distance import cdist, pdist, squareform
import bpy, bmesh
from mathutils import Vector
from goo.cell import Cell


class Handler:
    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        self.get_cells = get_cells
        self.dt = dt

    def run(self, scene, depsgraph):
        raise NotImplementedError("Subclasses must implement run() method.")


# TODO: remeshing seems to interfere with motion
class RemeshHandler(Handler):
    def __init__(self, freq=1, voxel_size=0.25, sphere_factor=0):
        self.freq = freq
        self.voxel_size = voxel_size
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
            bm.to_mesh(cell.obj.data)
            bm.free()

            # Perform remeshing operations
            if self.sphere_factor:
                self.cast_to_sphere(cell, self.sphere_factor)
            if self.voxel_size:
                cell.remesh(self.voxel_size)

            # Recenter and re-enable physics
            cell.recenter()
            cell.enable_physics()
            cell.cloth_mod.point_cache.frame_start = scene.frame_current

    def cast_to_sphere(self, cell, factor):
        with bpy.context.temp_override(active_object=cell.obj, object=cell.obj):
            cast_modifier = cell.obj.modifiers.new(name="Cast", type="CAST")
            cast_modifier.factor = factor
            bpy.ops.object.modifier_apply(modifier=cast_modifier.name)


class AdhesionLocationHandler(Handler):
    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            cell_size = cell.major_axis().length() / 2

            for force in cell.adhesion_forces:
                if not force.enabled():
                    continue
                force.loc = cell.COM()
                force.min_dist = cell_size - 0.4
                force.max_dist = cell_size + 0.4


Growth = Enum("Growth", ["LINEAR", "EXPONENTIAL", "LOGISTIC"])


class GrowthPIDHandler(Handler):
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

    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        super(GrowthPIDHandler, self).setup(get_cells, dt)
        for cell in self.get_cells():
            self.initialize_PID(cell)

    def initialize_PID(self, cell: Cell):
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


ForceDist = Enum("ForceDist", ["UNIFORM", "GAUSSIAN"])


class MotionHandler(Handler):
    def __init__(
        self,
        distribution: ForceDist = ForceDist.UNIFORM,
    ):
        self.distribution = distribution

    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        super(MotionHandler, self).setup(get_cells, dt)

    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            if not cell.physics_enabled:
                continue
            if not cell.motion_force.enabled:
                cell.motion_force.enable()

            match self.distribution:
                case ForceDist.UNIFORM:
                    dir = Vector(np.random.uniform(low=-1, high=1, size=(3,)))
                case ForceDist.GAUSSIAN:
                    raise NotImplementedError(
                        "Gaussian distribution not yet implemented"
                    )
                case _:
                    raise ValueError(
                        "Motion noise distribution must be one of UNIFORM or GAUSSIAN."
                    )
            cell.move_towards(dir)


Colorizer = Enum("Colorizer", ["PRESSURE", "VOLUME", "RANDOM"])


class ColorizeHandler(Handler):
    def __init__(self, colorizer: Colorizer = Colorizer.PRESSURE):
        self.colorizer = colorizer

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
            cell.recolor(color)


class SceneExtensionHandler(Handler):
    def __init__(self, end):
        self.end = end

    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            if cell.cloth_mod and cell.cloth_mod.point_cache.frame_end < self.end:
                cell.cloth_mod.point_cache.frame_end = self.end


def get_divisions(cells: list[Cell]):
    divisions = set()
    for cell in cells:
        if "divided" in cell and cell["divided"]:
            divisions.add(
                (cell.name[:-2], cell.name[:-2] + ".0", cell.name[:-2] + ".1")
            )
    return list(divisions)


def contact_area(cell1: Cell, cell2: Cell, threshold=0.1):
    faces1 = cell1.obj.data.polygons
    faces2 = cell2.obj.data.polygons

    centers1 = [cell1.obj.matrix_world @ f.center for f in faces1]
    centers2 = [cell2.obj.matrix_world @ f.center for f in faces2]

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


def contact_areas(cells, threshold):
    coms = [cell.COM() for cell in cells]
    dists = squareform(pdist(coms, "euclidean"))
    print(dists)

    mask = dists < threshold
    mask = np.triu(mask, k=1)

    pairs = np.where(mask)

    areas = {cell.name: [] for cell in cells}
    ratios = {cell.name: [] for cell in cells}
    for i, j in zip(pairs[0], pairs[1]):
        contact_area_i, contact_area_j, ratio_i, ratio_j = contact_area(
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
    TIMES = auto()
    DIVISIONS = auto()
    MOTION_PATH = auto()
    FORCE_PATH = auto()
    VOLUMES = auto()
    PRESSURES = auto()
    CONTACT_AREAS = auto()

    ALL = _all()


class DataExporter(Handler):
    def __init__(self, path, options: DataFlag):
        self.path = path
        self.options = options

    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        super(DataExporter, self).setup(get_cells, dt)
        self.time_start = datetime.now()
        self.out = {"seed": bpy.context.scene["seed"], "frames": []}

    def run(self, scene, depsgraph):
        frame_out = {}
        self.out["frames"].append(frame_out)

        if self.options & DataFlag.TIMES:
            frame_out["time"] = datetime.now() - self.time_start
        if self.options & DataFlag.DIVISIONS:
            frame_out["divisions"] = get_divisions(self.get_cells())

        frame_out["cells"] = []
        for cell in self.get_cells():
            cell_out = {"name": cell.name}
            frame_out["cells"].append(cell_out)

            if self.options & DataFlag.MOTION_PATH:
                cell_out["loc"] = cell.loc
            if self.options & DataFlag.FORCE_PATH:
                cell_out["motion_loc"] = cell.motion_force.loc
            if self.options & DataFlag.VOLUMES:
                cell_out["volume"] = cell.volume()
            if self.options & DataFlag.PRESSURES:
                cell_out["pressure"] = cell.pressure

        if self.options & DataFlag.CONTACT_AREAS:
            areas, ratios = contact_areas(self.get_cells())
            frame_out["contact_areas"] = areas
            frame_out["contact_ratios"] = ratios

        print(self.out)
