from typing import Callable
from enum import Enum
import time

import numpy as np
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
        target_volume=30,
    ):
        self.growth_type = growth_type
        self.growth_rate = growth_rate  # in cubic microns per frame
        self.Kp = 0.05
        self.Ki = 0.000001
        self.Kd = 0.5
        self.PID_scale = 60
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

        initial_volume = cell.volume()
        initial_pressure = cell.pressure

        cell["previous_pressure"] = initial_pressure
        cell["volume"] = initial_volume
        cell["next_volume"] = initial_volume
        cell["target_volume"] = self.target_volume

    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            if not cell.physics_enabled:
                continue
            if "target_volume" not in cell:
                self.initialize_PID(cell)
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
            volume_deviation = 1 - cell.volume() / cell["next_volume"]

            # print(
            #     f"Target volume: {cell['target_volume']}; "
            #     f"Volume: {cell.volume()}; "
            #     f"Next volume: {cell['next_volume']}; "
            #     f"Volume deviation: {volume_deviation}; "
            # )

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

            # print(f"New pressure for {cell.name}: {cell.pressure}")


ForceDist = Enum("ForceDist", ["UNIFORM", "GAUSSIAN"])


class MotionHandler(Handler):
    def __init__(
        self, distribution: ForceDist = ForceDist.UNIFORM, distribution_size=0.5
    ):
        self.distribution = distribution
        self.distribution_size = distribution_size

    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        super(MotionHandler, self).setup(get_cells, dt)
        for cell in self.get_cells():
            self.initialize_motion(cell)

    def initialize_motion(self, cell: Cell):
        if not cell.motion_force.enabled:
            cell.motion_force.enable()
        cell.motion_force.loc = cell.COM()
        cell["current_position"] = cell.COM()
        cell["previous_position"] = cell.COM()

    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            if not cell.physics_enabled:
                continue
            if "current_position" not in cell:
                self.initalize_motion(cell)
            cell["current_position"] = cell.COM()
            disp = Vector(cell["current_position"]) - Vector(cell["previous_position"])
            dist = disp.length

            match self.distribution:
                case ForceDist.UNIFORM:
                    noise = Vector(
                        np.random.uniform(
                            low=-self.distribution_size,
                            high=self.distribution_size,
                            size=(3,),
                        )
                    )
                case ForceDist.GAUSSIAN:
                    noise = Vector(
                        np.random.normal(loc=0, scale=self.distribution_size, size=(3,))
                    )
                case _:
                    raise ValueError(
                        "Motion noise distribution must be one of UNIFORM or GAUSSIAN."
                    )
            cell.motion_force.loc = cell.COM() + noise


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


class TimingHandler(Handler):
    def __init__(self, start_time):
        self.start_time = start_time

    def run(self):
        pass
