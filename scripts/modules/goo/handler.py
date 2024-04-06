from typing import Callable
from enum import Enum
import time

import bpy, bmesh
from goo.cell import Cell


# TODO: this is probably not necessary, handlers can exist within motion, divider submodules
class Handler:
    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        self.get_cells = get_cells
        self.dt = dt

    def run(self, scene, depsgraph):
        raise NotImplementedError("Subclasses must implement run() method.")


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
            cell.remesh(self.voxel_size)

            # Recenter and re-enable physics
            cell.recenter()
            cell.enable_physics()
            cell.cloth_mod.point_cache.frame_start = scene.frame_current

    def cast_to_sphere(self, cell, factor):
        with bpy.context.temp_override(active_object=cell.obj):
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
        self.target_volume = target_volume

    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        super(GrowthPIDHandler, self).setup(get_cells, dt)
        for cell in self.get_cells():
            self.initialize_PID(cell)

    def initialize_PID(self, cell: Cell):
        cell["Kp"] = self.Kp
        cell["Ki"] = self.Ki
        cell["Kd"] = self.Kd
        cell["PID_scale"] = 60
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
        time = scene.frame_current * self.dt
        for cell in self.get_cells():
            if not cell.physics_enabled:
                continue
            if "target_volume" not in cell:
                self.initialize_PID(cell)
            cell["volume"] = cell.volume()

            if self.growth_type == Growth.LINEAR:
                cell["next_volume"] += cell["growth_rate"] * self.dt
            elif self.growth_type == Growth.EXPONENTIAL:
                cell["next_volume"] *= 1 + cell["growth_rate"] * self.dt
            elif self.growth_type == Growth.LOGISTIC:
                cell["next_volume"] = cell["next_volume"] * (
                    1
                    + cell["growth_rate"]
                    * (1 - cell["next_volume"] / cell["target_volume"])
                    * self.dt
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


class MotionHandler(Handler):
    pass


class TimingHandler(Handler):
    def __init__(self, start_time):
        self.start_time = start_time

    def run(self):
        pass
