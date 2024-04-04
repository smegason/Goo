from typing import Callable

import bpy, bmesh
from goo.cell import Cell


# TODO: this is probably not necessary, handlers can exist within motion, divider submodules
class Handler:
    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        self.get_cells = get_cells
        self.dt = dt

    def run(self, scene, depsgraph):
        raise NotImplementedError("Subclasses must implement run() method.")


# TODO: handler to update mesh origins to COM
class DivisionHandler(Handler):
    pass


class MotionHandler(Handler):
    pass


class RemeshHandler(Handler):
    def __init__(self, freq):
        self.freq = freq

    def run(self, scene, depsgraph):
        if scene.frame_current % self.freq != 0:
            return
        for cell in self.get_cells():
            if cell.physics_enabled:
                bm = bmesh.new()
                bm.from_mesh(cell.obj_eval.to_mesh())

                cell.disable_physics()

                bm.to_mesh(cell.obj.data)
                bm.free()
                cell.remesh()
                cell.recenter()

                cell.enable_physics()
                cell.cloth_mod.point_cache.frame_start = scene.frame_current


class ForceUpdateHandler(Handler):
    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            cell_size = cell.get_major_axis().length() / 2
            for force in cell.adhesion_forces:
                if not force.enabled():
                    continue
                force.loc = cell.get_COM()
                force.min_dist = cell_size - 0.4
                force.max_dist = cell_size + 0.4


class TimingHandler(Handler):
    def __init__(self, start_time):
        self.start_time = start_time

    def run(self):
        pass
