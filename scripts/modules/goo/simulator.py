import bpy

from goo.handler import TimingHandler
from datetime import datetime


class Simulator:
    def __init__(self, celltypes=None, physics_dt=1):
        self.handlers = []
        self.celltypes = celltypes
        self.physics_dt = physics_dt

    def get_cells(self):
        return [cell for celltype in self.celltypes for cell in celltype.cells]

    def add_handler(self, handler):
        handler.setup(self.get_cells, self.physics_dt)
        self.handlers.append(handler.run)

    def run_simulation(self, start=1, end=250):
        start_time = datetime.now()
        bpy.context.scene.frame_set(start)
        bpy.context.scene.frame_start = start
        bpy.context.scene.frame_end = end
        bpy.app.handlers.frame_change_post.extend(self.handlers)

    def toggle_gravity(self, on):
        bpy.context.scene.use_gravity = on
