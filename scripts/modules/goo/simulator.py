import os
from datetime import datetime
import numpy as np
import bpy

from goo.handler import Handler


class Simulator:
    def __init__(self, celltypes=[], physics_dt=1):
        self.celltypes = celltypes
        self.physics_dt = physics_dt
        self.addons = ["add_mesh_extra_objects"]

    def setup_world(self, seed=1):
        # Enable addons
        for addon in self.addons:
            self.enable_addon(addon)

        # Set random seed
        np.random.seed(seed)
        bpy.context.scene["seed"] = seed

        # Turn off gravity
        bpy.context.scene.use_gravity = False

        # Set units to the metric system
        bpy.context.scene.unit_settings.system = "METRIC"
        bpy.context.scene.unit_settings.scale_length = 0.0001
        bpy.context.scene.unit_settings.system_rotation = "DEGREES"
        bpy.context.scene.unit_settings.length_unit = "MICROMETERS"
        bpy.context.scene.unit_settings.mass_unit = "MILLIGRAMS"
        bpy.context.scene.unit_settings.time_unit = "SECONDS"
        bpy.context.scene.unit_settings.temperature_unit = "CELSIUS"

    def enable_addon(self, addon):
        if addon not in bpy.context.preferences.addons:
            bpy.ops.preferences.addon_enable(module=addon)
            print(f"Addon '{addon}' has been enabled.")
        else:
            print(f"Addon '{addon}' is already enabled.")

    def toggle_gravity(self, on):
        bpy.context.scene.use_gravity = on

    def add_celltype(self, celltype):
        self.celltypes.append(celltype)

    def add_celltypes(self, celltypes):
        self.celltypes.extend(celltypes)

    def get_cells(self):
        return [cell for celltype in self.celltypes for cell in celltype.cells]

    def add_handler(self, handler: Handler):
        handler.setup(self.get_cells, self.physics_dt)
        bpy.app.handlers.frame_change_post.append(handler.run)

    def add_handlers(self, handlers: list[Handler]):
        for handler in handlers:
            self.add_handler(handler)

    def run_simulation(self, start=1, end=250):
        start_time = datetime.now()
        bpy.context.scene.frame_set(start)
        bpy.context.scene.frame_start = start
        bpy.context.scene.frame_end = end

    def render(self, start=1, end=250, save=True, path=None, camera=False):
        bpy.context.scene.frame_start = start
        bpy.context.scene.frame_end = end

        bpy.context.scene.render.image_settings.file_format = "PNG"
        bpy.context.scene.render.ffmpeg.format = "MPEG4"

        if not path:
            path = os.path.dirname(bpy.context.scene.render.filepath)
        elif path and not save:
            print("Save path set but render will not be saved!")

        for i in range(1, end + 1):
            bpy.context.scene.frame_set(i)
            bpy.context.scene.render.filepath = os.path.join(path, f"{i:04d}")
            if camera:
                bpy.ops.render.render(write_still=save)
            else:
                bpy.ops.render.opengl(write_still=save)
        bpy.context.scene.render.filepath = path
