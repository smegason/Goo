import sys
import os
import numpy as np
import bpy

from goo.handler import Handler
from goo.cell import CellType


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

        # Set units to the metric system
        bpy.context.scene.unit_settings.system = "METRIC"
        bpy.context.scene.unit_settings.scale_length = 0.0001
        bpy.context.scene.unit_settings.system_rotation = "DEGREES"
        bpy.context.scene.unit_settings.length_unit = "MICROMETERS"
        bpy.context.scene.unit_settings.mass_unit = "MILLIGRAMS"
        bpy.context.scene.unit_settings.time_unit = "SECONDS"
        bpy.context.scene.unit_settings.temperature_unit = "CELSIUS"

        # Turn off gravity
        self.toggle_gravity(False)

        # Set up rendering environment
        node_tree = bpy.context.scene.world.node_tree
        tree_nodes = node_tree.nodes
        tree_nodes.clear()

        # Add background node
        node_background = tree_nodes.new(type="ShaderNodeBackground")
        node_environment = tree_nodes.new("ShaderNodeTexEnvironment")
        scripts_paths = bpy.utils.script_paths()
        try:
            node_environment.image = bpy.data.images.load(
                scripts_paths[-1] + "/modules/goo/missile_launch_facility_01_4k.hdr"
            )
        except Exception:
            print(sys.exc_info())
            print(
                """WARNING FROM GOO: To enable proper rendering you must have
                /modules/goo/missile_launch_facility_01_4k.hdr in the right location"""
            )
        node_environment.location = -300, 0

        # Add output node
        node_output = tree_nodes.new(type="ShaderNodeOutputWorld")
        node_output.location = 200, 0
        # Link all nodes
        links = node_tree.links
        links.new(node_environment.outputs["Color"], node_background.inputs["Color"])
        links.new(node_background.outputs["Background"], node_output.inputs["Surface"])

        # # set film to transparent to hide background
        bpy.context.scene.render.film_transparent = True

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

    def get_cells_func(self, celltypes=None):
        celltypes = celltypes if celltypes is not None else self.celltypes

        def get_cells():
            return [cell for celltype in celltypes for cell in celltype.cells]

        return get_cells

    def add_handler(self, handler: Handler, celltypes: list[CellType] = None):
        handler.setup(self.get_cells_func(celltypes), self.physics_dt)
        bpy.app.handlers.frame_change_post.append(handler.run)

    def add_handlers(self, handlers: list[Handler], celltypes: list[CellType] = None):
        for handler in handlers:
            self.add_handler(handler, celltypes)

    def render(self, start=1, end=250, save=True, path=None, camera=False):
        bpy.context.scene.frame_start = start
        bpy.context.scene.frame_end = end

        bpy.context.scene.render.image_settings.file_format = "PNG"
        bpy.context.scene.render.ffmpeg.format = "MPEG4"

        if not path:
            path = os.path.dirname(bpy.context.scene.render.filepath)
        elif path and not save:
            print("Save path set but render will not be saved!")

        print("----- SIMULATION START -----")
        for i in range(1, end + 1):
            bpy.context.scene.frame_set(i)
            bpy.context.scene.render.filepath = os.path.join(path, f"{i:04d}")
            if camera:
                bpy.ops.render.render(write_still=save)
            else:
                bpy.ops.render.opengl(write_still=save)
        bpy.context.scene.render.filepath = path
        print("\n----- SIMULATION END -----")

    def run(self, end=250):
        print("----- SIMULATION START -----")
        for i in range(1, end + 1):
            print(i, end=" ")
            bpy.context.scene.frame_set(i)
        print("\n----- SIMULATION END -----")
