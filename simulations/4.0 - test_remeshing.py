from importlib import reload
import goo
from goo.division import *
from goo.handler import *

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.create_celltype("A")
force = goo.create_force("Force", (0, 0, 0), strength=-500)
celltype.homo_adhesions.add_force(force)

cell = celltype.create_cell("cellA", (0, 0, 0))
# logic = BisectDivisionLogic(margin=0.025)
# cells = cell.divide(logic)
# logic.flush()
# for cell in cells:
#     cell.loc = cell.loc * 2
#     cell.enable_physics(forces=False)
#     for force in cell.adhesion_forces:
#         force.disable()

sim = goo.Simulator(celltypes=[celltype])
sim.toggle_gravity(False)


# def remesh_handler(scene, depsgraph):
#     for cell in cells:
#         print(f"{scene.frame_current}\tRemeshing: {cell.name}")
#         if not cell.physics_enabled:
#             continue
#         print(f"{scene.frame_current}\tRemeshing: {cell.name}")

#         # Update mesh and disable physics
#         bm = bmesh.new()
#         bm.from_mesh(cell.obj_eval.to_mesh())
#         cell.disable_physics()
#         bm.to_mesh(cell.obj.data)
#         bm.free()

#         # Perform remeshing operations
#         cell.remesh()

#         # Recenter and re-enable physics
#         cell.recenter()
#         cell.enable_physics(forces=False)
#         cell.cloth_mod.point_cache.frame_start = scene.frame_current


# bpy.app.handlers.frame_change_post.append(remesh_handler)

sim.run_simulation(start=1, end=250)


def store_settings(mod: bpy.types.bpy_struct):
    settings = {}
    for p in mod.bl_rna.properties:
        id = p.identifier
        if not p.is_readonly:
            settings[id] = getattr(mod, id)
        elif id in ["settings", "collision_settings", "effector_weights"]:
            settings[id] = store_settings(getattr(mod, id))
    return settings


def declare_settings(mod: bpy.types.bpy_struct, settings: dict):
    for id, setting in settings.items():
        if isinstance(setting, dict):
            declare_settings(getattr(mod, id), settings[id])
        else:
            setattr(mod, id, setting)


settings = store_settings(cell.cloth_mod)
cell.obj.modifiers.remove(cell.cloth_mod)
cell.obj.modifiers.new(name="Cloth", type="CLOTH")
declare_settings(cell.cloth_mod, settings)
print(settings)
