# under development

import bpy
from goo import goo
from importlib import reload
reload(goo)

# clear stuff from before
bpy.app.handlers.frame_change_pre.clear()
bpy.app.handlers.frame_change_post.clear()
cell_objs = [obj for obj in bpy.data.objects if obj.name.startswith("Cell_")]
for cell_obj in cell_objs:
    bpy.data.objects.remove(cell_obj)

# keep frame handler from crashing interface
bpy.types.RenderSettings.use_lock_interface = True

goo.setup_world()

# create first cell
cell = goo.Cell(name_string="Cell_", loc=(0, 0, 0))
goo.make_cell(cell)
obj = cell.get_blender_object()

# obj.modifiers["Cloth"].settings.shrink_min = 0
# obj.modifiers["Cloth"].settings.keyframe_insert(data_path="shrink_min", frame=1)
# obj.modifiers["Cloth"].settings.shrink_min = -1
# obj.modifiers["Cloth"].settings.keyframe_insert(data_path="shrink_min", frame=40)
handler = goo.handler_class()
handler.set_division_rate('sphere', 30)
handler.set_growth_rate('sphere', 0.5)
handler.active_cell_types = ['sphere']

print(handler.division_rates)
bpy.app.handlers.frame_change_post.append(handler.growth_handler)
bpy.app.handlers.frame_change_post.append(handler.div_handler)