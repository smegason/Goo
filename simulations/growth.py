# under development

import bpy
from Goo import Goo
from importlib import reload
reload(Goo)

# clear stuff from before
bpy.app.handlers.frame_change_pre.clear()
bpy.app.handlers.frame_change_post.clear()
cell_objs = [obj for obj in bpy.data.objects if obj.name.startswith("Cell_")]
for cell_obj in cell_objs:
    bpy.data.objects.remove(cell_obj)

# keep frame handler from crashing interface
bpy.types.RenderSettings.use_lock_interface = True

Goo.setup_world()

# create first cell
cell = Goo.Cell(name_string="Cell_", loc=(0, 0, 0))
Goo.make_cell(cell)
obj = cell.get_blender_object()

# obj.modifiers["Cloth"].settings.shrink_min = 0
# obj.modifiers["Cloth"].settings.keyframe_insert(data_path="shrink_min", frame=1)
# obj.modifiers["Cloth"].settings.shrink_min = -1
# obj.modifiers["Cloth"].settings.keyframe_insert(data_path="shrink_min", frame=40)

def growth_handler(scene, depsgraph): # WIP
    print("Frame:",scene.frame_current)
    cell_objs = [obj for obj in scene.objects if obj.name.startswith("Cell_")] # change between scene and depsgragh here

    num_cells = len(bpy.data.collections["Cells"].objects)
    for cell_obj in cell_objs:
        # cell_name = bpy.data.collections["Cells"].objects[i].name
        cell_name = cell_obj.name
        # cell_obj = bpy.data.objects[cell_name]
        # bpy.data.objects[cell_name].select_set(True)
        # bpy.context.view_layer.objects.active = bpy.data.objects[cell_name]
        # bpy.ops.object.modifier_apply(modifier="CLOTH")
        
        print(cell_obj.name," Volume:",Goo.calculate_volume(cell_obj)," Shrinking Factor:",cell_obj.modifiers["Cloth"].settings.shrink_min)
        cell_obj.modifiers["Cloth"].settings.shrink_min -= 0.01 # constantly changing shrink_min

bpy.app.handlers.frame_change_post.append(growth_handler)
bpy.app.handlers.frame_change_post.append(Goo.div_handler)
