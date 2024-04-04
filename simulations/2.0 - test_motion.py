from importlib import reload
import bpy
import goo
import goo.handler as handler

reload(goo)

goo.reset_modules()
goo.reset_scene()

celltype = goo.create_celltype("default", physics_on=True)
cell = celltype.create_cell("cell", (1, 1, 0), size=2)
cell.motion_force.loc = (0, 0, 0)
cell.motion_force.strength = 10000

sim = goo.Simulator(celltypes=[celltype])
origin_handler = handler.ForceUpdateHandler()
sim.add_handler(origin_handler)

sim.toggle_gravity(False)
sim.run_simulation(start=1, end=50)

for i in range(1, 10):
    bpy.context.scene.frame_set(i)
    print(i, cell.obj_eval.location, cell.get_COM())
