from importlib import reload
import goo
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()


celltype = goo.OpaqueType("cells", physics_on=True)
cell = celltype.create_cell("cell", (1, 1, 0), mesh_kwargs={"size": 2})
cell.motion_force.loc = (0, 0, 0)
cell.motion_force.strength = 10000

sim = goo.Simulator(celltypes=[celltype])
sim.setup_world()
sim.add_handlers([AdhesionLocationHandler()])

for i in range(1, 10):
    bpy.context.scene.frame_set(i)
    print(i, cell.obj_eval.location, cell.COM())
