from importlib import reload
import goo
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()


celltype = goo.CellType("cellsA", pattern="standard")
celltype.motion_strength = 200
cell = celltype.create_cell("cell", (0, 0, 0), size=2)
cell.stiffness = 15
cell.pressure = 1

sim = goo.Simulator([celltype])
sim.setup_world()
sim.add_handlers(
    [
        RecenterHandler(),
        # RemeshHandler(),
    ]
)

for i in range(1, 21):
    bpy.context.scene.frame_set(i)
    cell.move(Vector((1, 0, 0)))
    print(cell.COM()[0])
