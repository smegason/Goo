from importlib import reload
import goo
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()


celltype = goo.CellType("A")
# goo.create_force("Force", (0, 0, 0), 30000)
celltype.set_motion(60000)
cell = celltype.create_cell("cellA", (1, 1, 0))
cell.motion_force.loc = (0, 0, 0)

sim = goo.Simulator([celltype])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        # RemeshHandler(freq=10, voxel_size=0.25),
        AdhesionLocationHandler(),
        MotionHandler(),
    ]
)
